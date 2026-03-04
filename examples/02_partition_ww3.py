"""
STEP 2: Partition WW3 Spectra

This script processes 2D WW3 model spectra, applies the partitioning algorithm,
and saves the results to CSV.

Workflow:
1. Read coordinates/timestamps from auxiliary file
2. For each point, load WW3 spectrum at closest time or multiple times
3. Apply watershed partitioning algorithm
4. Save Hs, Tp, Dp for each identified partition

Data Type Support (automatically detected):
- SAR/NDBC with timestamp column: Processes closest WW3 time to specified timestamp
  - CSV must have: 'ref'+'sar_time' (SAR) or 'station_id'+'ndbc_time/obs_time' (NDBC)
  - One output file per row in CSV
  - Format: ww3_XXX_YYYYMMDD-HHMMSS.csv
  
- NDBC without timestamp: Processes multiple WW3 times based on TIME_INTERVAL_HOURS
  - CSV must have: 'station_id' (without time column)
  - Multiple output files per station
  - Format: ww3_STATION_YYYYMMDD-HHMMSS.csv
  
The script automatically detects the data type and processing mode based on CSV columns.

Usage:
------
python 02_partition_ww3.py --config config.yaml

Arguments:
  --config: Path to configuration YAML file (default: config.yaml)
"""

import os
import argparse
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from pathlib import Path

# Import from wasp package
from wasp.io_ww3 import find_closest_time, load_ww3_spectrum
from wasp.wave_params import calculate_wave_parameters
from wasp.partition import partition_spectrum
from wasp.utils import load_config

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='WW3 Spectral Partitioning',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--config',
        type=str,
        default='config.yaml',
        help='Path to configuration YAML file (default: config.yaml)'
    )
    
    return parser.parse_args()

# ============================================================================
# LOAD CONFIGURATION
# ============================================================================

args = parse_arguments()
CONFIG = load_config(args.config)

# ============================================================================
# CONFIGURATION
# ============================================================================

#case = 'lee'
# case = 'freddy'
# case = 'surigae'
case = 'all'

# ============================================================================
# DATA SOURCE SELECTION
# ============================================================================
# Choose: 'sar' or 'ndbc'
DATA_SOURCE = 'sar'  # <-- Change this to switch between SAR and NDBC

# ============================================================================
# PATHS AND PARAMETERS (from config.yaml)
# ============================================================================

# Partitioning parameters (from config.yaml)
MIN_ENERGY_THRESHOLD_FRACTION = CONFIG['partitioning']['ww3']['min_energy_fraction']
MAX_PARTITIONS = CONFIG['partitioning']['ww3']['max_partitions']
THRESHOLD_PERCENTILE = CONFIG['partitioning']['ww3']['threshold_percentile']
MERGE_FACTOR = CONFIG['partitioning']['ww3']['merge_factor']

# Set paths based on data source
if DATA_SOURCE == 'ndbc':
    CSV_PATH = f'../auxdata/ndbc_ww3_matches_{case}.csv'
    WW3_DATA_PATH = CONFIG['paths']['ww3_ndbc']
    OUTPUT_DIR = f'../data/{case}/partition-ww3-ndbc-{THRESHOLD_PERCENTILE}-{MERGE_FACTOR}'
elif DATA_SOURCE == 'sar':
    CSV_PATH = f'../auxdata/sar_matches_{case}_track.csv'
    WW3_DATA_PATH = CONFIG['paths']['ww3_sar']
    OUTPUT_DIR = f'../data/{case}/partition-ww3-{THRESHOLD_PERCENTILE}-{MERGE_FACTOR}'
elif DATA_SOURCE == 'cfosat':
    # Use the new CFOSAT-WW3 matching file generated from points.list
    CSV_PATH = f'../auxdata/cfosat_ww3_matches_complete_for_ww3.csv'
    WW3_DATA_PATH = CONFIG['paths']['ww3_cfosat']
    OUTPUT_DIR = f'../data/{case}/partition-ww3-cfosat-{THRESHOLD_PERCENTILE}-{MERGE_FACTOR}'
else:
    raise ValueError(f"Unknown DATA_SOURCE: {DATA_SOURCE}. Use 'sar' or 'ndbc'.")

# Time sampling for NDBC (from config.yaml)
TIME_INTERVAL_HOURS = CONFIG['processing']['time_interval_hours']

# Figure generation flag
GENERATE_FIGURES = CONFIG['processing'].get('generate_figures', True)

# Plotting parameters (from config.yaml)
PLOT_VMIN = CONFIG['plotting'].get('spectrum_vmin', 0.5)
PLOT_VMAX = CONFIG['plotting'].get('spectrum_vmax', 60.0)
PLOT_STEP = CONFIG['plotting'].get('spectrum_step', 0.5)
PLOT_PERIOD_MAX = CONFIG['plotting'].get('period_max', 25)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def detect_data_type(df):
    """
    Detect if CSV is for SAR, NDBC, or CFOSAT based on column names
    
    Returns:
    --------
    str: 'sar', 'ndbc', or 'cfosat'
    """
    if 'station_id' in df.columns:
        return 'ndbc'
    elif 'sar_ref' in df.columns and 'cfosat_filename' in df.columns:
        return 'cfosat'
    elif 'ref' in df.columns:
        return 'sar'
    else:
        raise ValueError("Cannot determine data type. Expected 'station_id' (NDBC), 'ref' (SAR), or 'sar_ref' + 'cfosat_filename' (CFOSAT).")


def count_significant_partitions(results, min_energy_threshold):
    """Count partitions with energy above of threshold"""
    count = 0
    for i in range(1, len(results['Hs'])):
        if results['energy'][i] > min_energy_threshold:
            count += 1
    return count


def plot_partition_spectrum(E2d_partition, freq, dirs_deg, title='',
                            hs=None, tp=None, dp=None, selected_time=None,
                            partition_num=None, energy_fraction=None,
                            partition_results=None, min_energy_threshold=None):
    """
    Plot 2D directional spectrum in polar coordinates with partition information.
    
    Parameters:
    -----------
    E2d_partition : ndarray
        2D directional spectrum [m²·s·rad⁻¹] for the partition
    freq : ndarray
        Frequency array [Hz]
    dirs_deg : ndarray
        Direction array [degrees]
    title : str, optional
        Plot title
    hs : float, optional
        Significant wave height [m]
    tp : float, optional
        Peak period [s]
    dp : float, optional
        Peak direction [degrees]
    selected_time : datetime, optional
        Timestamp of data
    partition_num : int, optional
        Partition number
    energy_fraction : float, optional
        Energy fraction relative to total (%)
    partition_results : dict, optional
        Results from partitioning (to show partition info on total spectrum)
    min_energy_threshold : float, optional
        Minimum energy threshold for significant partitions
    
    Returns:
    --------
    fig : matplotlib.figure.Figure
        Generated figure
    ax : matplotlib.axes.Axes
        Polar axes
    """
    Eplot = np.nan_to_num(E2d_partition, nan=0.0, neginf=0.0, posinf=0.0)

    # Ensure 1D arrays
    freq_plot = np.asarray(freq).flatten()
    dirs_plot = np.asarray(dirs_deg).flatten()

    # Convert directions to radians and sort
    dirs_rad_plot = np.radians(dirs_plot)
    sort_idx = np.argsort(dirs_rad_plot)
    dirs_sorted = dirs_rad_plot[sort_idx]
    Eplot_sorted = Eplot[:, sort_idx]

    # Ensure periodic continuity (0 to 2π)
    if not np.isclose(dirs_sorted[0], 0.0):
        dirs_sorted = np.insert(dirs_sorted, 0, 0.0)
        Eplot_sorted = np.insert(Eplot_sorted, 0, Eplot_sorted[:, 0], axis=1)
    if not np.isclose(dirs_sorted[-1], 2*np.pi):
        dirs_sorted = np.append(dirs_sorted, 2*np.pi)
        Eplot_sorted = np.concatenate([Eplot_sorted, Eplot_sorted[:, 0:1]], axis=1)

    # Radial axis = period (s)
    with np.errstate(divide='ignore', invalid='ignore'):
        period = np.where(freq_plot > 0, 1.0 / freq_plot, 0)

    theta, r = np.meshgrid(dirs_sorted, period)

    # Line contours from 0 to 60 m²·s·rad⁻¹
    vmin = 0.0
    vmax = 60.0
    step = 0.5
    levels = np.arange(vmin + step, vmax + step*0.51, step)

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='polar')

    # Use contour for line contours
    cs = ax.contour(theta, r, Eplot_sorted, levels, cmap='rainbow', vmin=vmin, vmax=vmax)

    # Axis styling
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rticks([5, 10, 15, 20])
    ax.set_yticklabels(['5s', '10s', '15s', '20s'], color='gray', fontsize=7.5)
    ax.set_rlim(0, PLOT_PERIOD_MAX)
    ax.set_rlabel_position(30)
    ax.tick_params(axis='y', colors='gray', labelsize=16)
    ticks = ['N','NE','E','SE','S','SW','W','NW']
    tick_angles = np.deg2rad(np.linspace(0, 315, 8))
    ax.set_xticks(tick_angles)
    ax.set_xticklabels(ticks)
    ax.tick_params(axis='x', colors='k', labelsize=16)
    ax.set_title(title, fontsize=16, color='k', pad=30)

    # Statistics box
    show_stats = selected_time is not None or hs is not None or tp is not None or dp is not None
    show_partitions = partition_results is not None and min_energy_threshold is not None
    
    if show_stats:
        # Adjust box size if showing partitions
        if show_partitions:
            stats_ax = fig.add_axes([0.75, 0.50, 0.20, 0.38], facecolor='white')
        else:
            stats_ax = fig.add_axes([0.75, 0.70, 0.20, 0.18], facecolor='white')
        
        stats_ax.patch.set_alpha(0.8)
        stats_ax.patch.set_edgecolor('black')
        stats_ax.patch.set_linewidth(1.5)
        stats_ax.axis('off')
        
        y_pos = 0.95
        
        # Title
        if partition_num is not None:
            title_text = f'Partition {partition_num}'
        else:
            title_text = 'Total Spectrum'
        stats_ax.text(0.5, y_pos, title_text, fontsize=13, color='k', 
                     ha='center', va='top', weight='bold', transform=stats_ax.transAxes)
        y_pos -= 0.08
        
        # Date
        if selected_time is not None:
            date_str = selected_time.strftime('%Y-%m-%d %H:%M:%S')
            stats_ax.text(0.5, y_pos, f'Date: {date_str}', fontsize=10, color='k', 
                         ha='center', va='top', transform=stats_ax.transAxes)
            y_pos -= 0.06
        
        # Wave parameters
        if hs is not None:
            stats_ax.text(0.5, y_pos, f'Hs: {hs:.2f} m', fontsize=11, color='k', 
                         ha='center', va='top', transform=stats_ax.transAxes)
            y_pos -= 0.06
        
        if tp is not None:
            stats_ax.text(0.5, y_pos, f'Tp: {tp:.1f} s', fontsize=11, color='k', 
                         ha='center', va='top', transform=stats_ax.transAxes)
            y_pos -= 0.06
        
        if dp is not None:
            stats_ax.text(0.5, y_pos, f'Dp: {dp:.1f}°', fontsize=11, color='k', 
                         ha='center', va='top', transform=stats_ax.transAxes)
            y_pos -= 0.06
        
        if energy_fraction is not None:
            stats_ax.text(0.5, y_pos, f'Energy: {energy_fraction:.1f}%', fontsize=11, 
                         color='k', ha='center', va='top', transform=stats_ax.transAxes)
            y_pos -= 0.06
        
        # Add partition information if available (for total spectrum plots)
        if show_partitions:
            y_pos -= 0.03
            # Add separator line
            stats_ax.plot([0.1, 0.9], [y_pos, y_pos], 'k-', linewidth=1.2, 
                         transform=stats_ax.transAxes)
            y_pos -= 0.06
            
            # Calculate total energy from significant partitions
            total_sig_energy = sum(partition_results['energy'][i] for i in range(1, len(partition_results['energy'])) 
                                  if partition_results['energy'][i] > min_energy_threshold)
            
            # Show each significant partition
            for p in range(1, min(4, len(partition_results['Hs']))):
                if partition_results['energy'][p] > min_energy_threshold:
                    hs_p = partition_results['Hs'][p]
                    tp_p = partition_results['Tp'][p]
                    dp_p = partition_results['Dp'][p]
                    energy_frac = (partition_results['energy'][p] / total_sig_energy) * 100
                    
                    stats_ax.text(0.5, y_pos, f'P{p}', fontsize=10, color='k',
                                 ha='center', va='top', weight='bold', transform=stats_ax.transAxes)
                    y_pos -= 0.05
                    stats_ax.text(0.5, y_pos, f'Hs={hs_p:.2f}m Tp={tp_p:.1f}s', 
                                 fontsize=9, color='k', ha='center', va='top', transform=stats_ax.transAxes)
                    y_pos -= 0.05
                    stats_ax.text(0.5, y_pos, f'Dp={dp_p:.0f}° E={energy_frac:.1f}%', 
                                 fontsize=9, color='k', ha='center', va='top', transform=stats_ax.transAxes)
                    y_pos -= 0.06

    # Colorbar (horizontal) - Linear scale 0 to 60
    colorbar_label = 'm²·s·rad⁻¹'
    norm = mpl.colors.Normalize(vmin=0, vmax=60)
    sm = ScalarMappable(cmap='rainbow', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, orientation='horizontal', fraction=0.046, pad=0.08, 
                        ax=ax, extend='max', aspect=30)
    cbar.set_label(colorbar_label, fontsize=12)
    cbar.ax.tick_params(labelsize=11)
    # Ticks every 10 units
    cbar.set_ticks(np.arange(0, 70, 10))

    # Manual adjustment
    fig.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.12)

    return fig, ax


# ============================================================================
# PRINTING FUNCTIONS
# ============================================================================

def print_case_header(idx, total_cases, ref_id, target_time_dt, data_type='sar'):
    """Print header of case being processed"""
    print(f"\n{'='*60}")
    print(f"Processing case {idx + 1}/{total_cases}")
    print(f"{'='*60}")
    if data_type == 'ndbc':
        print(f"Station ID: {ref_id}")
        print(f"Target time: {target_time_dt}")
    else:
        print(f"Reference ID: {ref_id}")
        print(f"Target SAR time: {target_time_dt}")


def print_time_match_info(closest_time, itime, time_diff_hours):
    """Print information of matching temporal"""
    print(f"Closest WW3 time: {closest_time}")
    print(f"Time index (itime): {itime}")
    print(f"Time difference: {time_diff_hours:.2f} hours")


def print_partitioning_summary(n_peaks_initial, n_partitions_final):
    """Print summary of process of partitioning"""
    print("\n" + "="*70)
    print(" SPECTRAL PARTITIONING - PROCESS SUMMARY")
    print("="*70)
    print(f"🔍 Spectral peaks initially identified: {n_peaks_initial}")
    print(f"🔗 After merging nearby systems: {n_partitions_final} partition(s)")
    print("="*70)


def print_partitioning_results(results, min_energy_threshold):
    """Print results detailed of partitioning"""
    n_partitions = count_significant_partitions(results, min_energy_threshold)
    
    print("\n" + "="*70)
    print(" PARTITIONING RESULTS")
    print("="*70)
    print(f"Number of partitions found: {n_partitions}")
    print("─"*70)
    
    # Show each partition
    partition_count = 0
    for i in range(1, len(results['Hs'])):
        if results['energy'][i] > min_energy_threshold:
            partition_count += 1
            energy_fraction = (results['energy'][i] / results['total_m0']) * 100
            
            print(f"\nPartition {partition_count}:")
            print(f"  Hs = {results['Hs'][i]:.2f} m")
            print(f"  Tp = {results['Tp'][i]:.2f} s")
            print(f"  Dp = {results['Dp'][i]:.0f}°")
            print(f"  Energy: {results['energy'][i]:.4f} m²")
            print(f"  Energy fraction: {energy_fraction:.1f}%")
    
    # Show integrated total
    print("\n" + "─"*70)
    print(f"Integrated total:")
    print(f"  Hs = {results['total_Hs']:.2f} m")
    print(f"  Tp = {results['total_Tp']:.2f} s")
    print(f"  Dp = {results['total_Dp']:.0f}°")
    print("="*70)


def print_save_confirmation(output_path, df_results):
    """Print confirmation of saving"""
    print(f"\n✓ Results saved to: {output_path}")
    print(f"\nColumns in CSV: {list(df_results.columns)}")
    print(f"\nPreview:")
    print(df_results.T)


# ============================================================================
# PROCESSAMENTO DE DADOS
# ============================================================================

def create_partition_data_dict(ref, selected_time, lon, lat, wnd, wnddir, file_path, 
                                results, min_energy_threshold):
    """
    Create dictionary with results of partitioning
    
    Returns:
    --------
    dict: Dictionary ready for conversion in DataFrame
    """
    moments = results['moments']
    m0_total = moments['total'][0]
    m1_total = moments['total'][1]
    m2_total = moments['total'][2]
    
    data = {
        'reference_id': ref,
        'obs_time': selected_time,
        'longitude': float(lon),
        'latitude': float(lat),
        'source_file': os.path.basename(file_path),
        
        # Spectrum total
        'total_energy': results['total_m0'],
        'total_Hs': results['total_Hs'],
        'total_Tp': results['total_Tp'],
        'total_Dp': results['total_Dp'],
        'total_m0': m0_total,
        'total_m1': m1_total,
        'total_m2': m2_total,
    }
    
    # Add partition data (up to 3 partitions)
    for p in range(1, 4):
        if p < len(results['Hs']) and results['energy'][p] > min_energy_threshold:
            data[f'P{p}_energy'] = results['energy'][p]
            data[f'P{p}_Hs'] = results['Hs'][p]
            data[f'P{p}_Tp'] = results['Tp'][p]
            data[f'P{p}_Dp'] = results['Dp'][p]
            data[f'P{p}_m0'] = moments['m0'][p]
            data[f'P{p}_m1'] = moments['m1'][p]
            data[f'P{p}_m2'] = moments['m2'][p]
        else:
            # Fill with zeros if partition doesn't exist
            data[f'P{p}_energy'] = 0.0
            data[f'P{p}_Hs'] = 0.0
            data[f'P{p}_Tp'] = 0.0
            data[f'P{p}_Dp'] = 0.0
            data[f'P{p}_m0'] = 0.0
            data[f'P{p}_m1'] = 0.0
            data[f'P{p}_m2'] = 0.0
    
    # Add wind data at the end
    data['wnd'] = float(wnd) if wnd is not None else np.nan
    data['wnddir'] = float(wnddir) if wnddir is not None else np.nan
    
    return data


def save_partition_results(ref_id, selected_time, data, output_dir, data_type='sar'):
    """
    Save results of partitioning in file CSV
    
    Returns:
    --------
    tuple: (output_path, df_results)
    """
    # Create nome of file
    date_time_formatted = selected_time.strftime('%Y%m%d-%H%M%S')
    
    if data_type == 'ndbc':
        # For NDBC: use station_id as string
        output_filename = f'ww3_{ref_id}_{date_time_formatted}.csv'
    else:
        # For SAR: use reference_id as integer with padding
        output_filename = f'ww3_{ref_id:03d}_{date_time_formatted}.csv'
    
    output_path = os.path.join(output_dir, output_filename)
    
    # Create DataFrame e save
    df_results = pd.DataFrame([data])
    df_results.to_csv(output_path, index=False, float_format='%.6f')
    
    return output_path, df_results


def process_single_case(row, idx, total_cases, output_dir, data_type='sar'):
    """
    Process a single case WW3
    
    Parameters:
    -----------
    row : pandas.Series
        Row of CSV with information of case
    idx : int
        Index of case current
    total_cases : int
        Number total of cases a process
    output_dir : str
        Directory of output for results
    data_type : str
        Type of data: 'sar' or 'ndbc'
    """
    # Extract information of case based on data type
    if data_type == 'ndbc':
        ref_id = str(row['station_id'])
        # Check if CSV has time column (ndbc_time, ww3_time, obs_time, etc.)
        time_col = None
        for col in ['ndbc_time', 'ww3_time', 'obs_time', 'time']:
            if col in row.index and pd.notna(row[col]):
                time_col = col
                break
        
        if time_col:
            target_time_dt = pd.to_datetime(row[time_col])
        else:
            target_time_dt = None  # Will process multiple times
        
        file_path = f"{WW3_DATA_PATH}/ww3_{ref_id}.nc"
    elif data_type == 'cfosat':
        # CFOSAT data uses sar_ref and cfosat_time columns
        ref_id = int(row['sar_ref'])
        target_time_str = row['cfosat_time']
        target_time_dt = pd.to_datetime(target_time_str)
        file_path = f'{WW3_DATA_PATH}/ww3_cfosat{ref_id:05d}_2020_spec.nc'
    else:  # sar
        ref_id = int(row['ref'])
        target_time_str = row['time']
        target_time_dt = pd.to_datetime(target_time_str)
        file_path = f'{WW3_DATA_PATH}/ww3_sar{ref_id:04d}_2020_spec.nc'

    # Print header
    print_case_header(idx, total_cases, ref_id, target_time_dt if target_time_dt else "Multiple times", data_type)
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"⚠ WW3 file not found: {file_path}")
        return
    
    # If target_time_dt is None, process multiple times based on TIME_INTERVAL_HOURS
    # Otherwise, find and process the closest time
    if target_time_dt is None:
        # Open file to get all times
        with xr.open_dataset(file_path) as ds:
            times = pd.to_datetime(ds.time.values)
            
            # Filter times based on configured interval
            # Keep times where hour is divisible by TIME_INTERVAL_HOURS
            filtered_times = [t for t in times if t.hour % TIME_INTERVAL_HOURS == 0]
            filtered_indices = [i for i, t in enumerate(times) if t.hour % TIME_INTERVAL_HOURS == 0]
            
            interval_str = f"{TIME_INTERVAL_HOURS}-hourly" if TIME_INTERVAL_HOURS > 1 else "hourly"
            print(f"Processing {len(filtered_times)} time steps ({interval_str}) for {data_type.upper()} {ref_id}")
            print(f"Total available: {len(times)} (sampling every {TIME_INTERVAL_HOURS}h)")
            
            # Process each filtered time
            for idx, (itime, obs_time) in enumerate(zip(filtered_indices, filtered_times)):
                print(f"\n  Time {idx+1}/{len(filtered_times)}: {obs_time}")
                
                # Load spectrum for this time
                E2d, freq, dirs, dirs_rad, lon, lat, wnd, wnddir = load_ww3_spectrum(file_path, itime)
                
                # Process this time step
                process_time_step(ref_id, obs_time, E2d, freq, dirs, dirs_rad, 
                                lon, lat, wnd, wnddir, file_path, output_dir, data_type)
        return
    else:
        # Find closest time to target
        itime, closest_time, time_diff_hours = find_closest_time(file_path, target_time_dt)
        print_time_match_info(closest_time, itime, time_diff_hours)
        
        # Load spectrum
        E2d, freq, dirs, dirs_rad, lon, lat, wnd, wnddir = load_ww3_spectrum(file_path, itime)
        
        # Process single time step
        process_time_step(ref_id, closest_time, E2d, freq, dirs, dirs_rad,
                        lon, lat, wnd, wnddir, file_path, output_dir, data_type)


def process_time_step(ref_id, selected_time, E2d, freq, dirs, dirs_rad,
                     lon, lat, wnd, wnddir, file_path, output_dir, data_type='sar'):
    """
    Process a single time step (common for SAR and NDBC)
    """
    
    # Apply partitioning with WW3-specific parameters
    # WW3 has moderate resolution, so we use:
    # - Conservative threshold (99.0%) to avoid detecting noise peaks
    # - Higher merge factor (0.6) to merge nearby systems more aggressively
    results = partition_spectrum(
        E2d, freq, dirs_rad,
        threshold_mode='adaptive',
        threshold_percentile=THRESHOLD_PERCENTILE,
        merge_factor=MERGE_FACTOR,
        max_partitions=MAX_PARTITIONS
    )
    
    if results is None:
        print("    ⚠ No spectral peaks identified!")
        return
    
    # Calculate threshold and count partitions
    min_energy_threshold = MIN_ENERGY_THRESHOLD_FRACTION * results['total_m0']
    n_peaks_initial = len(results['peaks'])
    n_partitions_final = count_significant_partitions(results, min_energy_threshold)
    
    # Print results (condensed for NDBC multi-time processing)
    if data_type == 'ndbc':
        print(f"    Peaks: {n_peaks_initial} → Partitions: {n_partitions_final} | "
              f"Hs={results['total_Hs']:.2f}m Tp={results['total_Tp']:.1f}s")
    else:
        print_partitioning_summary(n_peaks_initial, n_partitions_final)
        print_partitioning_results(results, min_energy_threshold)
    
    # Create and save results
    data = create_partition_data_dict(ref_id, selected_time, lon, lat, wnd, wnddir, file_path,
                                      results, min_energy_threshold)
    output_path, df_results = save_partition_results(ref_id, selected_time, data, output_dir, data_type)
    
    # Generate and save figures (if enabled)
    if GENERATE_FIGURES:
        # Only generate figures for SAR/CFOSAT (not for NDBC multi-time to avoid too many plots)
        if data_type in ['sar', 'cfosat']:
            print(f"\n📊 Generating 2D spectrum figure...")
            
            # Plot total spectrum with partition information
            fig_title = f'WW3 Total Spectrum - {data_type.upper()} {ref_id}'
            fig, ax = plot_partition_spectrum(
                E2d, freq, np.degrees(dirs_rad), title=fig_title,
                hs=results['total_Hs'], tp=results['total_Tp'], dp=results['total_Dp'],
                selected_time=selected_time,
                partition_results=results,
                min_energy_threshold=min_energy_threshold
            )
            
            # Save figure alongside CSV
            dt = pd.to_datetime(selected_time)
            date_time_formatted = dt.strftime('%Y%m%d-%H%M%S')
            
            if data_type == 'cfosat':
                fig_filename = f'ww3_{ref_id}_{date_time_formatted}.png'
            else:  # sar
                fig_filename = f'ww3_{ref_id}_{date_time_formatted}.png'
            
            fig_path = os.path.join(output_dir, fig_filename)
            fig.savefig(fig_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"  ✓ Saved spectrum figure: {fig_filename}")
    
    # Print confirmation (condensed for NDBC)
    if data_type == 'sar':
        print_save_confirmation(output_path, df_results)
    else:
        print(f"    ✓ Saved: {os.path.basename(output_path)}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    # Setup
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load cases
    df = pd.read_csv(CSV_PATH)
    total_cases = len(df)
    
    # Detect data type (SAR or NDBC)
    data_type = detect_data_type(df)
    
    print(f"{'='*60}")
    print(f"WW3 SPECTRAL PARTITIONING PROCESSOR")
    print(f"{'='*60}")
    print(f"Config file: {args.config}")
    print(f"Data type: {data_type.upper()}")
    print(f"Total cases to process: {total_cases}")
    print(f"WW3 data path: {WW3_DATA_PATH}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"CSV file: {CSV_PATH}")
    print(f"Parameters: threshold={THRESHOLD_PERCENTILE}%, merge={MERGE_FACTOR}")
    
    # Process each case
    for idx, row in df.iterrows():
        try:
            process_single_case(row, idx, total_cases, OUTPUT_DIR, data_type)
        except Exception as e:
            print(f"\n❌ Error processing case {idx + 1}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n{'='*60}")
    print(f"✓ Processing complete!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
