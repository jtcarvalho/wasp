"""
NDBC Buoy Spectral Partitioning

Process all available NDBC buoy spectral data and apply partitioning.
This script processes all stations, files, and timesteps in the NDBC directory.

Usage:
------
python 03_partition_ndbc.py --config config.yaml

Arguments:
  --config: Path to configuration YAML file (default: config.yaml)

Configuration:
  - Edit TIME_FREQUENCY_HOURS below to control processing frequency
  - 0 or 1 = process all timesteps
  - 3 = process every 3 hours
  - 6 = process every 6 hours
  - 12 = process every 12 hours
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import xarray as xr
from pathlib import Path
from glob import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable

# Import from wasp package
from wasp.wave_params import calculate_wave_parameters
from wasp.partition import partition_spectrum
from wasp.io_ndbc import load_ndbc_spectrum
from wasp.utils import load_config

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='NDBC Spectral Partitioning with configurable frequency',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default config
  python 03_partition_ndbc.py
  
  # Use custom config file
  python 03_partition_ndbc.py --config my_config.yaml
        """
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

CASE_NAME = 'all'

# ====================================================================
# TIME FREQUENCY CONFIGURATION - EDIT THIS VALUE
# ====================================================================
# Process every N hours:
#   0 or 1 = process all timesteps (every hour)
#   3 = process every 3 hours
#   6 = process every 6 hours (RECOMMENDED for faster processing)
#   12 = process every 12 hours
# ====================================================================
TIME_FREQUENCY_HOURS = 6  
# ====================================================================

# Directories
NDBC_BASE_DIR = CONFIG['paths']['ndbc_data']  # Base directory with station folders

# Partitioning parameters (from config.yaml)
MAX_PARTITIONS = CONFIG['partitioning']['ndbc']['max_partitions']
MIN_ENERGY_FRACTION = CONFIG['partitioning']['ndbc']['min_energy_fraction']
MAX_TIME_DIFF_HOURS = CONFIG['partitioning']['ndbc']['max_time_diff_hours']
THRESHOLD_PERCENTILE = CONFIG['partitioning']['ndbc']['threshold_percentile']
MERGE_FACTOR = CONFIG['partitioning']['ndbc']['merge_factor']

# Year filter (None = all years, or specify year like 2020)
YEAR_FILTER = 2020  # Only process 2020 files

# Create output directory with parameters in name
OUTPUT_DIR = f'../data/{CASE_NAME}/partition-ndbc-{THRESHOLD_PERCENTILE}-{MERGE_FACTOR}'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Figure generation flag
GENERATE_FIGURES = CONFIG.get('processing', {}).get('generate_figures', True)

# Plotting parameters (from config.yaml)
PLOT_VMIN = CONFIG['plotting'].get('spectrum_vmin', 0.5)
PLOT_VMAX = CONFIG['plotting'].get('spectrum_vmax', 60.0)
PLOT_STEP = CONFIG['plotting'].get('spectrum_step', 0.5)
PLOT_PERIOD_MAX = CONFIG['plotting'].get('period_max', 25)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def create_partition_data_dict(station_id, timestamp, lon, lat, file_path, 
                                results, min_energy_threshold):
    """
    Create dictionary with partitioning results
    
    Returns:
    --------
    dict: Dictionary ready for DataFrame conversion
    """
    moments = results['moments']
    m0_total = moments['total'][0]
    m1_total = moments['total'][1]
    m2_total = moments['total'][2]
    
    data = {
        'station_id': station_id,
        'obs_time': timestamp,
        'longitude': float(lon),
        'latitude': float(lat),
        'source_file': os.path.basename(file_path),
        
        # Total spectrum
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
    
    return data


def save_partition_results(station_id, timestamp, data, output_dir):
    """
    Save partitioning results to individual CSV file
    
    Returns:
    --------
    str: Path to saved file
    """
    # Create filename
    date_time_formatted = timestamp.strftime('%Y%m%d-%H%M%S')
    output_filename = f'ndbc_{station_id}_{date_time_formatted}.csv'
    output_path = os.path.join(output_dir, output_filename)
    
    # Create DataFrame and save
    df_results = pd.DataFrame([data])
    df_results.to_csv(output_path, index=False, float_format='%.6f')
    
    return output_path


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
# MAIN PROCESSING
# ============================================================================

def process_ndbc_file(file_path, station_id, output_dir, time_freq_hours=1):
    """
    Process all timesteps in a single NDBC file and apply partitioning.
    Saves individual CSV and PNG files for each timestep.
    
    Parameters:
    -----------
    file_path : Path
        Path to NDBC NetCDF file
    station_id : str
        Station identifier
    output_dir : str
        Output directory for results
    time_freq_hours : int
        Process every N hours (default: 1 = all timesteps)
    
    Returns:
    --------
    tuple: (success_count, failed_count, total_timesteps, processed_timesteps)
    """
    print(f"\n{'='*70}")
    print(f"Station: {station_id}")
    print(f"File: {file_path.name}")
    print(f"{'='*70}")
    
    try:
        # Open NDBC dataset
        ds = xr.open_dataset(file_path)
        
        # Get station location
        lon = float(ds.longitude.values.item() if ds.longitude.size == 1 else ds.longitude.values[0])
        lat = float(ds.latitude.values.item() if ds.latitude.size == 1 else ds.latitude.values[0])
        
        # Get all times
        times = pd.to_datetime(ds.time.values)
        n_times = len(times)
        
        print(f"Location: {lon:.3f}°E, {lat:.3f}°N")
        print(f"Time range: {times[0]} to {times[-1]}")
        print(f"Total timesteps: {n_times}")
        
        # Filter timesteps based on frequency
        if time_freq_hours > 1:
            # Select timesteps at the specified frequency
            time_indices = []
            for i, t in enumerate(times):
                if i == 0 or (t.hour % time_freq_hours == 0 and t.minute == 0):
                    time_indices.append(i)
            
            if not time_indices:
                # If no exact matches, sample every N timesteps
                time_indices = list(range(0, n_times, time_freq_hours))
            
            print(f"Processing every {time_freq_hours} hours: {len(time_indices)} timesteps selected")
        else:
            # 0 or 1 = process all timesteps
            time_indices = list(range(n_times))
            print(f"Processing all timesteps")
        
        # Process each selected timestep
        success_count = 0
        failed_count = 0
        
        for itime in time_indices:
            timestamp = times[itime]
            
            try:
                # Load spectrum
                spectrum_result = load_ndbc_spectrum(ds, itime)
                
                if spectrum_result is None:
                    failed_count += 1
                    continue
                
                E2d, freq, dirs, dirs_rad, _, _ = spectrum_result
                
                # Apply partitioning with NDBC-specific parameters
                results = partition_spectrum(
                    E2d, freq, dirs_rad,
                    threshold_mode='adaptive',
                    threshold_percentile=THRESHOLD_PERCENTILE,
                    merge_factor=MERGE_FACTOR,
                    max_partitions=MAX_PARTITIONS
                )
                
                if results is None:
                    failed_count += 1
                    continue
                
                # Calculate threshold
                min_energy_threshold = MIN_ENERGY_FRACTION * results['total_m0']
                
                # Count significant partitions
                n_partitions = sum(1 for i in range(1, len(results['Hs'])) 
                                  if results['energy'][i] > min_energy_threshold)
                
                if n_partitions == 0:
                    failed_count += 1
                    continue
                
                # Skip if output already exists (resume support)
                date_time_formatted = timestamp.strftime('%Y%m%d-%H%M%S')
                output_filename = f'ndbc_{station_id}_{date_time_formatted}.csv'
                output_path = os.path.join(output_dir, output_filename)
                if os.path.exists(output_path):
                    success_count += 1
                    continue

                # Create and save partition data
                data = create_partition_data_dict(
                    station_id, timestamp, lon, lat, file_path,
                    results, min_energy_threshold
                )
                
                output_path = save_partition_results(
                    station_id, timestamp, data, output_dir
                )
                
                # Generate and save figure if enabled
                if GENERATE_FIGURES:
                    try:
                        fig_title = f'NDBC {station_id} - Total Spectrum'
                        fig, ax = plot_partition_spectrum(
                            E2d, freq, dirs, title=fig_title,
                            hs=results['total_Hs'], 
                            tp=results['total_Tp'], 
                            dp=results['total_Dp'],
                            selected_time=timestamp,
                            partition_results=results,
                            min_energy_threshold=min_energy_threshold
                        )
                        
                        # Save figure alongside CSV
                        fig_filename = f'ndbc_{station_id}_{date_time_formatted}.png'
                        fig_path = os.path.join(output_dir, fig_filename)
                        fig.savefig(fig_path, dpi=150, bbox_inches='tight')
                        plt.close(fig)
                    except Exception as fig_err:
                        print(f"  ⚠ Error generating figure: {fig_err}")
                
                success_count += 1
                    
            except Exception as e:
                failed_count += 1
                if itime == 0:  # Print error for first timestep
                    print(f"  ⚠ Error at timestep {itime}: {e}")
                continue
        
        ds.close()
        
        print(f"\n✓ Processed: {success_count}/{len(time_indices)} selected timesteps")
        print(f"✗ Failed: {failed_count}/{len(time_indices)} selected timesteps")
        print(f"  (Total timesteps in file: {n_times})")
        print(f"{'='*70}")
        
        return success_count, failed_count, n_times, len(time_indices)
        
    except Exception as e:
        print(f"  ✖ Error opening file: {e}")
        import traceback
        traceback.print_exc()
        return 0, 0, 0, 0


def find_all_ndbc_files(base_dir, year_filter=None):
    """
    Find all NDBC NetCDF files in station subdirectories.
    
    Parameters:
    -----------
    base_dir : str or Path
        Base directory containing station folders
    year_filter : int or None
        If provided, only include files for this year (e.g., 2020)
    
    Returns:
    --------
    list of tuple
        List of (station_id, file_path) tuples
    """
    base_path = Path(base_dir)
    station_files = []
    
    print(f"Scanning base directory: {base_path}")
    print(f"Directory exists: {base_path.exists()}")
    if year_filter:
        print(f"Year filter: {year_filter}")
    
    if not base_path.exists():
        print(f"ERROR: Base directory does not exist!")
        return station_files
    
    # Look for station subdirectories
    subdirs = [d for d in base_path.iterdir() if d.is_dir()]
    print(f"Found {len(subdirs)} subdirectories")
    
    for station_dir in sorted(subdirs):
        station_id = station_dir.name
        
        # Find all NetCDF files in station directory
        nc_files = list(station_dir.glob('*.nc'))
        nc4_files = list(station_dir.glob('*.nc4'))
        all_files = nc_files + nc4_files
        
        # Apply year filter if specified
        if year_filter:
            year_pattern = f"w{year_filter}"
            all_files = [f for f in all_files if year_pattern in f.stem]
        
        if all_files:
            print(f"  {station_id}: {len(all_files)} file(s)")
        
        for file_path in sorted(all_files):
            station_files.append((station_id, file_path))
    
    return station_files


# ============================================================================
# MAIN
# ============================================================================

def main():
    """
    Process all available NDBC files and timesteps, saving individual CSV and PNG files.
    """
    print("="*70)
    print("NDBC SPECTRA PARTITIONING - ALL FILES")
    print("="*70)
    print(f"Case: {CASE_NAME}")
    print(f"Config file: {args.config}")
    if TIME_FREQUENCY_HOURS <= 1:
        print(f"Time frequency: All timesteps (every hour)")
    else:
        print(f"Time frequency: Every {TIME_FREQUENCY_HOURS} hour(s)")
    print(f"NDBC directory: {NDBC_BASE_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Generate figures: {GENERATE_FIGURES}")
    print("="*70)
    
    # Find all NDBC files
    print("\nScanning for NDBC files...")
    station_files = find_all_ndbc_files(NDBC_BASE_DIR, year_filter=YEAR_FILTER)
    
    if not station_files:
        print(f"\n⚠️  No NDBC files found in: {NDBC_BASE_DIR}")
        print("\nPlease check that the directory contains station subdirectories with .nc or .nc4 files.")
        return
    
    print(f"\n✓ Found {len(station_files)} NDBC files")
    
    # Group by station
    stations = {}
    for station_id, file_path in station_files:
        if station_id not in stations:
            stations[station_id] = []
        stations[station_id].append(file_path)
    
    print(f"  Unique stations: {len(stations)}")
    for station_id, files in sorted(stations.items()):
        print(f"    {station_id}: {len(files)} file(s)")
    
    # Process each file
    total_files = len(station_files)
    total_success = 0
    total_failed = 0
    total_timesteps = 0
    total_processed = 0
    
    for idx, (station_id, file_path) in enumerate(station_files, 1):
        print(f"\n[{idx}/{total_files}]")
        
        try:
            success, failed, n_times, n_processed = process_ndbc_file(
                file_path, station_id, OUTPUT_DIR, 
                time_freq_hours=TIME_FREQUENCY_HOURS
            )
            total_success += success
            total_failed += failed
            total_timesteps += n_times
            total_processed += n_processed
                
        except Exception as e:
            print(f"  ✖ Unexpected error: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Summary
    print(f"\n{'='*70}")
    print("PROCESSING SUMMARY")
    print(f"{'='*70}")
    print(f"Files processed: {total_files}")
    print(f"Total timesteps in files: {total_timesteps}")
    if TIME_FREQUENCY_HOURS <= 1:
        print(f"Timesteps selected: {total_processed} (all)")
    else:
        print(f"Timesteps selected (every {TIME_FREQUENCY_HOURS}h): {total_processed}")
    print(f"✓ Successful: {total_success}")
    print(f"✗ Failed: {total_failed}")
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"  Individual CSV files: {total_success}")
    if GENERATE_FIGURES:
        print(f"  Individual PNG files: {total_success}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
