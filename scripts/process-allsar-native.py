"""
SAR Native Partition Extraction Script

This script extracts native partition parameters (Hs, Tp, Dp) from SAR L3 files
using the partition_domain_spec variable and integrated spectral parameters.
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
from utils import load_sar_spectrum, calculate_wave_parameters

case = 'surigae'
case = 'lee'
case = 'freddy'
#case = 'all'

# Configuration
OUTPUT_DIR = f'../data/{case}/partition-native'
SAR_DATA_PATH = f'/Users/jtakeo/data/sentinel1ab/{case}'
CSV_PATH = f'../auxdata/sar_matches_{case}_track_3day.csv'
GROUP_NAME = "obs_params"


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def extract_native_partitions(ds, E2d, freq, dirs_rad, index):
    """
    Extract native SAR partition parameters from partition_domain_spec
    
    Returns:
        list of dicts with partition info (partition_id, Hs, Tp, Dp, energy, fraction)
    """
    # Get partition map from SAR L3 file
    partition_map = ds.partition_domain_spec[:, :, index].values
    
    # Find unique partition IDs (excluding 0 = background)
    unique_partitions = np.unique(partition_map)
    partition_ids = unique_partitions[unique_partitions > 0]
    
    # Calculate total energy for fraction computation
    hs_total, tp_total, dp_total, m0_total, _, _, _, _ = calculate_wave_parameters(E2d, freq, dirs_rad)
    
    partitions = []
    
    for pid in partition_ids:
        # Create mask for this partition
        mask = (partition_map == pid)
        
        # Create spectrum with only this partition
        E2d_partition = np.zeros_like(E2d)
        E2d_partition[mask] = E2d[mask]
        
        # Calculate integrated parameters for this partition
        hs_part, tp_part, dp_part, m0_part, _, _, _, _ = calculate_wave_parameters(
            E2d_partition, freq, dirs_rad
        )
        
        energy_part = m0_part
        energy_fraction = (energy_part / m0_total) if m0_total > 0 else 0
        
        # Count spectral points in this partition
        n_points = np.sum(mask)
        pct_points = 100 * n_points / partition_map.size
        
        partitions.append({
            'partition_id': int(pid),
            'Hs': hs_part,
            'Tp': tp_part,
            'Dp': dp_part,
            'energy': energy_part,
            'energy_fraction': energy_fraction,
            'n_spectral_points': n_points,
            'pct_spectral_points': pct_points
        })
    
    return partitions, m0_total, hs_total, tp_total, dp_total


# ============================================================================
# DATA PROCESSING FUNCTIONS
# ============================================================================

def load_sar_data(file_path, index):
    """Load SAR dataset and extract spectrum"""
    ds = xr.open_dataset(file_path, group=GROUP_NAME)
    E2d, freq, dirs, dirs_rad, selected_time = load_sar_spectrum(ds, date_time=None, index=index)
    return ds, E2d, freq, dirs, dirs_rad, selected_time


def extract_sar_metadata(ds, index):
    """Extract metadata from SAR dataset"""
    quality_flag = ds.L2_partition_quality_flag[index].values
    actual_time = pd.to_datetime(ds.time[index].values)
    lon = ds.longitude[index].values if 'longitude' in ds else None
    lat = ds.latitude[index].values if 'latitude' in ds else None
    
    # Extract SAR native partition parameters (VAVH, VTPK, VPED)
    vavh_sar = ds.VAVH[index].values
    vtpk_sar = ds.VTPK[index].values
    vped_sar = ds.VPED[index].values
    
    return quality_flag, actual_time, lon, lat, vavh_sar, vtpk_sar, vped_sar


def create_partition_data_dict(ref, index, quality_flag, date_time, lon, lat, 
                                file_path, m0_total, hs_total, tp_total, dp_total,
                                partitions, vavh_sar, vtpk_sar, vped_sar):
    """Create data dictionary with native partition results"""
    
    data = {
        'reference_id': ref,
        'obs_index': index,
        'quality_flag': quality_flag,
        'obs_time': date_time,
        'longitude': float(lon),
        'latitude': float(lat),
        'source_file': os.path.basename(file_path),
        'total_energy': m0_total,
        'total_Hs': hs_total,
        'total_Tp': tp_total,
        'total_Dp': dp_total,
        'n_partitions': len(partitions),
        # SAR L3 native values (VAVH, VTPK, VPED)
        'sar_native_Hs': vavh_sar,
        'sar_native_Tp': vtpk_sar,
        'sar_native_Dp': vped_sar,
    }
    
    # Add up to 3 partitions
    for p in range(1, 4):
        if p <= len(partitions):
            part = partitions[p-1]
            data.update({
                f'P{p}_id': part['partition_id'],
                f'P{p}_Hs': part['Hs'],
                f'P{p}_Tp': part['Tp'],
                f'P{p}_Dp': part['Dp'],
                f'P{p}_energy': part['energy'],
                f'P{p}_energy_fraction': part['energy_fraction'],
                f'P{p}_n_points': part['n_spectral_points'],
                f'P{p}_pct_points': part['pct_spectral_points'],
            })
        else:
            data.update({
                f'P{p}_id': 0,
                f'P{p}_Hs': 0.0,
                f'P{p}_Tp': 0.0,
                f'P{p}_Dp': 0.0,
                f'P{p}_energy': 0.0,
                f'P{p}_energy_fraction': 0.0,
                f'P{p}_n_points': 0,
                f'P{p}_pct_points': 0.0,
            })
    
    return data


def save_partition_results(ref, index, date_time, data, output_dir):
    """Save partition results to CSV file"""
    dt = pd.to_datetime(date_time)
    date_time_formatted = dt.strftime('%Y%m%d-%H%M%S')
    output_filename = f'sar_native_{ref:03d}_{index}_{date_time_formatted}.csv'
    output_path = os.path.join(output_dir, output_filename)
    
    df_results = pd.DataFrame([data])
    df_results.to_csv(output_path, index=False, float_format='%.6f')
    
    return output_path, df_results


# ============================================================================
# PRINTING FUNCTIONS
# ============================================================================

def print_progress(idx, total, ref, index, file_name):
    """Print progress header"""
    print(f"\n{'='*60}")
    print(f"Processing case {idx + 1}/{total}")
    print(f"{'='*60}")
    print(f"File: {file_name}")
    print(f"Index: {index}, Ref: {ref}")


def print_partition_results(partitions, vavh_sar, vtpk_sar, vped_sar):
    """Print native partition results"""
    n_partitions = len(partitions)
    
    print(f"\n{'='*70}")
    print(" SAR NATIVE PARTITIONS (from partition_domain_spec)")
    print(f"{'='*70}")
    print(f"Native partitions identified: {n_partitions}")
    
    for i, part in enumerate(partitions, 1):
        pct = part['energy_fraction'] * 100
        print(f"\nPartition {i} (ID={part['partition_id']}):")
        print(f"  Hs={part['Hs']:.2f}m, Tp={part['Tp']:.2f}s, Dp={part['Dp']:.0f}°")
        print(f"  Energy={part['energy']:.6f}m² ({pct:.1f}%)")
        print(f"  Spectral points: {part['n_spectral_points']} ({part['pct_spectral_points']:.1f}%)")
    
    print(f"\n{'─'*70}")
    print("SAR L3 file values (VAVH, VTPK, VPED):")
    print(f"  Hs={vavh_sar:.2f}m, Tp={vtpk_sar:.2f}s, Dp={vped_sar:.0f}°")
    print(f"{'='*70}")


def print_save_info(output_path, n_partitions):
    """Print save confirmation"""
    print(f"\n✓ Saved: {output_path}")
    print(f"✓ Native partitions: {n_partitions}")


# ============================================================================
# MAIN PROCESSING
# ============================================================================

def process_single_case(row, idx, total_cases, output_dir):
    """Process a single SAR case"""
    # Extract parameters
    file_path = os.path.join(SAR_DATA_PATH, row['filename'])
    index = int(row['obs_index'])
    date_time = row['time']
    ref = int(row['ref'])
    
    print_progress(idx, total_cases, ref, index, os.path.basename(file_path))
    
    # Load data
    ds, E2d, freq, dirs, dirs_rad, selected_time = load_sar_data(file_path, index)
    
    # Extract metadata
    quality_flag, actual_time, lon, lat, vavh_sar, vtpk_sar, vped_sar = extract_sar_metadata(ds, index)
    
    print(f"Location: ({lon:.2f}°E, {lat:.2f}°N), Quality: {quality_flag}")
    
    # Extract native partitions
    partitions, m0_total, hs_total, tp_total, dp_total = extract_native_partitions(
        ds, E2d, freq, dirs_rad, index
    )
    
    if len(partitions) == 0:
        print("⚠ No native partitions identified!")
        ds.close()
        return
    
    # Print results
    print_partition_results(partitions, vavh_sar, vtpk_sar, vped_sar)
    
    # Save results
    data = create_partition_data_dict(ref, index, quality_flag, date_time, lon, lat,
                                      file_path, m0_total, hs_total, tp_total, dp_total,
                                      partitions, vavh_sar, vtpk_sar, vped_sar)
    output_path, _ = save_partition_results(ref, index, date_time, data, output_dir)
    
    print_save_info(output_path, len(partitions))
    ds.close()


def main():
    """Main execution function"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(CSV_PATH)
    
    print(f"{'='*60}")
    print("SAR NATIVE PARTITION EXTRACTOR")
    print(f"{'='*60}")
    print(f"Total cases: {len(df)}")
    
    for idx, row in df.iterrows():
        try:
            process_single_case(row, idx, len(df), OUTPUT_DIR)
        except Exception as e:
            print(f"\n❌ Error in case {idx + 1}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n{'='*60}")
    print("✓ Processing complete!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
