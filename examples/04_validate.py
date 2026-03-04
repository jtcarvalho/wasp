"""
STEP 4: Validate CFOSAT vs WW3

This script compares CFOSAT and WW3 partitioning results, matching
similar partitions and generating validation metrics and plots.

Workflow:
1. Load CFOSAT and WW3 partitioning CSVs
2. Match partitions based on Tp and Dp
3. Generate paired partition files (partition1.csv, partition2.csv, partition3.csv)
4. Create dispersion plots comparing CFOSAT vs WW3
5. Calculate statistical metrics (bias, RMSE, correlation)
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving figures

# ============================================================================
# CONFIGURATION
# ============================================================================

# case = 'surigae'
#case = 'lee'
# case = 'freddy'
case = 'all'

# Directories
WW3_DIR = f'../data/{case}/partition-ww3/'
CFOSAT_DIR = f'../data/cfosat/partition-cfosat/'
SAR_DIR = CFOSAT_DIR  # Alias for compatibility
OUTPUT_DIR = f'../output/cfosat/'

# Filtros
QUALITY_FLAG_OPTIONS = [0]  # Apenas data SAR of high quality (0 = best)

# Partition matching criteria
TP_MIN_THRESHOLD = 10.0  # SAR not reliable for Tp < 10s (wind sea)
MAX_TIME_DIFF_HOURS = 1.0  # Maximum acceptable temporal difference [hours]

# Physical constraints (reject impossible matches - different wave systems)
TP_MAX_DIFF = 4.0   # Maximum Tp difference [s] - systems must be similar
DP_MAX_DIFF = 45.0  # Maximum Dp difference [degrees] - systems must come from similar direction

# Score-based ranking (among physically plausible matches)
MAX_SCORE = None  # Optional: maximum acceptable score (None = no limit)

# Matching weights (relative importance of each parameter)
# Score = TP_WEIGHT*(tp_diff/tp_sar) + DP_WEIGHT*(dp_diff/180°) + HS_WEIGHT*(hs_diff/hs_sar)
# Lower score = better match. All parameters use relative/normalized errors.
TP_WEIGHT = 1.0  # Weight for Tp difference (1.0 = standard importance)
DP_WEIGHT = 1.0  # Weight for Dp difference (1.0 = standard importance) 
HS_WEIGHT = 0.3  # Weight for Hs difference (0.3 = moderate, 0=disabled, 1=equal to Tp)

# Wind speed filtering (set None to disable filter)
WND_MIN = 0  # Minimum wind speed [m/s] (e.g., 0)
WND_MAX = 10  # Maximum wind speed [m/s] (e.g., 10)

# Peak period filtering for P1 (set None to disable filter)
TP_MIN = 12  # Minimum Tp for P1 [s] (e.g., 10)
TP_MAX = 18  # Maximum Tp for P1 [s] (e.g., 20)

# Debug options
DEBUG_SCORE = False  # Print detailed score calculations for first few matches

# Plot limits
PLOT_LIMITS = {
    'Hs': (0, 8),
    'Tp': (4, 20),
    'Dp': (0, 360)
}

# Variables to compare
COMPARISON_VARIABLES = [
    'total_Hs', 'total_Tp', 'total_Dp',
    'P1_Hs', 'P1_Tp', 'P1_Dp',
    'P2_Hs', 'P2_Tp', 'P2_Dp',
    'P3_Hs', 'P3_Tp', 'P3_Dp'
]

# ============================================================================
# I/O AND FILTERING FUNCTIONS
# ============================================================================

def load_and_filter_sar(sar_file, quality_flags=None):
    """
    Load data SAR e filtra por quality_flag
    
    Parameters:
    -----------
    sar_file : Path
        Path of file SAR
    quality_flags : list, optional
        Lista of quality flags aceitos (standard: [0])
    
    Returns:
    --------
    pd.DataFrame: Dados SAR filtrados
    """
    if quality_flags is None:
        quality_flags = QUALITY_FLAG_OPTIONS
    
    df = pd.read_csv(sar_file)
    
    if 'quality_flag' in df.columns:
        df_filtered = df[df['quality_flag'].isin(quality_flags)]
    else:
        df_filtered = df
    
    return df_filtered


def load_ww3_files_dict(ww3_dir):
    """
    Load todos os files WW3 e cria dictionary indexado por reference_id
    
    IMPORTANT: Each reference_id can have multiple simulated timestamps.
    Armazenamos uma lista of DataFrames for fazer matching temporal posterior.
    
    Parameters:
    -----------
    ww3_dir : str or Path
        Directory contendo files WW3
    
    Returns:
    --------
    dict: {reference_id: list of DataFrames}
    """
    ww3_files = list(Path(ww3_dir).glob('ww3_*.csv'))
    ww3_dict = {}
    
    for ww3_file in ww3_files:
        df_ww3 = pd.read_csv(ww3_file)
        if 'reference_id' in df_ww3.columns and len(df_ww3) > 0:
            ref_id = df_ww3['reference_id'].iloc[0]
            
            # Store as list to support multiple timestamps
            if ref_id not in ww3_dict:
                ww3_dict[ref_id] = []
            ww3_dict[ref_id].append(df_ww3)
    
    return ww3_dict


# ============================================================================
# PARTITION MATCHING FUNCTIONS
# ============================================================================

def validatete_temporal_match(sar_time, ww3_time, max_diff_hours=MAX_TIME_DIFF_HOURS):
    """
    Validate se o matching temporal between SAR e WW3 é acceptable
    
    Parameters:
    -----------
    sar_time : str or pd.Timestamp
        Timestamp of the observation SAR
    ww3_time : str or pd.Timestamp
        Timestamp of the simulation WW3
    max_diff_hours : float
        Difference temporal maximum acceptable [hours]
    
    Returns:
    --------
    tuple: (is_valid, time_diff_hours)
        is_valid: bool - True se difference temporal é acceptable
        time_diff_hours: float - Difference temporal in hours
    """
    sar_dt = pd.to_datetime(sar_time)
    ww3_dt = pd.to_datetime(ww3_time)
    
    time_diff_hours = abs((sar_dt - ww3_dt).total_seconds() / 3600.0)
    is_valid = time_diff_hours <= max_diff_hours
    
    return is_valid, time_diff_hours


def compute_angular_difference(angle1, angle2):
    """
    Calculates difference angular minimum considering circularity (0-360°)
    
    Parameters:
    -----------
    angle1, angle2 : float
        Angles in degrees
    
    Returns:
    --------
    float: Difference angular minimum [0, 180]
    """
    return abs((angle1 - angle2 + 180) % 360 - 180)


def extract_partitions_from_row(row, prefix=''):
    """
    Extracts data of partitions (Hs, Tp, Dp) of uma line of DataFrame
    
    Parameters:
    -----------
    row : pd.Series
        Row of DataFrame
    prefix : str
        Prefixo of colunas (vazio por standard)
    
    Returns:
    --------
    list: Lista of dicts with data of partitions
    """
    partitions = []
    
    for p in [1, 2, 3]:
        hs = row.get(f'{prefix}P{p}_Hs', np.nan)
        tp = row.get(f'{prefix}P{p}_Tp', np.nan)
        dp = row.get(f'{prefix}P{p}_Dp', np.nan)
        
        # Add only if Hs is not NaN (partition exists)
        if not np.isnan(hs):
            partitions.append({
                'partition': p,
                'hs': hs,
                'tp': tp,
                'dp': dp
            })
    
    return partitions


def find_best_match(sar_partitions, ww3_partitions, sar_pnum, 
                    tp_min=TP_MIN_THRESHOLD,
                    tp_max_diff=TP_MAX_DIFF, dp_max_diff=DP_MAX_DIFF,
                    tp_weight=TP_WEIGHT, dp_weight=DP_WEIGHT, hs_weight=HS_WEIGHT,
                    max_score=MAX_SCORE,
                    debug=False):
    """
    Find a best partition WW3 correspondente a uma partition SAR using:
    1) Physical constraints to reject impossible matches (different wave systems)
    2) Score-based ranking among physically plausible candidates
    
    Score = tp_weight*(tp_diff/tp_sar) + dp_weight*(dp_diff/180°) + hs_weight*(hs_diff/hs_sar)
    
    Note: 
    - Physical filter: tp_diff <= tp_max_diff AND dp_diff <= dp_max_diff
    - Score uses relative/normalized errors (scale-independent)
    - Tp and Hs: relative error (% of observed value)
    - Dp: normalized by 180° (maximum angular difference)
    - Partitions with Tp < tp_min are rejected (SAR unreliable for wind sea)
    
    Parameters:
    -----------
    sar_partitions : list
        Lista of partitions SAR
    ww3_partitions : list
        Lista of partitions WW3
    sar_pnum : int
        Number of the partition SAR (1, 2, ou 3)
    tp_min : float
        Tp minimum for considerar [s]
    tp_max_diff : float
        Maximum Tp difference to accept [s] - physical constraint
    dp_max_diff : float
        Maximum Dp difference to accept [degrees] - physical constraint
    tp_weight : float
        Weight for Tp difference (1.0 = standard)
    dp_weight : float
        Weight for Dp difference (1.0 = standard)
    hs_weight : float
        Weight for Hs penalty (0=disabled, 0.3=moderate, 1=high)
    max_score : float or None
        Maximum acceptable score (None = no limit)
    
    Returns:
    --------
    tuple: (sar_data, ww3_data) ou (sar_data, None) se not houver match
    """
    # Findr partition SAR
    sar_data = next((p for p in sar_partitions if p['partition'] == sar_pnum), None)
    if not sar_data:
        return None, None
    
    # Reject SAR partitions with Tp < tp_min (SAR unreliable for wind sea)
    if not np.isnan(sar_data['tp']) and sar_data['tp'] < tp_min:
        return sar_data, None
    
    # Buscar best match WW3
    best_ww3 = None
    min_score = None
    candidates = []  # Store all candidates with their scores
    
    for ww3 in ww3_partitions:
        # Skip if values are NaN
        if (np.isnan(ww3['tp']) or np.isnan(sar_data['tp']) or
            np.isnan(ww3['dp']) or np.isnan(sar_data['dp'])):
            continue
        
        # Reject WW3 partitions with Tp < tp_min (cannot be validateted with SAR)
        if ww3['tp'] < tp_min:
            continue
        
        # Calculate differences
        tp_diff = abs(ww3['tp'] - sar_data['tp'])
        dp_diff = compute_angular_difference(ww3['dp'], sar_data['dp'])
        
        # STEP 1: SCORE CALCULATION for ALL candidates (no filter yet)
        # 
        # Each parameter contributes to the total score:
        # - Tp penalty: relative error (tp_diff/tp_sar) × TP_WEIGHT
        # - Dp penalty: normalized by max angular diff (dp_diff/180°) × DP_WEIGHT  
        # - Hs penalty: relative error (hs_diff/hs_sar) × HS_WEIGHT
        #
        # Lower score = better match
        
        tp_penalty = tp_diff / max(sar_data['tp'], 1.0)  # relative error
        dp_penalty = dp_diff / 180.0  # normalized by max angular difference (180°)
        
        score = tp_weight * tp_penalty + dp_weight * dp_penalty
        
        # Add Hs penalty if enabled and valid
        hs_diff = None
        if hs_weight > 0 and not np.isnan(sar_data['hs']) and not np.isnan(ww3['hs']):
            hs_diff = abs(ww3['hs'] - sar_data['hs'])
            hs_penalty = hs_diff / max(sar_data['hs'], 0.1)  # relative error
            score += hs_weight * hs_penalty
        
        # Store candidate with score and differences
        candidates.append({
            'ww3': ww3,
            'score': score,
            'tp_diff': tp_diff,
            'dp_diff': dp_diff,
            'hs_diff': hs_diff,
            'tp_penalty': tp_penalty,
            'dp_penalty': dp_penalty
        })
    
    # Sort candidates by score (lower = better)
    candidates.sort(key=lambda x: x['score'])
    
    # STEP 2: PHYSICAL CONSTRAINT FILTER
    # Among best-scored candidates, pick first that passes physical constraints
    for idx, candidate in enumerate(candidates):
        ww3 = candidate['ww3']
        tp_diff = candidate['tp_diff']
        dp_diff = candidate['dp_diff']
        score = candidate['score']
        
        # Debug: print all candidates in rank order
        if debug:
            rank = idx + 1
            print(f"    Rank {rank}: WW3 P{ww3['partition']} (score={score:.4f}):")
            print(f"      Tp: {ww3['tp']:.1f}s vs SAR {sar_data['tp']:.1f}s → diff={tp_diff:.2f}s → penalty={candidate['tp_penalty']:.4f}")
            print(f"      Dp: {ww3['dp']:.0f}° vs SAR {sar_data['dp']:.0f}° → diff={dp_diff:.1f}° → penalty={candidate['dp_penalty']:.4f}")
            if hs_weight > 0 and candidate['hs_diff'] is not None:
                hs_penalty = candidate['hs_diff'] / max(sar_data['hs'], 0.1)
                print(f"      Hs: {ww3['hs']:.2f}m vs SAR {sar_data['hs']:.2f}m → diff={candidate['hs_diff']:.2f}m → penalty={hs_penalty:.4f}")
        
        # Check physical constraints
        if tp_diff <= tp_max_diff and dp_diff <= dp_max_diff:
            # Check optional max_score threshold
            if max_score is None or score <= max_score:
                best_ww3 = ww3
                min_score = score
                if debug:
                    print(f"      ✓ SELECTED (passes constraints: Tp≤{tp_max_diff}s, Dp≤{dp_max_diff}°)")
                break
            elif debug:
                print(f"      ✗ Rejected: score {score:.4f} > max_score {max_score}")
        elif debug:
            print(f"      ✗ Rejected: Tp_diff={tp_diff:.2f}s>{tp_max_diff}s OR Dp_diff={dp_diff:.1f}°>{dp_max_diff}°")
    
    return sar_data, best_ww3


def create_match_record(ref_id, sar_row, ww3_row, sar_match, ww3_match, time_diff_hours):
    """
    Create registro of partitions paired
    
    Returns:
    --------
    dict: Registro with data SAR, WW3 e differences
    """
    tp_diff = abs(sar_match['tp'] - ww3_match['tp'])
    dp_diff = compute_angular_difference(sar_match['dp'], ww3_match['dp'])
    
    return {
        'reference_id': ref_id,
        'sar_obs_time': sar_row.get('obs_time', ''),
        'ww3_obs_time': ww3_row.get('obs_time', ''),
        'time_diff_hours': time_diff_hours,
        'longitude': sar_row.get('longitude', np.nan),
        'latitude': sar_row.get('latitude', np.nan),
        'quality_flag': sar_row.get('quality_flag', np.nan),
        'sar_partition': sar_match['partition'],
        'ww3_partition': ww3_match['partition'],
        'sar_Hs': sar_match['hs'],
        'sar_Tp': sar_match['tp'],
        'sar_Dp': sar_match['dp'],
        'ww3_Hs': ww3_match['hs'],
        'ww3_Tp': ww3_match['tp'],
        'ww3_Dp': ww3_match['dp'],
        'tp_diff': tp_diff,
        'dp_diff': dp_diff
    }


# ============================================================================
# CREATING PAIRED PARTITION FILES
# ============================================================================

def create_partition_matches(quality_flags=None, wnd_min=WND_MIN, wnd_max=WND_MAX,
                             tp_min=TP_MIN, tp_max=TP_MAX):
    """
    Create files of partitions paired (partition1.csv, partition2.csv, partition3.csv)
    
    Parameters:
    -----------
    quality_flags : list or None
        Quality flags to accept (None = use QUALITY_FLAG_OPTIONS)
    wnd_min : float or None
        Minimum wind speed to consider [m/s]. None = no minimum filter
    wnd_max : float or None
        Maximum wind speed to consider [m/s]. None = no maximum filter
    tp_min : float or None
        Minimum Tp for P1 to consider [s]. None = no minimum filter
    tp_max : float or None
        Maximum Tp for P1 to consider [s]. None = no maximum filter
    
    Returns:
    --------
    dict: {partition_num: list of matches}
    """
    if quality_flags is None:
        quality_flags = QUALITY_FLAG_OPTIONS
    
    # Look for CFOSAT files (cfosat_*.csv)
    cfosat_files = list(Path(CFOSAT_DIR).glob('cfosat_*.csv'))
    sar_files = cfosat_files  # Alias for compatibility with rest of code
    ww3_dict = load_ww3_files_dict(WW3_DIR)
    
    print(f"Found {len(cfosat_files)} CFOSAT files and {len(ww3_dict)} WW3 files")
    
    # Check if wnd column exists in WW3 files (debug)
    if len(ww3_dict) > 0:
        first_ref = list(ww3_dict.keys())[0]
        first_ww3_df = ww3_dict[first_ref][0]
        print(f"\nDEBUG: Columns in WW3 file: {list(first_ww3_df.columns)}")
        if 'wnd' in first_ww3_df.columns:
            wnd_sample = first_ww3_df['wnd'].iloc[0]
            print(f"DEBUG: Sample wnd value: {wnd_sample}")
            print(f"DEBUG: wnd range in sample: [{first_ww3_df['wnd'].min():.2f}, {first_ww3_df['wnd'].max():.2f}] m/s")
        else:
            print("WARNING: 'wnd' column NOT found in WW3 files!")
        print()
    
    # Storage for paired partitions and statistics
    partition_matches = {1: [], 2: [], 3: []}
    total_sar_files = 0
    temporal_match_valid = 0
    temporal_match_rejected = 0
    spatial_match_not_found = 0
    wind_filter_rejected = 0
    tp_filter_rejected = 0
    
    # Process each file SAR
    for sar_file in sar_files:
        df_sar = load_and_filter_sar(sar_file, quality_flags=quality_flags)
        
        if len(df_sar) == 0:
            continue
        
        total_sar_files += 1
        
        # Obter reference_id of file SAR
        if 'reference_id' not in df_sar.columns:
            continue
        
        ref_id = df_sar['reference_id'].iloc[0]
        
        # Findr data WW3 correspondentes (matching espacial)
        if ref_id not in ww3_dict:
            spatial_match_not_found += 1
            continue
        
        ww3_list = ww3_dict[ref_id]  # Lista of DataFrames WW3 for este ref_id
        
        # Extract time SAR
        if len(df_sar) == 0:
            continue
            
        sar_row = df_sar.iloc[0]
        sar_time = sar_row.get('obs_time', '')
        sar_time_dt = pd.to_datetime(sar_time)
        
        # Find closest WW3 temporally (with optional Tp filtering)
        best_ww3 = None
        best_time_diff = float('inf')
        
        for df_ww3 in ww3_list:
            if len(df_ww3) == 0:
                continue
            ww3_row = df_ww3.iloc[0]
            
            # Apply Tp filtering if configured (check P1 Tp BEFORE temporal matching)
            if tp_min is not None or tp_max is not None:
                tp_p1_value = ww3_row.get('P1_Tp', np.nan)
                
                # Debug first few cases
                if total_sar_files <= 3:
                    print(f"    WW3 candidate: P1_Tp={tp_p1_value:.2f}s, range=[{tp_min}, {tp_max}]")
                
                # Skip if Tp is NaN or outside the specified range
                if np.isnan(tp_p1_value):
                    if total_sar_files <= 3:
                        print(f"      → Skipped: P1_Tp is NaN")
                    continue
                
                if tp_min is not None and tp_p1_value < tp_min:
                    if total_sar_files <= 3:
                        print(f"      → Skipped: P1_Tp ({tp_p1_value:.2f}) < tp_min ({tp_min})")
                    continue
                
                if tp_max is not None and tp_p1_value > tp_max:
                    if total_sar_files <= 3:
                        print(f"      → Skipped: P1_Tp ({tp_p1_value:.2f}) > tp_max ({tp_max})")
                    continue
                
                if total_sar_files <= 3:
                    print(f"      → Accepted: P1_Tp within range")
            
            ww3_time = ww3_row.get('obs_time', '')
            ww3_time_dt = pd.to_datetime(ww3_time)
            
            time_diff = abs((sar_time_dt - ww3_time_dt).total_seconds() / 3600.0)
            
            if time_diff < best_time_diff:
                best_time_diff = time_diff
                best_ww3 = (df_ww3, ww3_row, time_diff)
        
        # Validatete if best temporal match is acceptable
        if best_ww3 is None:
            # No WW3 found (could be filtered by Tp or just not available)
            if tp_min is not None or tp_max is not None:
                tp_filter_rejected += 1
            else:
                spatial_match_not_found += 1
            continue
            
        df_ww3, ww3_row, time_diff_hours = best_ww3
        
        if time_diff_hours > MAX_TIME_DIFF_HOURS:
            temporal_match_rejected += 1
            continue
        
        # Apply wind speed filter if configured
        if wnd_min is not None or wnd_max is not None:
            wnd_value = ww3_row.get('wnd', np.nan)
            
            # Debug: Print first few wind values
            if total_sar_files <= 3:
                print(f"  DEBUG ref_id={ref_id}: wnd={wnd_value:.2f}, range=[{wnd_min}, {wnd_max}]")
            
            # Skip if wnd is NaN or outside the specified range
            if np.isnan(wnd_value):
                wind_filter_rejected += 1
                if total_sar_files <= 3:
                    print(f"    → Rejected: wnd is NaN")
                continue
            
            if wnd_min is not None and wnd_value < wnd_min:
                wind_filter_rejected += 1
                if total_sar_files <= 3:
                    print(f"    → Rejected: wnd ({wnd_value:.2f}) < wnd_min ({wnd_min})")
                continue
            
            if wnd_max is not None and wnd_value > wnd_max:
                wind_filter_rejected += 1
                if total_sar_files <= 3:
                    print(f"    → Rejected: wnd ({wnd_value:.2f}) > wnd_max ({wnd_max})")
                continue
            
            if total_sar_files <= 3:
                print(f"    → Accepted: wnd within range")
        
        temporal_match_valid += 1
        
        # Extract partitions
        sar_partitions = extract_partitions_from_row(sar_row)
        ww3_partitions = extract_partitions_from_row(ww3_row)
        
        # Findr melhores matches for each partition SAR
        for sar_pnum in [1, 2, 3]:
            # Enable debug for first 2 files if DEBUG_SCORE is True
            show_debug = DEBUG_SCORE and total_sar_files <= 2
            if show_debug:
                sar_p_hs = sar_row.get(f'P{sar_pnum}_Hs', np.nan)
                sar_p_tp = sar_row.get(f'P{sar_pnum}_Tp', np.nan)
                sar_p_dp = sar_row.get(f'P{sar_pnum}_Dp', np.nan)
                print(f"\n  DEBUG: Matching SAR P{sar_pnum} (Hs={sar_p_hs:.2f}m, Tp={sar_p_tp:.1f}s, Dp={sar_p_dp:.0f}°)")
            
            sar_match, ww3_match = find_best_match(
                sar_partitions, ww3_partitions, sar_pnum,
                tp_max_diff=TP_MAX_DIFF, dp_max_diff=DP_MAX_DIFF,
                tp_weight=TP_WEIGHT, dp_weight=DP_WEIGHT, hs_weight=HS_WEIGHT,
                max_score=MAX_SCORE,
                debug=show_debug
            )
            
            if show_debug and ww3_match:
                print(f"  ✓ Best match: WW3 P{ww3_match['partition']} (Hs={ww3_match['hs']:.2f}m, Tp={ww3_match['tp']:.1f}s, Dp={ww3_match['dp']:.0f}°)\n")
            
            if sar_match and ww3_match:
                # Apply Tp filter to the MATCHED PARTITIONS (not just the file)
                # Only save pairs where BOTH partitions have Tp in the specified range
                if tp_min is not None or tp_max is not None:
                    sar_tp = sar_match['tp']
                    ww3_tp = ww3_match['tp']
                    
                    # Skip if either Tp is outside the range
                    if not np.isnan(sar_tp) and not np.isnan(ww3_tp):
                        if tp_min is not None and (sar_tp < tp_min or ww3_tp < tp_min):
                            continue  # Skip this match
                        if tp_max is not None and (sar_tp > tp_max or ww3_tp > tp_max):
                            continue  # Skip this match
                
                match_record = create_match_record(
                    ref_id, sar_row, ww3_row, sar_match, ww3_match, time_diff_hours
                )
                partition_matches[sar_pnum].append(match_record)
    
    # Print matching statistics
    print(f"\n{'='*70}")
    print(f" TEMPORAL MATCHING STATISTICS")
    print(f"{'='*70}")
    print(f"Total SAR files processed: {total_sar_files}")
    print(f"Spatial matches found (same reference_id): {temporal_match_valid + temporal_match_rejected + wind_filter_rejected + tp_filter_rejected}")
    print(f"Spatial matches NOT found: {spatial_match_not_found}")
    print(f"\nTemporal validatetion (max diff = {MAX_TIME_DIFF_HOURS} hour):")
    print(f"  ✓ Valid temporal matches: {temporal_match_valid + wind_filter_rejected + tp_filter_rejected}")
    print(f"  ✗ Rejected (time diff > {MAX_TIME_DIFF_HOURS}h): {temporal_match_rejected}")
    
    if wnd_min is not None or wnd_max is not None:
        wnd_range = f"[{wnd_min if wnd_min is not None else 0}, {wnd_max if wnd_max is not None else '∞'}] m/s"
        print(f"\nWind speed filtering (range: {wnd_range}):")
        print(f"  ✓ Valid wind speed matches: {temporal_match_valid + tp_filter_rejected}")
        print(f"  ✗ Rejected (outside wind range): {wind_filter_rejected}")
    
    if tp_min is not None or tp_max is not None:
        tp_range = f"[{tp_min if tp_min is not None else 0}, {tp_max if tp_max is not None else '∞'}] s"
        print(f"\nPeak period (P1) filtering (range: {tp_range}):")
        print(f"  ✓ Valid Tp matches: {temporal_match_valid}")
        print(f"  ✗ Rejected (outside Tp range): {tp_filter_rejected}")
    
    total_valid = temporal_match_valid
    total_rejected = temporal_match_rejected + wind_filter_rejected + tp_filter_rejected
    
    if total_valid + total_rejected > 0:
        valid_pct = 100 * total_valid / (total_valid + total_rejected)
        print(f"\nOverall success rate: {valid_pct:.1f}%")
    print(f"{'='*70}")
    
    # Save partitions paired
    save_partition_matches(partition_matches)
    
    # Print detailed matching statistics
    print_matching_statistics(partition_matches)
    
    return partition_matches


def print_matching_statistics(partition_matches):
    """Print detailed statistics about partition matching"""
    print(f"\n{'='*70}")
    print(f" PARTITION MATCHING STATISTICS")
    print(f"{'='*70}")
    
    for sar_pnum in [1, 2, 3]:
        if len(partition_matches[sar_pnum]) == 0:
            print(f"\nSAR P{sar_pnum}: No matches found")
            continue
        
        df = pd.DataFrame(partition_matches[sar_pnum])
        
        # Count how many times each WW3 partition was matched
        ww3_counts = df['ww3_partition'].value_counts().sort_index()
        
        print(f"\nSAR P{sar_pnum} ({len(df)} matches):")
        print(f"  Matched with WW3 partitions:")
        for ww3_p, count in ww3_counts.items():
            pct = 100 * count / len(df)
            print(f"    P{ww3_p}: {count:4d} times ({pct:5.1f}%)")
        
        # Statistics of differences
        print(f"  Match quality:")
        print(f"    Mean Tp diff: {df['tp_diff'].mean():5.2f} ± {df['tp_diff'].std():4.2f} s")
        print(f"    Mean Dp diff: {df['dp_diff'].mean():5.1f} ± {df['dp_diff'].std():4.1f}°")
        
        # Hs statistics
        hs_diff = df['ww3_Hs'] - df['sar_Hs']
        hs_diff_pct = 100 * hs_diff / df['sar_Hs']
        print(f"    Mean Hs diff: {hs_diff.mean():5.2f} ± {hs_diff.std():4.2f} m ({hs_diff_pct.mean():5.1f}%)")
    
    print(f"{'='*70}")
    
    return partition_matches


def save_partition_matches(partition_matches):
    """Save partitions paired in files CSV e print summary"""
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for pnum in [1, 2, 3]:
        if len(partition_matches[pnum]) > 0:
            df_partition = pd.DataFrame(partition_matches[pnum])
            output_file = output_dir / f'partition{pnum}.csv'
            df_partition.to_csv(output_file, index=False)
            
            print(f"\nPartition {pnum}: {len(df_partition)} matches found")
            print(f"  Saved to: {output_file}")
            print(f"  SAR Hs range: [{df_partition['sar_Hs'].min():.2f}, {df_partition['sar_Hs'].max():.2f}] m")
            print(f"  WW3 Hs range: [{df_partition['ww3_Hs'].min():.2f}, {df_partition['ww3_Hs'].max():.2f}] m")
            print(f"  Mean Tp diff: {df_partition['tp_diff'].mean():.2f} s")
            print(f"  Mean Dp diff: {df_partition['dp_diff'].mean():.2f}°")
        else:
            print(f"\nPartition {pnum}: No matches found")


# ============================================================================
# METRICS CALCULATION FUNCTIONS
# ============================================================================

def compute_metrics(obs, model):
    """
    Calculates metrics of comparison between observations e model
    
    Parameters:
    -----------
    obs : array-like
        Valores observados (SAR)
    model : array-like
        Valores of model (WW3)
    
    Returns:
    --------
    dict: Dictionary with metrics (nbias, nrmse, pearson_r, n_points)
    """
    obs = np.array(obs)
    model = np.array(model)
    
    # Remover valores NaN
    mask = ~(np.isnan(obs) | np.isnan(model))
    obs = obs[mask]
    model = model[mask]
    
    if len(obs) == 0:
        return {'nbias': np.nan, 'nrmse': np.nan, 'pearson_r': np.nan, 'n_points': 0}
    
    # Bias normalizado
    bias = np.mean(model - obs)
    nbias = bias / np.mean(obs) if np.mean(obs) != 0 else np.nan
    
    # RMSE normalizado
    rmse = np.sqrt(np.mean((model - obs)**2))
    nrmse = rmse / np.mean(obs) if np.mean(obs) != 0 else np.nan
    
    # Pearson correlation
    if len(obs) > 1:
        pearson_r, _ = pearsonr(obs, model)
    else:
        pearson_r = np.nan
    
    return {
        'nbias': nbias,
        'nrmse': nrmse,
        'pearson_r': pearson_r,
        'n_points': len(obs)
    }


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def setup_plot_axis(ax, var, axis_limits):
    """Configure limits e appearance dos axes"""
    axis_min, axis_max = axis_limits
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3)
    return axis_min, axis_max


def plot_dispersion_and_lines(ax, sar_clean, ww3_clean, axis_min, axis_max):
    """Plot points of dispersion, line 1:1 e regression"""
    # Scatter plot
    ax.scatter(sar_clean, ww3_clean, alpha=0.6, s=100, 
              edgecolors='black', linewidth=0.5)
    
    # Row 1:1
    ax.plot([axis_min, axis_max], [axis_min, axis_max], 
           'r--', linewidth=2, label='1:1 line')
    
    # Row of regression
    if len(sar_clean) > 1:
        z = np.polyfit(sar_clean, ww3_clean, 1)
        p = np.poly1d(z)
        x_line = np.linspace(axis_min, axis_max, 100)
        ax.plot(x_line, p(x_line), 'b-', linewidth=1.5, alpha=0.7,
               label=f'Fit: y={z[0]:.2f}x+{z[1]:.2f}')


def add_metrics_textbox(ax, sar_clean, ww3_clean):
    """Add text box with metrics to plot"""
    # Calculate metrics
    bias = np.mean(ww3_clean - sar_clean)
    nbias = bias / np.mean(sar_clean) if np.mean(sar_clean) != 0 else np.nan
    rmse = np.sqrt(np.mean((ww3_clean - sar_clean)**2))
    nrmse = rmse / np.mean(sar_clean) if np.mean(sar_clean) != 0 else np.nan
    
    if len(sar_clean) > 1:
        pearson_r, _ = pearsonr(sar_clean, ww3_clean)
    else:
        pearson_r = np.nan
    
    # Create texto with metrics
    metrics_text = (
        f'n = {len(sar_clean)}\n'
        f'Bias = {bias:.3f}\n'
        f'NBias = {nbias:.3f}\n'
        f'RMSE = {rmse:.3f}\n'
        f'NRMSE = {nrmse:.3f}\n'
        f'R = {pearson_r:.3f}'
    )
    
    # Posicionar text box
    y_pos = 0.05 if np.mean(ww3_clean) > np.mean(sar_clean) else 0.95
    va = 'bottom' if y_pos == 0.05 else 'top'
    
    ax.text(0.95, y_pos, metrics_text,
           transform=ax.transAxes,
           fontsize=10,
           verticalalignment=va,
           horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))


def plot_single_variable(ax, df, var, var_label):
    """Plot comparison for uma single variable"""
    sar_col = f'sar_{var}'
    ww3_col = f'ww3_{var}'
    
    sar_data = df[sar_col].values
    ww3_data = df[ww3_col].values
    
    # Remover valores NaN
    mask = ~(np.isnan(sar_data) | np.isnan(ww3_data))
    sar_clean = sar_data[mask]
    ww3_clean = ww3_data[mask]
    
    if len(sar_clean) == 0:
        ax.text(0.5, 0.5, 'No valid data',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title(var_label)
        return
    
    # Obter limits dos axes
    axis_limits = PLOT_LIMITS.get(var, (0, max(sar_clean.max(), ww3_clean.max())))
    axis_min, axis_max = setup_plot_axis(ax, var, axis_limits)
    
    # Plotar data e linhas
    plot_dispersion_and_lines(ax, sar_clean, ww3_clean, axis_min, axis_max)
    
    # Add labels e legenda
    ax.set_xlabel(f'SAR {var}', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'WW3 {var}', fontsize=12, fontweight='bold')
    ax.set_title(var_label, fontsize=13, fontweight='bold')
    ax.legend(loc='upper left', fontsize=9)
    
    # Add metrics
    add_metrics_textbox(ax, sar_clean, ww3_clean)


def plot_partition_comparisons():
    """Create dispersion plots for each partition comparing CFOSAT vs WW3"""
    output_dir = Path(OUTPUT_DIR)
    
    # Variables to plot
    variables = [
        ('Hs', 'Significant Wave Height (m)'),
        ('Tp', 'Peak Period (s)'),
        ('Dp', 'Peak Direction (deg)')
    ]
    
    for pnum in [1, 2, 3]:
        partition_file = output_dir / f'partition{pnum}.csv'
        
        if not partition_file.exists():
            print(f"Partition {pnum} file not found, skipping...")
            continue
        
        df = pd.read_csv(partition_file)
        
        if len(df) == 0:
            print(f"Partition {pnum} has in the data, skipping...")
            continue
        
        # Create figura with 3 subplots
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle(f'Partition {pnum} - CFOSAT vs WW3 Comparison (n={len(df)})',
                     fontsize=16, fontweight='bold')
        
        # Plotar each variable
        for idx, (var, var_label) in enumerate(variables):
            plot_single_variable(axes[idx], df, var, var_label)
        
        plt.tight_layout()
        
        # Save figura
        output_file = output_dir / f'partition{pnum}_dispersion.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")
        plt.close()
    
    print("\nAll dispersion plots created successfully!")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main(create_files=True, create_plots=True):
    """
    Main execution function
    
    Parameters:
    -----------
    create_files : bool
        Se True, cria files of partitions paired (partition*.csv)
    create_plots : bool
        Se True, cria dispersion plots
    """
    if create_files:
        # Create files of partitions paired
        print("="*80)
        print("CREATING MATCHED PARTITION FILES (CFOSAT vs WW3)")
        print("="*80)
        partition_matches = create_partition_matches()
    
    if create_plots:
        # Create dispersion plots
        print("\n" + "="*80)
        print("CREATING SCATTER PLOTS")
        print("="*80)
        plot_partition_comparisons()
    
    print("\n" + "="*80)
    print("✓ Validatetion complete!")
    print("="*80)


if __name__ == '__main__':
    # ========================================================================
    # EXECUTION OPTIONS
    # ========================================================================
    
    RUN_CREATE_FILES = True   # Create files partition*.csv (matching SAR/WW3)
    RUN_CREATE_PLOTS = True   # Create dispersion plots
    
    # ========================================================================
    
    print("\n" + "="*80)
    print("VALIDATION CONFIGURATION (CFOSAT vs WW3)")
    print("="*80)
    print(f"Case: {case}")
    print(f"Create partition files: {RUN_CREATE_FILES}")
    print(f"Create dispersion plots:   {RUN_CREATE_PLOTS}")
    print(f"Quality flags: {QUALITY_FLAG_OPTIONS}")
    print(f"Tp min threshold: {TP_MIN_THRESHOLD} s")
    print(f"Physical constraints: Tp_max_diff={TP_MAX_DIFF}s, Dp_max_diff={DP_MAX_DIFF}°")
    print(f"Max score threshold: {MAX_SCORE if MAX_SCORE is not None else 'DISABLED (no limit)'}")
    print(f"Matching weights: Tp={TP_WEIGHT}, Dp={DP_WEIGHT}, Hs={HS_WEIGHT}")
    
    # Print wind filter configuration
    if WND_MIN is not None or WND_MAX is not None:
        wnd_range = f"[{WND_MIN if WND_MIN is not None else 0}, {WND_MAX if WND_MAX is not None else '∞'}] m/s"
        print(f"Wind speed filter: {wnd_range}")
    else:
        print(f"Wind speed filter: DISABLED")
    
    # Print Tp filter configuration
    if TP_MIN is not None or TP_MAX is not None:
        tp_range = f"[{TP_MIN if TP_MIN is not None else 0}, {TP_MAX if TP_MAX is not None else '∞'}] s"
        print(f"Peak period (P1) filter: {tp_range}")
    else:
        print(f"Peak period (P1) filter: DISABLED")
    
    print("="*80 + "\n")
    
    main(create_files=RUN_CREATE_FILES, create_plots=RUN_CREATE_PLOTS)
