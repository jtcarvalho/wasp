"""Partition matching and observation-association helpers.

The first section implements the current descriptor and Hungarian-assignment
interfaces. The latter section contains legacy spatial/temporal helpers for
SAR, NDBC, and WW3 records; these remain in this module for compatibility.
"""

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from scipy.optimize import linear_sum_assignment


def compute_partition_descriptors(E2d, frequencies, directions_rad, mask,
                                  partition_labels=None):
    """Compute physical descriptors for each labelled spectral partition.

    The returned descriptors contain only independent spectral information:
    peak frequency, peak direction, integrated energy, frequency bandwidth,
    directional spreading, and the two-dimensional frequency-direction spreading
    used to normalize inter-system distance. Significant wave height is
    deliberately absent because it is derived from integrated energy.
    """
    E2d = np.asarray(E2d, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    directions_rad = np.asarray(directions_rad, dtype=float)
    mask = np.asarray(mask)

    if E2d.ndim != 2 or mask.shape != E2d.shape:
        raise ValueError("E2d and mask must be two-dimensional arrays with the same shape")
    if E2d.shape != (len(frequencies), len(directions_rad)):
        raise ValueError("E2d shape must be (len(frequencies), len(directions_rad))")
    if len(frequencies) < 2 or len(directions_rad) < 2:
        raise ValueError("at least two frequency and direction bins are required")

    clean_energy = np.where(np.isfinite(E2d) & (E2d >= 0), E2d, 0.0)
    frequency_weights = np.empty_like(frequencies)
    frequency_weights[0] = (frequencies[1] - frequencies[0]) / 2
    frequency_weights[-1] = (frequencies[-1] - frequencies[-2]) / 2
    frequency_weights[1:-1] = (frequencies[2:] - frequencies[:-2]) / 2
    ddir = 2 * np.pi / len(directions_rad)
    cell_energy = clean_energy * frequency_weights[:, np.newaxis] * ddir

    if partition_labels is None:
        partition_labels = np.unique(mask[mask > 0])

    freq_grid, direction_grid = np.meshgrid(frequencies, directions_rad, indexing="ij")
    x_grid = freq_grid * np.cos(direction_grid)
    y_grid = freq_grid * np.sin(direction_grid)
    descriptors = []

    for label in partition_labels:        
        region = (mask == label)

        # Ignore labels that no longer exist after merge/renumbering
        if not np.any(region):
            continue

        energy = float(np.sum(cell_energy[region]))

        # Ignore empty or numerically invalid partitions
        if (not np.isfinite(energy)) or (energy <= 0):
            continue

        peak_i, peak_j = np.unravel_index(
            np.argmax(np.where(region, clean_energy, -np.inf)), clean_energy.shape
        )
        mean_frequency = np.sum(cell_energy[region] * freq_grid[region]) / energy
        bandwidth = np.sqrt(np.sum(
            cell_energy[region] * (freq_grid[region] - mean_frequency)**2
        ) / energy)

        mean_cos = np.sum(cell_energy[region] * np.cos(direction_grid[region])) / energy
        mean_sin = np.sum(cell_energy[region] * np.sin(direction_grid[region])) / energy
        resultant_length = np.clip(np.hypot(mean_cos, mean_sin), 0.0, 1.0)
        directional_spreading = np.sqrt(max(0.0, 2 * (1 - resultant_length)))

        mean_x = np.sum(cell_energy[region] * x_grid[region]) / energy
        mean_y = np.sum(cell_energy[region] * y_grid[region]) / energy
        spectral_spreading = (
            np.sum(cell_energy[region] * (x_grid[region] - mean_x)**2) / energy
            + np.sum(cell_energy[region] * (y_grid[region] - mean_y)**2) / energy
        )

        descriptors.append({
            "partition": int(label),
            "peak_frequency": float(frequencies[peak_i]),
            "peak_direction": float(np.degrees(directions_rad[peak_j]) % 360),
            "energy": energy,
            "bandwidth": float(bandwidth),
            "directional_spreading": float(directional_spreading),
            "spectral_spreading": float(spectral_spreading),
        })

    return descriptors


def _descriptor_value(descriptor, name):
    """Read and validate one required scalar descriptor."""
    try:
        value = float(descriptor[name])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"descriptor must contain numeric '{name}'") from exc
    if not np.isfinite(value):

        print("\n==============================")
        print("INVALID DESCRIPTOR")
        print("field =", name)
        print(descriptor)
        print("==============================")



        raise ValueError(f"descriptor field '{name}' must be finite")
    return value


def _physical_matching_cost(observed, modeled, alpha, beta, gamma, delta, method='probabilistic'):
    """Return the currently implemented cost for one descriptor pair.

    Input descriptors use ``tp`` (seconds), ``dp`` (degrees), positive energy,
    bandwidth, directional spreading, and spectral spreading. The active cost
    is the sum of log-period, circular-direction, log-energy, log-bandwidth,
    and log-directional-spread penalties. ``spectral_spreading`` is validated
    but the calculated Cartesian normalization is not part of the active sum.

    ``alpha``, ``beta``, ``gamma``, and ``delta`` are retained for API
    compatibility with the earlier weighted formulation; they do not currently
    alter the returned value. ``method='probabilistic'`` additionally returns a
    component breakdown.
    """
    obs_tp = _descriptor_value(observed, "tp")
    mod_tp = _descriptor_value(modeled, "tp")
    obs_frequency = 1.0/obs_tp
    mod_frequency = 1.0/mod_tp

    obs_direction = np.radians(_descriptor_value(observed, "dp"))
    mod_direction = np.radians(_descriptor_value(modeled, "dp"))
    obs_energy = _descriptor_value(observed, "energy")
    mod_energy = _descriptor_value(modeled, "energy")
    obs_bandwidth = _descriptor_value(observed, "bandwidth")
    mod_bandwidth = _descriptor_value(modeled, "bandwidth")
    obs_directional_spread = _descriptor_value(observed, "directional_spreading")
    mod_directional_spread = _descriptor_value(modeled, "directional_spreading")
    obs_spectral_spread = _descriptor_value(observed, "spectral_spreading")
    mod_spectral_spread = _descriptor_value(modeled, "spectral_spreading")

    if min(obs_energy, mod_energy) <= 0:
        raise ValueError("descriptor energy must be positive")
    if min(obs_bandwidth, mod_bandwidth) < 0:
        raise ValueError("descriptor bandwidth must be non-negative")
    if min(obs_directional_spread, mod_directional_spread) < 0:
        raise ValueError("descriptor directional_spreading must be non-negative")
    if min(obs_spectral_spread, mod_spectral_spread) < 0:
        raise ValueError("descriptor spectral_spreading must be non-negative")

    obs_x, obs_y = obs_frequency * np.cos(obs_direction), obs_frequency * np.sin(obs_direction)
    mod_x, mod_y = mod_frequency * np.cos(mod_direction), mod_frequency * np.sin(mod_direction)
    distance_squared = (obs_x - mod_x)**2 + (obs_y - mod_y)**2

    # This is the Hanson & Phillips-style physical normalization: distance in
    # frequency-direction space relative to the systems' combined spreading.
    numerical_floor = 1.0e-3  # physical floor
    dnorm = np.sqrt(distance_squared / max(obs_spectral_spread + mod_spectral_spread,
                                            numerical_floor))
    energy_term = abs(np.log(obs_energy / mod_energy))
    bandwidth_term = abs(np.log(max(obs_bandwidth, numerical_floor)
                                / max(mod_bandwidth, numerical_floor)))
    directional_spread_term = abs(np.log(max(obs_directional_spread, numerical_floor)
                                         / max(mod_directional_spread, numerical_floor)))

    method == "probabilistic"
    # Relative peak-period difference (logarithmic)
    sigma_tp = 0.1      # configurable parameter
    dtheta = np.arctan2(np.sin(obs_direction-mod_direction), np.cos(obs_direction-mod_direction))
    sigma_theta = np.clip(np.sqrt(obs_directional_spread**2 +mod_directional_spread**2),np.deg2rad(5),np.deg2rad(35))
    sigma_theta2 = sigma_theta**2

    tp_cost = (0.5 * (np.log(obs_tp / mod_tp) / sigma_tp) ** 2)
    direction_cost = 0.5*(dtheta**2)/sigma_theta2
    # Penalização para diferenças angulares muito grandes
    theta_max = np.deg2rad(40.0)

    direction_cost *= (1.0 + (np.abs(dtheta) / theta_max)**4)
    # sigma_f2 = max(obs_bandwidth**2 + mod_bandwidth**2, numerical_floor)

    # sigma_theta2 = max(obs_directional_spread**2 +
    #                    mod_directional_spread**2,
    #                    numerical_floor)

    # frequency_cost = 0.5*(obs_frequency-mod_frequency)**2/sigma_f2
    # direction_cost = 0.5*(dtheta**2)/sigma_theta2

    energy_cost = 0.5*(np.log(obs_energy/mod_energy))**2
    bandwidth_cost = 0.5*(np.log(max(obs_bandwidth,numerical_floor)/
                                max(mod_bandwidth,numerical_floor)))**2
    spread_cost = 0.5*(np.log(max(obs_directional_spread,numerical_floor)/
                                max(mod_directional_spread,numerical_floor)))**2

    #cost = (frequency_cost + direction_cost + energy_cost + bandwidth_cost + spread_cost)
    cost = (tp_cost + direction_cost + energy_cost + bandwidth_cost + spread_cost)
    breakdown={
        "tp_cost":float(tp_cost),
         "direction_cost":float(direction_cost),
        "energy_cost":float(energy_cost),
        "bandwidth_cost":float(bandwidth_cost),
        "spread_cost":float(spread_cost),
    }


    return (float(cost), breakdown) if method=='probabilistic' else float(cost)
    #return (alpha * dnorm + beta * energy_term + gamma * bandwidth_term
    #        + delta * directional_spread_term)




def get_tp_tolerance(obs_tp,
                     tp_break_1=10.0,
                     tp_break_2=14.0,
                     tp_rel_low=0.20,
                     tp_rel_mid=0.25,
                     tp_rel_high=0.30):
    """Return the adaptive physical tolerance for peak period matching."""
    if obs_tp < tp_break_1:
        return tp_rel_low * obs_tp
    elif tp_break_1 <= obs_tp < tp_break_2:
        return tp_rel_mid * obs_tp
    else:
        return tp_rel_high * obs_tp


def match_partition_properties(observed_partitions, modeled_partitions,
                               alpha=1.0, beta=1.0, gamma=1.0, delta=1.0,
                               method='probabilistic', probability_threshold=0.05,
                               physical_filter=False,
                               max_direction_deg=60.0,
                               tp_break_1=10.0,
                               tp_break_2=14.0,
                               tp_rel_low=0.20,
                               tp_rel_mid=0.25,
                               tp_rel_high=0.30,
                               max_direction_cost=2.0):
    """Associate property dictionaries with the current PCSPM implementation.

    Each descriptor must contain numeric ``tp``, ``dp``, ``energy``,
    ``bandwidth``, ``directional_spreading``, and ``spectral_spreading`` values;
    ``partition`` is optional metadata. A cost matrix is built with
    :func:`_physical_matching_cost`, augmented with dummy assignments, and
    solved by :func:`scipy.optimize.linear_sum_assignment`.

    In probabilistic mode, the dummy cost is the greater of the probability
    threshold cost and the 95th percentile of finite candidate costs. Legacy
    non-probabilistic mode uses the 75th percentile. Candidate pairs can also be
    excluded by the optional period/direction physical filter. Systems not
    accepted by the assignment are returned unchanged in ``unmatched_observed``
    and ``unmatched_modeled``.

    The four weight arguments and ``max_direction_cost`` are compatibility
    parameters and do not affect the current active cost.
    """
    # weights = np.asarray([alpha, beta, gamma, delta], dtype=float)
    # if np.any(~np.isfinite(weights)) or np.any(weights < 0):
    #     raise ValueError("alpha, beta, gamma, and delta must be finite and non-negative")

    observed = list(observed_partitions)
    modeled = list(modeled_partitions)
    costs = np.empty((len(observed), len(modeled)), dtype=float)
    candidate_costs=[]
    matching_ranking=[]
    for obs_index, observed_system in enumerate(observed):

        obs_tp = observed_system["tp"]
        tp_tol = get_tp_tolerance(
            obs_tp,
            tp_break_1=tp_break_1,
            tp_break_2=tp_break_2,
            tp_rel_low=tp_rel_low,
            tp_rel_mid=tp_rel_mid,
            tp_rel_high=tp_rel_high,
        )

        for mod_index, modeled_system in enumerate(modeled):

            dtheta = float(abs(np.degrees(np.arctan2(
                np.sin(np.radians(observed_system["dp"]-modeled_system["dp"])),
                np.cos(np.radians(observed_system["dp"]-modeled_system["dp"]))
            ))))

            dtp = abs(observed_system["tp"]-modeled_system["tp"])

            if physical_filter and (
                dtheta > max_direction_deg or
                dtp > tp_tol
            ):
                costs[obs_index, mod_index] = 1.0e12
                candidate_costs.append({
                    "observed_index": int(obs_index),
                    "modeled_index": int(mod_index),
                    "observed_partition": observed_system.get("partition", obs_index),
                    "modeled_partition": modeled_system.get("partition", mod_index),
                    "cost": np.nan,
                    "probability": np.nan,
                    "direction_difference_deg": dtheta,
                    "tp_difference": dtp,
                    "filter_reason": "PRE_PHYSICAL_FILTER",
                })
                continue

            res = _physical_matching_cost(
                observed_system, modeled_system, alpha, beta, gamma, delta, method
            )

            if method=="probabilistic":
                c,b=res
            else:
                c=res; b={}

            costs[obs_index, mod_index]=c

            candidate_costs.append({
                "observed_index":int(obs_index),
                "modeled_index":int(mod_index),
                "observed_partition":observed_system.get("partition",obs_index),
                "modeled_partition":modeled_system.get("partition",mod_index),
                "cost":float(c),
                "probability":float(np.exp(-c)) if method=="probabilistic" else np.nan,
                "direction_difference_deg":dtheta,
                "tp_difference":dtp,
                "filter_reason":"OK",
                **b
            })

    if costs.size == 0:
        return {
            "matched_pairs": [],
            "unmatched_observed": observed,
            "unmatched_modeled": modeled,
            "cost_matrix": costs,
        }

    if method == "probabilistic":
        finite = costs[np.isfinite(costs) & (costs < 1.0e11)]
        if finite.size:
            dummy_cost = max(
                -np.log(max(probability_threshold, np.finfo(float).eps)),
                np.percentile(finite,95)
            )
        else:
            dummy_cost = -np.log(max(probability_threshold, np.finfo(float).eps))

    else:
        finite = costs[np.isfinite(costs)]
        dummy_cost = np.percentile(finite, 75) if finite.size else 1.0

    n_obs, n_mod = costs.shape
    N = max(n_obs, n_mod)
    aug = np.full((N+n_obs, N+n_mod), dummy_cost, dtype=float)
    aug[:n_obs, :n_mod] = costs

    for i in range(n_obs, N+n_obs):
        aug[i, n_mod:] = 0.0

    rows, cols = linear_sum_assignment(aug)

    matched_pairs=[]
    matched_observed=set()
    matched_modeled=set()

    for r,c in zip(rows, cols):
        if r < n_obs and c < n_mod and costs[r,c] < dummy_cost:

            if method == "probabilistic":
                _, breakdown = _physical_matching_cost(
                    observed[r], modeled[c],
                    alpha, beta, gamma, delta, method
                )
            else:
                breakdown = {}

            obs_dir = observed[r]["dp"]
            mod_dir = modeled[c]["dp"]

            dtheta = abs(
                np.degrees(
                    np.arctan2(
                        np.sin(np.radians(obs_dir - mod_dir)),
                        np.cos(np.radians(obs_dir - mod_dir))
                    )
                )
            )

            # ---- Diagnostic quantities (Stage 1) ----
            obs_tp = observed[r]["tp"]
            mod_tp = modeled[c]["tp"]
            tp_tol = get_tp_tolerance(
                obs_tp,
                tp_break_1=tp_break_1,
                tp_break_2=tp_break_2,
                tp_rel_low=tp_rel_low,
                tp_rel_mid=tp_rel_mid,
                tp_rel_high=tp_rel_high,
            )
            dtp = abs(obs_tp - mod_tp)
            direction_ok = dtheta <= max_direction_deg
            tp_ok = dtp <= tp_tol

            if physical_filter:
                print(f"[PCSPM] physical_filter={physical_filter} dtheta={dtheta:.1f}")

                if dtheta > max_direction_deg:
                    print(
                        f"[PCSPM] REJECTED match: "
                        f"obs={r} mod={c} dtheta={dtheta:.1f}"
                    )
                    for cc in candidate_costs:
                        if cc["observed_index"]==r and cc["modeled_index"]==c:
                            cc["filter_reason"]="Direction > max_direction_deg"
                    continue

                # Adaptive physical filter for peak period
                obs_tp = observed[r]["tp"]
                mod_tp = modeled[c]["tp"]

                tp_tol = get_tp_tolerance(
                    obs_tp,
                    tp_break_1=tp_break_1,
                    tp_break_2=tp_break_2,
                    tp_rel_low=tp_rel_low,
                    tp_rel_mid=tp_rel_mid,
                    tp_rel_high=tp_rel_high,
                )

                dtp = abs(obs_tp - mod_tp)

                # print(
                #     f"[PCSPM] Tp filter: "
                #     f"Obs={obs_tp:.2f}s Mod={mod_tp:.2f}s "
                #     f"ΔTp={dtp:.2f}s Tol={tp_tol:.2f}s"
                # )

                if dtp > tp_tol:
                    print(
                        f"[PCSPM] REJECTED match: "
                        f"obs={r} mod={c} ΔTp={dtp:.2f}s > {tp_tol:.2f}s"
                    )
                    for cc in candidate_costs:
                        if cc["observed_index"]==r and cc["modeled_index"]==c:
                            cc["filter_reason"]="Tp tolerance"
                    continue

            matched_pairs.append({
                "observed": observed[r],
                "modeled": modeled[c],
                "observed_index": int(r),
                "modeled_index": int(c),
                "cost": float(costs[r,c]),
                "direction_difference_deg": float(dtheta),
                "tp_difference": float(dtp),
                "tp_tolerance": float(tp_tol),
                "direction_ok": bool(direction_ok),
                "tp_ok": bool(tp_ok),
                "physical_filter_passed": bool(direction_ok and tp_ok),
                **breakdown,
                "match_probability": float(np.exp(-costs[r,c])) if method=="probabilistic" else None,
            })
            matched_observed.add(r)
            matched_modeled.add(c)

    for i in range(n_obs):
        s=sorted([x for x in candidate_costs if x["observed_index"]==i],key=lambda x:x["cost"])
        if len(s)>=2:
            c1,c2=s[0],s[1]
            matching_ranking.append({
                "observed_index":i,
                "observed_partition":c1["observed_partition"],
                "best_partition":c1["modeled_partition"],
                "second_partition":c2["modeled_partition"],
                "best_cost":c1["cost"],
                "second_cost":c2["cost"],
                "delta_cost":c2["cost"]-c1["cost"],
                "confidence":(1-c1["cost"]/c2["cost"]) if c2["cost"]>0 else np.nan,
                "tp_cost":c1.get("tp_cost",np.nan),
                "direction_cost":c1.get("direction_cost",np.nan),
                "energy_cost":c1.get("energy_cost",np.nan),
                "bandwidth_cost":c1.get("bandwidth_cost",np.nan),
                "spread_cost":c1.get("spread_cost",np.nan),
            })
    matching_ranking=sorted(matching_ranking,key=lambda x:x["delta_cost"])
    return {
        "candidate_costs":candidate_costs,
        "matching_ranking":matching_ranking,
        "matched_pairs": matched_pairs,
        "unmatched_observed": [
            system for index, system in enumerate(observed) if index not in matched_observed
        ],
        "unmatched_modeled": [
            system for index, system in enumerate(modeled) if index not in matched_modeled
        ],
        "cost_matrix": costs,
    }


def match_spectral_partitions(observed_descriptors, modeled_descriptors,
                              alpha=1.0, beta=1.0, gamma=1.0, delta=1.0,
                              method='weighted'):
    """Compatibility interface for matching precomputed descriptors.

    Despite its historical name, this function does not compute descriptors
    from spectra. It delegates the two descriptor sequences directly to
    :func:`match_partition_properties` using legacy ``method='weighted'`` by
    default.
    """
    return match_partition_properties(
        observed_descriptors, modeled_descriptors,
        alpha=alpha, beta=beta, gamma=gamma, delta=delta, method=method,
    )


def haversine_distance(lon1, lat1, lon2, lat2):
    """
    Calculate great circle distance between two points in km.
    
    Parameters:
    -----------
    lon1, lat1 : float or array
        Longitude and latitude of first point(s) in degrees
    lon2, lat2 : float or array
        Longitude and latitude of second point(s) in degrees
    
    Returns:
    --------
    float or array
        Distance in kilometers
    """
    R = 6371  # Earth radius in km
    
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    return R * c


def check_ndbc_has_spectral_data(ndbc_file):
    """
    Check if NDBC file contains spectral directional data.
    
    Parameters:
    -----------
    ndbc_file : str or Path
        Path to NDBC NetCDF file
    
    Returns:
    --------
    bool
        True if file has required spectral variables
    """
    required_vars = [
        'spectral_wave_density',
        'wave_spectrum_r1',
        'wave_spectrum_r2',
        'mean_wave_dir',
        'principal_wave_dir',
        'frequency'
    ]
    
    try:
        with xr.open_dataset(ndbc_file) as ds:
            has_all = all(var in ds.variables for var in required_vars)
            return has_all
    except:
        return False


def scan_ndbc_stations(ndbc_base_dir, year_range=None):
    """
    Scan NDBC directory and identify stations with spectral data.
    
    Parameters:
    -----------
    ndbc_base_dir : str or Path
        Base directory containing NDBC station folders
    year_range : tuple of int, optional
        (start_year, end_year) to filter files
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: station_id, year, file_path, has_spectral
    """
    ndbc_path = Path(ndbc_base_dir)
    results = []
    
    # List all station directories
    station_dirs = sorted([d for d in ndbc_path.iterdir() if d.is_dir()])
    
    print(f"Scanning {len(station_dirs)} NDBC stations...")
    
    for station_dir in station_dirs:
        station_id = station_dir.name
        
        # Find all NetCDF files for this station
        nc_files = sorted(station_dir.glob('*.nc'))
        
        for nc_file in nc_files:
            # Extract year from filename (e.g., 41010w2020.nc)
            try:
                year = int(nc_file.stem[-4:])
            except:
                continue
            
            # Filter by year range if specified
            if year_range and (year < year_range[0] or year > year_range[1]):
                continue
            
            # Check if has spectral data
            has_spectral = check_ndbc_has_spectral_data(nc_file)
            
            results.append({
                'station_id': station_id,
                'year': year,
                'file_path': str(nc_file),
                'has_spectral': has_spectral
            })
    
    df = pd.DataFrame(results)
    
    # Get station coordinates from first file of each station
    coords = []
    for station_id in df['station_id'].unique():
        first_file = df[df['station_id'] == station_id].iloc[0]['file_path']
        try:
            with xr.open_dataset(first_file) as ds:
                lon = float(ds.longitude.values.item())
                lat = float(ds.latitude.values.item())
                coords.append({'station_id': station_id, 'lon': lon, 'lat': lat})
        except:
            coords.append({'station_id': station_id, 'lon': np.nan, 'lat': np.nan})
    
    df_coords = pd.DataFrame(coords)
    df = df.merge(df_coords, on='station_id', how='left')
    
    return df


def find_sar_ndbc_matches(sar_dir, ndbc_info_df, max_distance_km=50, 
                          max_time_diff_hours=3, year_range=None, limit_files=None):
    """
    Find SAR observations within distance and time of NDBC stations.
    
    Parameters:
    -----------
    sar_dir : str or Path
        Directory containing SAR NetCDF files
    ndbc_info_df : pd.DataFrame
        DataFrame from scan_ndbc_stations() with spectral stations
    max_distance_km : float
        Maximum distance for a match (km)
    max_time_diff_hours : float
        Maximum time difference for a match (hours)
    year_range : tuple of int, optional
        (start_year, end_year) to filter SAR files
    limit_files : int, optional
        Limit processing to first N files (for testing)
    
    Returns:
    --------
    pd.DataFrame
        Matches with columns: station_id, sar_file, sar_index, sar_lon, sar_lat,
        sar_time, ndbc_lon, ndbc_lat, distance_km, ndbc_file
    """
    sar_path = Path(sar_dir)
    
    # Filter only stations with spectral data
    spectral_stations = ndbc_info_df[ndbc_info_df['has_spectral']].copy()
    
    if len(spectral_stations) == 0:
        print("⚠️  No NDBC stations with spectral data found!")
        return pd.DataFrame()
    
    # Get unique stations with coordinates
    stations = spectral_stations[['station_id', 'lon', 'lat']].drop_duplicates()
    stations = stations.dropna(subset=['lon', 'lat'])
    
    print(f"\nSearching matches for {len(stations)} stations with spectral data")
    print(f"Max distance: {max_distance_km} km")
    print(f"Max time diff: {max_time_diff_hours} hours")
    
    # Get all SAR files
    sar_files = sorted(sar_path.glob('*.nc'))
    
    # Filter by year if specified
    if year_range:
        sar_files = [f for f in sar_files 
                     if year_range[0] <= int(f.stem.split('_')[1][:4]) <= year_range[1]]
    
    # Limit files if requested (for testing)
    if limit_files:
        sar_files = sar_files[:limit_files]
    
    print(f"Processing {len(sar_files)} SAR files...")
    
    matches = []
    
    for sar_idx, sar_file in enumerate(sar_files):
        if (sar_idx + 1) % 100 == 0:
            print(f"  Processed {sar_idx + 1}/{len(sar_files)} files, found {len(matches)} matches")
        
        try:
            with xr.open_dataset(sar_file) as ds_sar:
                sar_lons = ds_sar['longitude'].values
                sar_lats = ds_sar['latitude'].values
                sar_times = pd.to_datetime(ds_sar['time'].values)
                
                # Handle both 1D and 2D coordinate arrays
                # SAR files can have shape (nobs,) or (nobs, ndir)
                if sar_lons.ndim == 2:
                    # Take first direction column for coordinates
                    sar_lons = sar_lons[:, 0]
                    sar_lats = sar_lats[:, 0]
                
                # Check each SAR observation
                for obs_idx in range(len(sar_lons)):
                    sar_lon = sar_lons[obs_idx]
                    sar_lat = sar_lats[obs_idx]
                    sar_time = sar_times[obs_idx]
                    
                    # Skip invalid coordinates
                    if np.isnan(sar_lon) or np.isnan(sar_lat):
                        continue
                    
                    # Calculate distance to all stations
                    distances = haversine_distance(
                        sar_lon, sar_lat,
                        stations['lon'].values, stations['lat'].values
                    )
                    
                    # Find stations within max distance
                    close_stations = stations[distances <= max_distance_km].copy()
                    close_stations['distance_km'] = distances[distances <= max_distance_km]
                    
                    if len(close_stations) == 0:
                        continue
                    
                    # For each close station, check if we have data at this time
                    for _, station in close_stations.iterrows():
                        station_id = station['station_id']
                        year = sar_time.year
                        
                        # Find NDBC file for this year
                        ndbc_files = spectral_stations[
                            (spectral_stations['station_id'] == station_id) &
                            (spectral_stations['year'] == year)
                        ]
                        
                        if len(ndbc_files) == 0:
                            continue
                        
                        ndbc_file = ndbc_files.iloc[0]['file_path']
                        
                        # Check if NDBC has data within time window
                        try:
                            with xr.open_dataset(ndbc_file) as ds_ndbc:
                                ndbc_times = pd.to_datetime(ds_ndbc.time.values)
                                time_diffs = np.abs(ndbc_times - sar_time)
                                min_time_diff = time_diffs.min()
                                
                                if min_time_diff <= pd.Timedelta(hours=max_time_diff_hours):
                                    closest_idx = np.argmin(time_diffs)
                                    ndbc_time = ndbc_times[closest_idx]
                                    
                                    matches.append({
                                        'station_id': station_id,
                                        'station_lon': station['lon'],
                                        'station_lat': station['lat'],
                                        'sar_file': sar_file.name,
                                        'sar_index': obs_idx,
                                        'sar_lon': sar_lon,
                                        'sar_lat': sar_lat,
                                        'sar_time': sar_time,
                                        'ndbc_file': Path(ndbc_file).name,
                                        'ndbc_time': ndbc_time,
                                        'ndbc_time_index': closest_idx,
                                        'distance_km': station['distance_km'],
                                        'time_diff_hours': min_time_diff.total_seconds() / 3600
                                    })
                        except Exception as e:
                            continue
                            
        except Exception as e:
            print(f"  ⚠ Error processing {sar_file.name}: {e}")
            continue
    
    print(f"\n✓ Found {len(matches)} total matches")
    
    return pd.DataFrame(matches)


def add_ww3_info(matches_df, ww3_dir):
    """
    Add WW3 file information to matches DataFrame.
    
    Parameters:
    -----------
    matches_df : pd.DataFrame
        DataFrame with SAR-NDBC matches
    ww3_dir : str or Path
        Directory containing WW3 NetCDF files (format: ww3_STATION.nc)
    
    Returns:
    --------
    pd.DataFrame
        Updated DataFrame with ww3_file and ww3_available columns
    """
    ww3_path = Path(ww3_dir)
    
    # Create copy to avoid modifying original
    df = matches_df.copy()
    
    # Add WW3 columns
    df['ww3_file'] = None
    df['ww3_available'] = False
    
    # Get list of available WW3 files
    ww3_files = {f.stem.split('_')[1]: f.name for f in ww3_path.glob('ww3_*.nc')}
    
    print(f"\nChecking WW3 data availability...")
    print(f"Found WW3 files for {len(ww3_files)} stations")
    
    # Check each match
    for idx, row in df.iterrows():
        station_id = str(row['station_id'])
        
        if station_id in ww3_files:
            ww3_file = ww3_files[station_id]
            
            # Verify WW3 has data at this time
            try:
                ww3_path_full = ww3_path / ww3_file
                with xr.open_dataset(ww3_path_full) as ds_ww3:
                    ww3_times = pd.to_datetime(ds_ww3.time.values)
                    sar_time = pd.to_datetime(row['sar_time'])
                    
                    # Check if time exists in WW3 data
                    time_diffs = np.abs(ww3_times - sar_time)
                    min_diff = time_diffs.min()
                    
                    if min_diff <= pd.Timedelta(hours=3):
                        df.at[idx, 'ww3_file'] = ww3_file
                        df.at[idx, 'ww3_available'] = True
            except Exception as e:
                continue
    
    n_with_ww3 = df['ww3_available'].sum()
    print(f"✓ {n_with_ww3}/{len(df)} matches have WW3 data available")
    
    return df
