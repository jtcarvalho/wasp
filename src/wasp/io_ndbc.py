"""
Functions for reading and processing NDBC buoy spectral data
"""

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path


def find_station_file(ndbc_data_dir, station_id, year):
    """
    Find NDBC file for given station and year.
    
    Parameters:
    -----------
    ndbc_data_dir : str or Path
        Directory containing NDBC data files
    station_id : str or int or float
        NDBC station identifier
    year : int
        Year to search for
    
    Returns:
    --------
    Path or None
        Path to NDBC file if found, None otherwise
    """
    station_str = str(int(float(station_id)))
    
    patterns = [
        f"{station_str}w{year}.nc",
        f"{station_str}w{year}.nc4",
        f"{station_str}_{year}.nc",
    ]
    
    ndbc_path = Path(ndbc_data_dir)
    for pattern in patterns:
        file_path = ndbc_path / pattern
        if file_path.exists():
            return file_path
    
    return None


def find_closest_time(ds, target_time_dt):
    """
    Find the closest time index in NDBC dataset to target time.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        NDBC dataset
    target_time_dt : pd.Timestamp
        Target time to search for
    
    Returns:
    --------
    tuple
        (time_index, closest_time, time_diff_hours) or (None, None, None) if error
    """
    try:
        ndbc_times = pd.to_datetime(ds.time.values)
        
        time_diffs = np.abs(ndbc_times - target_time_dt)
        itime = np.argmin(time_diffs)
        closest_time = ndbc_times[itime]
        time_diff_hours = time_diffs[itime].total_seconds() / 3600
        
        return itime, closest_time, time_diff_hours
    except Exception as e:
        print(f"  ⚠ Error finding closest time: {e}")
        return None, None, None


def load_ndbc_spectrum(ds, time_index, direction_resolution=15):
    """
    Load NDBC spectral data and reconstruct 2D spectrum in m²·s·rad⁻¹.
    
    NDBC provides:
    - spectral_wave_density: 1D spectrum in m²/Hz
    - Fourier coefficients (r1, r2, alpha1, alpha2) for directional reconstruction
    
    This function reconstructs the 2D directional spectrum using the 
    Maximum Entropy Method (MEM) based on Fourier coefficients.
    
    UNITS:
    - NDBC spectral_wave_density: m²/Hz
    - Spreading function D(θ): rad⁻¹ (normalized: ∫D(θ)dθ = 1 over 0 to 2π)
    - E2D(f,θ) = E(f) × D(θ) in m²/Hz/rad = m²·s·rad⁻¹
    
    Parameters:
    -----------
    ds : xarray.Dataset
        NDBC dataset
    time_index : int
        Time index to extract
    direction_resolution : float, optional
        Direction resolution in degrees (default: 15°)
    
    Returns:
    --------
    tuple or None
        (E2d, freq, dirs, dirs_rad, lon, lat) where:
        - E2d : ndarray (NF, ND) - Spectrum in m²·s·rad⁻¹
        - freq : ndarray (NF,) - Frequencies in Hz
        - dirs : ndarray (ND,) - Directions in degrees (oceanographic)
        - dirs_rad : ndarray (ND,) - Directions in radians
        - lon : float - Station longitude
        - lat : float - Station latitude
        Returns None if error occurs
    """
    try:
        # Get data at specific time index
        spec_1d = ds['spectral_wave_density'].isel(time=time_index).squeeze()
        r1 = ds['wave_spectrum_r1'].isel(time=time_index).squeeze()
        r2 = ds['wave_spectrum_r2'].isel(time=time_index).squeeze()
        alpha1 = ds['mean_wave_dir'].isel(time=time_index).squeeze()
        alpha2 = ds['principal_wave_dir'].isel(time=time_index).squeeze()
        
        # Get frequencies
        freq = ds['frequency'].values
        
        # Create directional grid
        # NDBC uses oceanographic FROM convention (direction waves are coming FROM),
        # measured clockwise from North. theta_deg matches this convention directly.
        theta_deg = np.arange(0, 360, direction_resolution)
        theta_rad = np.radians(theta_deg)
        ntheta = len(theta_rad)
        
        # Initialize directional matrix
        E2d = np.zeros((len(freq), ntheta))
        
        # Reconstruct 2D spectrum using Fourier coefficients
        # Based on NDBC Technical Document: 
        # D(θ) = (1/2π) * [1 + 2*r1*cos(θ-α1) + 2*r2*cos(2*(θ-α2))]
        for i, f in enumerate(freq):
            # Convert directional moments to Cartesian components
            a1 = r1.values[i] * np.cos(np.deg2rad(alpha1.values[i]))
            b1 = r1.values[i] * np.sin(np.deg2rad(alpha1.values[i]))
            a2 = r2.values[i] * np.cos(2 * np.deg2rad(alpha2.values[i]))
            b2 = r2.values[i] * np.sin(2 * np.deg2rad(alpha2.values[i]))
            
            # Directional spreading function D(θ) [rad⁻¹]
            # Normalized: ∫D(θ)dθ = 1 over θ ∈ [0, 2π]
            spread = (1 / (2 * np.pi)) * (
                1
                + 2 * a1 * np.cos(theta_rad)
                + 2 * b1 * np.sin(theta_rad)
                + 2 * a2 * np.cos(2 * theta_rad)
                + 2 * b2 * np.sin(2 * theta_rad)
            )
            
            # 2D spectrum: E(f,θ) = E(f) × D(θ)
            # Units: [m²/Hz] × [rad⁻¹] = [m²/(Hz·rad)] = [m²·s·rad⁻¹]
            E2d[i, :] = spec_1d.values[i] * spread
        
        # DEBUG: Check integrated energy
        _trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
        dtheta = np.deg2rad(direction_resolution)
        m0_check = 0.0
        for j in range(ntheta):
            E_clean = np.where(np.isfinite(E2d[:, j]) & (E2d[:, j] >= 0), E2d[:, j], 0)
            m0_check += _trapz(E_clean, freq) * dtheta
        
        print(f"  [DEBUG] NDBC m0={m0_check:.6f} m², Hs={4*np.sqrt(m0_check):.2f} m")
        
        # Ensure no negative values (can occur due to numerical issues)
        E2d = np.maximum(E2d, 0.0)
        
        # Get station coordinates (extract scalar from array)
        lon = float(ds.longitude.values.item())
        lat = float(ds.latitude.values.item())
        
        return E2d, freq, theta_deg, theta_rad, lon, lat
        
    except Exception as e:
        print(f"  ⚠ Error loading spectrum: {e}")
        return None


def load_ndbc_at_time(ndbc_data_dir, station_id, target_time_dt, max_time_diff_hours=3.0):
    """
    Load NDBC spectrum at specific time (convenience function).
    
    This function combines finding the file, opening the dataset, and loading
    the spectrum at the closest time to the target.
    
    Parameters:
    -----------
    ndbc_data_dir : str or Path
        Directory containing NDBC data files
    station_id : str or int or float
        NDBC station identifier
    target_time_dt : pd.Timestamp
        Target time to search for
    max_time_diff_hours : float, optional
        Maximum acceptable time difference in hours (default: 3.0)
    
    Returns:
    --------
    dict or None
        Dictionary with keys:
        - 'E2d' : 2D spectrum
        - 'freq' : frequencies
        - 'dirs' : directions in degrees
        - 'dirs_rad' : directions in radians
        - 'lon' : longitude
        - 'lat' : latitude
        - 'time' : actual time of data
        - 'time_diff_hours' : time difference from target
        Returns None if file not found or time difference too large
    """
    year = target_time_dt.year
    
    # Find NDBC file
    ndbc_file = find_station_file(ndbc_data_dir, station_id, year)
    if ndbc_file is None:
        print(f"  ⚠ NDBC file not found for station {station_id}, year {year}")
        return None
    
    try:
        # Open dataset
        ds = xr.open_dataset(ndbc_file)
        
        # Find closest time
        itime, closest_time, time_diff_hours = find_closest_time(ds, target_time_dt)
        
        if itime is None or time_diff_hours > max_time_diff_hours:
            print(f"  ⚠ No suitable NDBC data within {max_time_diff_hours}h " 
                  f"(diff: {time_diff_hours:.2f}h)")
            ds.close()
            return None
        
        # Load spectrum
        spectrum_data = load_ndbc_spectrum(ds, itime)
        ds.close()
        
        if spectrum_data is None:
            return None
        
        E2d, freq, dirs, dirs_rad, lon, lat = spectrum_data
        
        return {
            'E2d': E2d,
            'freq': freq,
            'dirs': dirs,
            'dirs_rad': dirs_rad,
            'lon': lon,
            'lat': lat,
            'time': closest_time,
            'time_diff_hours': time_diff_hours,
            'source_file': ndbc_file.name
        }
        
    except Exception as e:
        print(f"  ✖ Error processing NDBC data: {str(e)}")
        return None
