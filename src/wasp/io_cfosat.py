"""
Functions for reading and processing CFOSAT SWIM spectral data
"""

import numpy as np
import pandas as pd
import datetime
import netCDF4 as nc
from math import pi


def k_to_wavelength(k):
    """
    Convert wavenumber to wavelength.
    
    Parameters:
    -----------
    k : float or array
        Wavenumber in rad/m
        
    Returns:
    --------
    wavelength : float or array
        Wavelength in meters
    """
    return 2 * pi / k


def k_to_frequency(k, gravity=9.81):
    """
    Convert wavenumber to frequency using deep water dispersion relation.
    omega^2 = g*k
    
    Parameters:
    -----------
    k : float or array
        Wavenumber in rad/m
    gravity : float
        Gravitational acceleration (m/s^2)
        
    Returns:
    --------
    frequency : float or array
        Frequency in Hz
    """
    omega = np.sqrt(gravity * k)  # angular frequency (rad/s)
    return omega / (2 * pi)  # frequency in Hz


def convert_cfosat_slope_to_elevation(spectrum_slope, k_spectra, frequencies, 
                                     normalization_factor=1.0):
    """
    Convert CFOSAT slope spectrum to elevation spectrum in m²·s/rad.
    
    CFOSAT provides "mean slope spectrum" in meters (m).
    This function converts to ELEVATION SPECTRUM in m²·s/rad for comparison with SAR/WW3.
    
    The conversion applies:
    1. Slope to elevation conversion: S_η(k) = S_slope(k) / k²
    2. k to f conversion with Jacobian: S(f) = S(k) × |dk/df|
    3. Direction conversion: degrees to radians (π/180)
    
    Parameters:
    -----------
    spectrum_slope : ndarray [n_directions, n_frequencies]
        Slope spectrum from CFOSAT
    k_spectra : ndarray [n_frequencies]
        Wavenumber array (rad/m)
    frequencies : ndarray [n_frequencies]
        Frequency array (Hz)
    normalization_factor : float
        Empirical normalization factor (default: 1.0)
    
    Returns:
    --------
    spectrum_elevation : ndarray [n_directions, n_frequencies]
        Elevation spectrum in m²·s/rad
    """
    g = 9.81  # gravity
    spectrum_elevation = spectrum_slope.copy()
    
    # Apply conversion for each frequency
    for i_f in range(len(frequencies)):
        k = k_spectra[i_f]
        f = frequencies[i_f]
        
        if k > 0 and f > 0:
            # Jacobian: dk/df = 8π²f/g from dispersion ω² = gk
            dkdf = 8 * pi**2 * f / g
            # Combined conversion factor (WITHOUT deg_to_rad - that's in the integration step)
            conversion = normalization_factor * dkdf / (k**2)
            spectrum_elevation[:, i_f] *= conversion
    
    return spectrum_elevation


def load_cfosat_variables(filepath):
    """
    Load basic variables from CFOSAT NetCDF file.
    
    Parameters:
    -----------
    filepath : str
        Path to CFOSAT NetCDF L2 file
    
    Returns:
    --------
    dict with:
        'cdf': NetCDF Dataset object (remember to close!)
        'time': time array
        'lon': longitude array
        'lat': latitude array
        'k_spectra': wavenumber array (rad/m)
        'phi_vector': direction array (degrees)
        'pp_mean': mean slope spectrum
        'file_type': 'L2' or 'L2PBOX'
        'has_beams': boolean
        'n_k': number of wavenumber bins
        'n_phi': number of direction bins
        'n_boxes': number of boxes
    """
    cdf = nc.Dataset(filepath)
    
    # Load coordinate and spatial variables
    time = cdf.variables['time_spec_l2'][:]
    lon = cdf.variables['lon_spec_l2'][:]
    lat = cdf.variables['lat_spec_l2'][:]
    
    # Load wavenumber and direction vectors
    k_spectra = cdf.variables['k_spectra'][:]  # Wavenumber (rad/m)
    phi_vector = cdf.variables['phi_vector'][:]  # Direction (degrees)
    
    # Check file format (L2 or L2PBOX)
    if 'p_combined' in cdf.variables:
        # L2 file - use p_combined which already combines all beams
        pp_mean = cdf.variables['p_combined'][:]
        file_type = 'L2'
        has_beams = False
    elif 'pp_mean' in cdf.variables:
        pp_mean_var = cdf.variables['pp_mean']
        pp_mean = pp_mean_var[:]
        # Check if has beam dimension
        if len(pp_mean_var.shape) == 5:
            file_type = 'L2'
            has_beams = True
        else:
            file_type = 'L2PBOX'
            has_beams = False
    else:
        cdf.close()
        raise ValueError("File does not contain pp_mean or p_combined")
    
    # Dimensions
    n_k = len(k_spectra)  # 32
    n_phi = len(phi_vector)  # 12 or 24
    n_boxes = cdf.dimensions['n_box'].size
    
    return {
        'cdf': cdf,
        'time': time,
        'lon': lon,
        'lat': lat,
        'k_spectra': k_spectra,
        'phi_vector': phi_vector,
        'pp_mean': pp_mean,
        'file_type': file_type,
        'has_beams': has_beams,
        'n_k': n_k,
        'n_phi': n_phi,
        'n_boxes': n_boxes
    }


def load_cfosat_spectrum(filepath, box, posneg=0, beam_index=None, 
                        apply_wavelength_limit=True, min_wavelength=500,
                        normalization_factor=0.0267):
    """
    Load and process CFOSAT SWIM 2D directional spectrum.
    
    CFOSAT SWIM measures waves in 6 tilted beams (6°, 8°, 10°) on both sides
    of nadir. Each beam measures in 12 azimuthal directions.
    
    The 180° ambiguity means CFOSAT cannot distinguish if a wave comes from 
    direction φ or φ+180°. The spectrum appears "mirrored" - replicated at 180°.
    This is an intrinsic characteristic of the SWIM instrument.
    
    Parameters:
    -----------
    filepath : str
        Path to CFOSAT NetCDF L2 file
    box : int
        Box number (0 to n_boxes-1)
    posneg : int
        Side of orbit (0 or 1)
    beam_index : int, optional
        Beam index (0-2). If None, uses combined data if available.
        Only for L2 files with multiple beams
    apply_wavelength_limit : bool
        If True, masks wavelengths > min_wavelength (unreliable data)
    min_wavelength : float
        Wavelength limit in meters (default: 500m)
    normalization_factor : float
        Empirical normalization factor for elevation spectrum (default: 0.0267)
        This value was calibrated to match CFOSAT Hs measurements.
        Factor ≈ 1/37.5 determined empirically comparing m0 vs observed Hs.
        
    Returns:
    --------
    dict with:
        'frequency': frequency array (Hz)
        'wavelength': wavelength array (m)
        'k': wavenumber array (rad/m)
        'direction': direction array (degrees, 0-360)
        'spectrum': spectrum matrix [n_directions x n_frequencies] in m²·s/rad
        'units': spectrum units ('m²·s/rad' - elevation spectrum)
        'spectrum_type': 'elevation'
        'lon': measurement longitude
        'lat': measurement latitude
        'time': measurement time (datetime)
        'wave_params': dict with Hs, Tp, Dp if available
    """
    # Load variables
    data = load_cfosat_variables(filepath)
    cdf = data['cdf']
    
    # Extract raw spectrum based on file type
    if data['has_beams']:
        # L2 file with beams: [nk, n_phi, n_posneg, n_box, n_beam]
        if beam_index is None:
            beam_index = 0  # Use first beam as default
        spec_raw = data['pp_mean'][:, :, posneg, box, beam_index]
    else:
        # L2PBOX or L2 with p_combined: [nk, n_phi, n_posneg, n_box]
        spec_raw = data['pp_mean'][:, :, posneg, box]
    
    # STEP 1: Create direction vector
    # Check if file already has 360° or needs expansion
    if data['n_phi'] >= 24:
        # File already has complete 0-360° spectrum
        directions_full = data['phi_vector']
        # Transpose spectrum from [n_k, n_phi] to [n_phi, n_k]
        spectrum_360 = spec_raw.T
    else:
        # File has only 12 directions, needs to apply 180° ambiguity
        n_dir_full = data['n_phi'] * 2  # 24 directions
        directions_full = np.linspace(0, 360, n_dir_full, endpoint=False)
        
        # STEP 2: Reorganize spectrum to 360° using 180° ambiguity logic
        # spec_raw has shape [n_k, n_phi] = [32, 12]
        # Create spectrum with phi=0° at North (oceanographic convention)
        spectrum_360 = np.ma.zeros((n_dir_full, data['n_k']))
        
        # IMPORTANT: 180° ambiguity means we don't know if wave comes from φ or φ+180°
        # We replicate the spectrum but the ENERGY is the same - we just don't know exact direction
        # The mirroring represents uncertainty, not energy division
        
        # First half (0-180°)
        spectrum_360[0:6, :] = spec_raw[:, 5::-1].T
        spectrum_360[6:12, :] = spec_raw[:, 11:5:-1].T
        
        # Second half (180-360°) - replication due to ambiguity
        spectrum_360[12:18, :] = spec_raw[:, 5::-1].T
        spectrum_360[18:24, :] = spec_raw[:, 11:5:-1].T
    
    # STEP 3: Convert wavenumber to frequency and wavelength
    frequencies = k_to_frequency(data['k_spectra'])
    wavelengths = k_to_wavelength(data['k_spectra'])
    
    # STEP 4: Apply mask for invalid values
    spectrum_360 = np.ma.masked_where(spectrum_360 < -1e8, spectrum_360)
    
    # STEP 5: Apply wavelength limit (> 500m unreliable)
    if apply_wavelength_limit:
        k_limit = 2*pi / min_wavelength
        for i in range(spectrum_360.shape[0]):
            spectrum_360[i, :] = np.ma.masked_where(
                data['k_spectra'] < k_limit, 
                spectrum_360[i, :]
            )
    
    # STEP 6: Convert slope spectrum to elevation spectrum
    spectrum_elevation = convert_cfosat_slope_to_elevation(
        spectrum_360, data['k_spectra'], frequencies, normalization_factor
    )
    
    # STEP 7: Apply 180° rotation to spectrum to match direction convention
    # CFOSAT has 180° ambiguity - the Dp is corrected by adding 180°
    # We need to rotate the spectrum data by 180° to match this correction
    # so that the plotted spectrum aligns with the corrected Dp value
    n_dir = spectrum_elevation.shape[0]
    shift_amount = n_dir // 2  # Rotate by 180° (half of directions)
    spectrum_elevation = np.roll(spectrum_elevation, shift_amount, axis=0)
    
    # Get location and time information
    # CFOSAT reference: 2009-01-01 00:00:00 (units: s+us since 2009-01-01)
    t0 = datetime.datetime(2009, 1, 1, 0, 0, 0)
    time_val = data['time'][posneg, box]
    if np.ndim(time_val) > 0:
        time_seconds = float(time_val[0])
        time_microseconds = float(time_val[1]) if len(time_val) > 1 else 0
    else:
        time_seconds = float(time_val)
        time_microseconds = 0
    measurement_time = t0 + datetime.timedelta(seconds=time_seconds, 
                                               microseconds=time_microseconds)
    
    # Extract wave parameters calculated by CFOSAT (Hs, Tp, Dp)
    wave_params = {}
    try:
        # Try wave_param_combined first (combined data from all beams)
        # Dimensions: (nparam, n_posneg, n_box)
        # Parameters: [0]=SWH, [1]=peak_wavelength, [2]=peak_direction
        if 'wave_param_combined' in cdf.variables:
            wp = cdf.variables['wave_param_combined']
            wp_shape = wp.shape
            # Check dimensions: (nparam, n_posneg, n_box)
            if len(wp_shape) >= 3 and wp_shape[0] >= 3 and posneg < wp_shape[1] and box < wp_shape[2]:
                Hs = float(wp[0, posneg, box])  # Significant wave height
                peak_wavelength = float(wp[1, posneg, box])  # Peak wavelength (m)
                Dp = float(wp[2, posneg, box])  # Peak direction
                
                # Convert peak wavelength to period using dispersion relation
                # omega^2 = g*k, where k = 2*pi/wavelength
                if peak_wavelength > 0:
                    k_peak = 2 * pi / peak_wavelength
                    omega_peak = np.sqrt(9.81 * k_peak)
                    Tp = 2 * pi / omega_peak
                else:
                    Tp = np.nan
                
                # CFOSAT has 180° ambiguity - add 180° to get correct direction
                # The instrument cannot distinguish between φ and φ+180°
                Dp_corrected = (Dp + 180) % 360
                
                wave_params = {
                    'Hs': Hs if Hs > -999 else np.nan,
                    'Tp': Tp if Tp > -999 else np.nan,
                    'Dp': Dp_corrected if Dp > -999 else np.nan,
                    'Dp_raw': Dp if Dp > -999 else np.nan,
                    'peak_wavelength': peak_wavelength if peak_wavelength > -999 else np.nan
                }
        elif 'wave_param' in cdf.variables:
            # Fallback to wave_param (with beams)
            wp = cdf.variables['wave_param']
            wp_shape = wp.shape
            bi = beam_index if beam_index is not None else 0
            
            # Try to read with proper dimension order
            # Need to check actual dimensions
            if len(wp_shape) == 4:  # (nparam, n_posneg, n_box, n_beam)
                if wp_shape[0] >= 3 and posneg < wp_shape[1] and box < wp_shape[2] and bi < wp_shape[3]:
                    Hs = float(wp[0, posneg, box, bi])
                    peak_wavelength = float(wp[1, posneg, box, bi])
                    Dp = float(wp[2, posneg, box, bi])
                    
                    if peak_wavelength > 0:
                        k_peak = 2 * pi / peak_wavelength
                        omega_peak = np.sqrt(9.81 * k_peak)
                        Tp = 2 * pi / omega_peak
                    else:
                        Tp = np.nan
                    
                    # CFOSAT has 180° ambiguity - add 180° to get correct direction
                    Dp_corrected = (Dp + 180) % 360
                    
                    wave_params = {
                        'Hs': Hs if Hs > -999 else np.nan,
                        'Tp': Tp if Tp > -999 else np.nan,
                        'Dp': Dp_corrected if Dp > -999 else np.nan,
                        'Dp_raw': Dp if Dp > -999 else np.nan,
                        'peak_wavelength': peak_wavelength if peak_wavelength > -999 else np.nan
                    }
    except Exception as e:
        # If fails, leave empty
        wave_params = {}
    
    # Close NetCDF file
    cdf.close()
    
    result = {
        'frequency': frequencies,
        'wavelength': wavelengths,
        'k': data['k_spectra'],
        'direction': directions_full,
        'spectrum': spectrum_elevation,
        'units': 'm²·s/rad',
        'spectrum_type': 'elevation',
        'lon': float(data['lon'][posneg, box]),
        'lat': float(data['lat'][posneg, box]),
        'time': measurement_time,
        'wave_params': wave_params
    }
    
    return result


def find_cfosat_boxes_in_region(filepath, lat_min, lat_max, lon_min, lon_max):
    """
    Find all CFOSAT boxes within a geographic region.
    
    Parameters:
    -----------
    filepath : str
        Path to CFOSAT NetCDF L2 file
    lat_min, lat_max : float
        Latitude bounds (degrees)
    lon_min, lon_max : float
        Longitude bounds (degrees)
    
    Returns:
    --------
    list of dict
        Each dict contains: {'box': int, 'posneg': int, 'lat': float, 'lon': float}
    """
    data = load_cfosat_variables(filepath)
    
    boxes_in_region = []
    
    for posneg in range(2):
        for box in range(data['n_boxes']):
            try:
                lat = float(data['lat'][posneg, box])
                lon = float(data['lon'][posneg, box])
                
                if (lat_min <= lat <= lat_max and 
                    lon_min <= lon <= lon_max):
                    boxes_in_region.append({
                        'box': box,
                        'posneg': posneg,
                        'lat': lat,
                        'lon': lon
                    })
            except:
                continue
    
    data['cdf'].close()
    
    return boxes_in_region


def find_closest_cfosat_box(filepath, target_lat, target_lon):
    """
    Find CFOSAT box closest to target location.
    
    Parameters:
    -----------
    filepath : str
        Path to CFOSAT NetCDF L2 file
    target_lat : float
        Target latitude (degrees)
    target_lon : float
        Target longitude (degrees)
    
    Returns:
    --------
    dict with:
        'box': int
        'posneg': int
        'lat': float
        'lon': float
        'distance_km': float (approximate distance)
    """
    data = load_cfosat_variables(filepath)
    
    min_distance = float('inf')
    closest_box = None
    
    for posneg in range(2):
        for box in range(data['n_boxes']):
            try:
                lat = float(data['lat'][posneg, box])
                lon = float(data['lon'][posneg, box])
                
                # Simple distance calculation (approximate)
                dlat = lat - target_lat
                dlon = lon - target_lon
                # Rough conversion to km (1 degree ~ 111 km)
                distance = np.sqrt((dlat * 111)**2 + 
                                 (dlon * 111 * np.cos(np.radians(lat)))**2)
                
                if distance < min_distance:
                    min_distance = distance
                    closest_box = {
                        'box': box,
                        'posneg': posneg,
                        'lat': lat,
                        'lon': lon,
                        'distance_km': distance
                    }
            except:
                continue
    
    data['cdf'].close()
    
    return closest_box
