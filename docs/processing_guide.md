# WASP Processing Guide

This document describes the spectral processing pipeline applied to each data source supported by WASP: **Sentinel-1 SAR**, **CFOSAT SWIM**, **NDBC buoys**, and **WaveWatch III (WW3)** model output. It covers how spectra are loaded, converted to a common physical unit, and partitioned into individual wave systems.

---

## Common Framework

All sources are ultimately converted to a **2D directional spectrum** $E(f, \theta)$ in units of **m²·s·rad⁻¹** before partitioning. The significant wave height is recovered via:

$$H_s = 4\sqrt{m_0}, \quad m_0 = \int_0^\infty \int_0^{2\pi} E(f,\theta)\, d\theta\, df$$

Directions follow the **oceanographic convention** (direction waves are *propagating toward*), measured clockwise from North.

---

## 1. Sentinel-1 SAR (`io_sar.py`)

### Data format

Sentinel-1 Level-2 OCN products (preprocessed via CMEMS) provide the wave spectrum as $E(k, \theta)$ in **m⁴** (wavenumber–direction space). The spectrum may also appear in older ESA format as `oswPolSpec`.

### Unit conversion

The deep-water dispersion relation $\omega^2 = gk$ (with $\omega = 2\pi f$) provides the Jacobian:

$$\left|\frac{dk}{df}\right| = \frac{8\pi^2 f}{g}$$

The conversion from wavenumber space to frequency–direction space is:

$$E(f,\theta)\,[\text{m}^2\cdot\text{s}\cdot\text{rad}^{-1}] = E(k,\theta)\,[\text{m}^4] \times \left|\frac{dk}{df}\right| \times \frac{\pi}{180}$$

The angular factor $\pi/180$ accounts for the direction bins being in degrees in the raw data.

### Direction convention

SAR spectra are already in the oceanographic ("going-to") convention. No additional rotation is applied.

### Partitioning parameters (recommended)

| Parameter | Value |
|---|---|
| `threshold_percentile` | 98 |
| `max_partitions` | 3 |
| `merge_factor` | 0.3 (conservative) |

---

## 2. CFOSAT SWIM (`io_cfosat.py`)

### Data format

CFOSAT SWIM provides a **slope spectrum** $S_\text{slope}(k, \phi)$ across 12 azimuthal directions, measured in 6 tilted beams (6°, 8°, 10°). The raw data format is either:
- **L2** — full orbit granule, may contain multiple beams
- **L2PBOX** — along-track box averages

### 180° ambiguity

The SWIM radar cannot distinguish waves coming from direction $\phi$ vs. $\phi + 180°$. This intrinsic limitation means the spectrum is symmetric and replicated every 180°. WASP handles this by:
1. Expanding the 12-direction grid to 24 directions (0–360°)
2. Mirroring the half-space spectrum — the energy is preserved, not doubled

### Unit conversion

The slope spectrum is converted to an elevation spectrum in two steps:

**Step 1 — Slope to elevation:**
$$S_\eta(k) = \frac{S_\text{slope}(k)}{k^2}$$

**Step 2 — Wavenumber to frequency with Jacobian:**
$$E(f,\theta)\,[\text{m}^2\cdot\text{s}\cdot\text{rad}^{-1}] = S_\eta(k) \times \left|\frac{dk}{df}\right| \times \frac{\pi}{180}$$

### Wavelength filter

Wavelengths > 500 m are masked by default (`apply_wavelength_limit=True`) because SWIM resolution degrades at very long wavelengths.

### Hs normalization

When `normalize_to_file_hs=True` (default), the computed spectrum is rescaled so that the integrated $H_s$ matches the value stored in the CFOSAT product `wave_params`. This corrects for empirical normalization differences between SWIM beam combinations.

### Partitioning parameters (recommended)

| Parameter | Value |
|---|---|
| `threshold_percentile` | 95 |
| `max_partitions` | 3 |
| `merge_factor` | 0.3 |
| `remove_mirrored` | True |
| `mirror_tolerance` | 15° |

---

## 3. NDBC Buoys (`io_ndbc.py`)

### Data format

NDBC NetCDF files provide:
- `spectral_wave_density`: 1D variance density spectrum $S(f)$ in **m²/Hz**
- Directional Fourier coefficients: `wave_spectrum_r1`, `wave_spectrum_r2`, `mean_wave_dir`, `principal_wave_dir`

### 2D spectrum reconstruction

The 2D directional spectrum is reconstructed from the Fourier coefficients using the standard truncated Fourier series expansion (NDBC Technical Document):

$$D(\theta) = \frac{1}{2\pi}\left[1 + 2r_1\cos(\theta - \alpha_1) + 2r_2\cos(2(\theta - \alpha_2))\right]$$

where $r_1, r_2$ are spectral spreading parameters and $\alpha_1, \alpha_2$ are mean and principal directions. The spreading function is normalized so that $\int_0^{2\pi} D(\theta)\,d\theta = 1$.

The 2D spectrum is then:
$$E(f,\theta)\,[\text{m}^2\cdot\text{s}\cdot\text{rad}^{-1}] = S(f) \times D(\theta)$$

### Time matching

`load_ndbc_at_time()` finds the closest timestamp in the yearly NetCDF file to a target time. A configurable `max_time_diff_hours` threshold (default: 3 h) guards against using temporally distant data.

### Partitioning parameters (recommended)

| Parameter | Value |
|---|---|
| `threshold_percentile` | 95 |
| `max_partitions` | 3 |
| `merge_factor` | 0.7 (more permissive, buoys are single-point) |

---

## 4. WaveWatch III Model (`io_ww3.py`)

### Data format

WW3 NetCDF output (`efth` variable) stores the directional spectrum in **m²·s·rad⁻¹** directly — no unit conversion is needed.

### Direction convention

WW3 stores propagation direction as *"going-to"* (oceanographic convention), so no correction is applied. The module reads `efth[time, station, frequency, direction]`.

### Co-location with observations

`find_closest_time()` identifies the WW3 timestamp nearest to a target observation time, enabling direct spectral comparison with SAR or buoy data.

### Partitioning parameters (recommended)

| Parameter | Value |
|---|---|
| `threshold_percentile` | 95 |
| `max_partitions` | 3 |
| `merge_factor` | 0.5 (moderate) |

---

## 5. Spectral Partitioning (`partition.py`)

The `partition_spectrum()` function implements the **Hanson & Phillips (2001)** watershed algorithm in seven steps:

| Step | Function | Description |
|---|---|---|
| 1 | `identify_spectral_peaks` | Find local maxima using 3×3 neighborhood analysis |
| 2 | `generate_mask` | Watershed propagation — assign each bin to the nearest peak |
| 3 | `calculate_peak_distances` | Squared Euclidean distances between peaks in $(f\cos\theta, f\sin\theta)$ space |
| 4 | `calculate_peak_spreading` | Variance of each partition in spectral Cartesian space |
| 5 | `merge_overlapping_systems` | Merge partitions if `dist(i,j) ≤ merge_factor × Eip(i)` |
| 6 | `calculate_peak_parameters` | Compute $T_p$ and $D_p$ for each surviving partition |
| 7 | `renumber_partitions_by_energy` | Sort partitions: partition 1 = highest $H_s$ |

### Energy threshold

The adaptive threshold (default mode) uses a configurable percentile of the spectrum energy to suppress noise peaks:

```python
results = partition_spectrum(
    E, freq, dirs_rad,
    threshold_mode='adaptive',
    threshold_percentile=95.0,
    max_partitions=3,
    merge_factor=0.5,
)
```

### Output

`partition_spectrum()` returns a dictionary with:
- `mask`: 2D array assigning each spectral bin to a partition
- `Hs`, `Tp`, `Dp`: arrays indexed by partition number (1-based)
- `energy`: spectral moment $m_0$ per partition
- `total_Hs`, `total_Tp`, `total_m0`: full-spectrum integrated values

---

## 6. Wave Parameters (`wave_params.py`)

The canonical wave parameter functions are:

### `calculate_wave_parameters(E2d, freq, dirs_rad)`

Returns $(H_s, T_p, D_p, m_0, \Delta f, \Delta\theta, i_\text{peak}, j_\text{peak})$ using trapezoidal integration:

$$m_0 = \sum_j \left(\int E(f,\theta_j)\, df\right) \Delta\theta$$

### `spectrum1d_from_2d(E2d, dirs_rad)`

Integrates the 2D spectrum over direction to produce $E(f)$:

$$E(f) = \sum_j E(f, \theta_j) \cdot \Delta\theta$$

### Unit conversions

`convert_spectrum_units()` handles conversions between:
- `m2_s_rad` ↔ `m2_Hz_rad` (factor: $2\pi$)
- `m2_Hz_rad` → `m2_Hz_deg` (factor: $180/\pi$)

---

## 7. Configuration (`config.yaml`)

All tunable parameters are controlled through `examples/config.yaml`:

```yaml
partitioning:
  sar:
    min_energy_fraction: 0.01
    max_partitions: 3
    threshold_percentile: 98
    merge_factor: 0.3
  ww3:
    threshold_percentile: 95.0
    merge_factor: 0.5
  ndbc:
    threshold_percentile: 95.0
    merge_factor: 0.7
  cfosat:
    threshold_percentile: 95.0
    merge_factor: 0.3
    remove_mirrored: true
    mirror_tolerance: 15
```

---

## 8. Visualization (`plotting.py`)

`plot_directional_spectrum()` renders the spectrum in polar coordinates with:
- Radial axis = wave period (s), up to 25 s
- Angular axis = propagation direction (oceanographic, clockwise from North)
- Optional overlay of partition wave parameters ($H_s$, $T_p$, $D_p$)

---

## References

- Hanson, J. L., & Phillips, O. M. (2001). *Automated analysis of ocean surface directional wave spectra.* Journal of Atmospheric and Oceanic Technology, 18(2), 277–293.
- NDBC Technical Document 96-01 — *Nondirectional and Directional Wave Data Analysis Procedures.*
- Hauser, D., et al. (2021). *CFOSAT, a new Chinese-French satellite for joint observations of ocean wind vector and directional spectra of ocean waves.* JGR Oceans.
