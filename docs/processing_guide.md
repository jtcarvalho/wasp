# WASP processing and I/O guide

This guide describes the current package adapters and the handoff to
`partition_spectrum()`. It does not imply that all products share one direction
convention automatically.

## Core input contract

The scientific core expects:

```text
E2d.shape == (len(frequencies), len(directions_rad))
```

- frequency is in Hz;
- direction is in radians on a uniform full-circle grid;
- intended spectral-density units are `m² s rad⁻¹`;
- frequency is non-periodic and direction is periodic;
- the direction convention is preserved, not converted, by the core.

Callers comparing sensors must explicitly verify coming-from versus going-to
conventions.

## Sentinel/CMEMS SAR

Module: `wasp.io_sar`

`load_sar_spectrum(ds, date_time=None, index=0)` accepts an open Xarray dataset.
It first searches these preprocessed variable families:

- spectrum: `wave_spec`, `obs_params/wave_spec`, `wave_spectrum`, or
  `obs_params/wave_spectrum`;
- wavenumber: `wavenumber_spec` or `obs_params/wavenumber_spec`;
- direction: `direction_spec` or `obs_params/direction_spec`;
- several common time names.

If that read fails and `oswPolSpec` exists, it falls back to the legacy ESA
variables `oswPolSpec`, `oswK`, and `oswPhi`.

`convert_sar_energy_units()` applies the deep-water dispersion relation and
`dk/df = 8π²f/g`, plus the implemented angular-density factor. It returns:

```text
E2d, frequency, directions_deg, directions_rad, selected_time
```

The output spectrum is `(frequency, direction)`. The adapter treats SAR
directions as going-to and performs no 180° rotation.

## NDBC

Module: `wasp.io_ndbc`

`load_ndbc_spectrum(ds, time_index, direction_resolution=15)` reconstructs a
directional spectrum from:

- `spectral_wave_density`;
- `wave_spectrum_r1` and `wave_spectrum_r2`;
- `mean_wave_dir` and `principal_wave_dir`;
- `frequency`.

The implemented reconstruction is a truncated first/second-order Fourier
series, not a maximum-entropy solver. Negative reconstructed density is clipped
to zero after a diagnostic integration.

The returned tuple is:

```text
E2d, frequency, directions_deg, directions_rad, longitude, latitude
```

Directions are created from 0° through 360° in the NDBC coming-from convention.

`load_ndbc_at_time()` additionally finds a station/year file, selects the
nearest record, applies `max_time_diff_hours` (default 3 h), and returns a
metadata dictionary or `None`.

## WaveWatch III

Module: `wasp.io_ww3`

`load_ww3_spectrum(file_path, time_index)` reads:

```text
efth[time_index, 0, frequency, direction]
```

along with station coordinates and optional `wnd`/`wnddir`. The adapter adds
180° to the stored WW3 propagation direction and therefore returns its
direction grid in coming-from convention:

```text
E2d, frequency, directions_deg, directions_rad,
longitude, latitude, wind_speed, wind_direction
```

`find_closest_time()` opens the file independently and returns the closest time
index, timestamp, and absolute offset in hours. It does not apply an acceptance
threshold.

## CFOSAT SWIM

Module: `wasp.io_cfosat`

`load_cfosat_spectrum()` accepts a NetCDF path, box, orbit side (`posneg`), and
optional beam. It supports `p_combined` and `pp_mean` structures.

Current processing is:

1. read wavenumber, direction, position, time, and slope spectrum;
2. expand 12 directions to a mirrored 24-direction grid when needed;
3. convert wavenumber to frequency/wavelength using deep-water dispersion;
4. optionally mask wavelengths greater than `min_wavelength` (default 500 m);
5. convert slope density with `dk/df / k²`;
6. roll the directional spectrum by 180°;
7. read available processor Hs, peak wavelength/Tp, and direction metadata;
8. optionally normalize amplitude to processor Hs.

The return value is a dictionary. Its `spectrum` has shape
`(direction, frequency)`, unlike the partitioning core, so examples transpose
it before use. Other fields include `frequency`, `wavelength`, `k`, `direction`,
`units`, location, time, and `wave_params`.

When Hs normalization is enabled, the loader targets twice the processor m0.
This reflects the current example workflow, which subsequently removes one of
the two mirrored lobes with a notebook-local ambiguity function. The combined
rotation/ambiguity policy requires scientific review before it can be safely
centralized.

## Partitioning

```python
from wasp import partition_spectrum

result = partition_spectrum(
    E2d,
    frequencies,
    directions_rad,
    threshold_mode="adaptive",
    threshold_percentile=98.0,
    max_partitions=3,
    merge_factor=0.315,
)
```

The actual pipeline is:

```text
threshold → peaks → watershed → primary merge → SPR
          → energy renumbering → final parameters/moments/descriptors
```

The complete current behavior is summarized in `README_archi.md`; SPR is
specified in `WASP_SPR_specification.md`.

## Configuration

`examples/config.yaml` stores sensor-specific notebook settings. It is loaded
with `wasp.load_config()` but is not automatically applied to any reader or to
`partition_spectrum()`.

Current partition settings in the example file are:

| Source | Percentile | Merge factor | Max primary peaks |
|---|---:|---:|---:|
| WW3 | 98 | 0.5 | 3 |
| SAR | 98 | 0.3 | 3 |
| NDBC | 95 | 0.7 | 3 |
| CFOSAT | 95 | 0.3 | 3 |

`min_energy_fraction=0.01` is applied by example CSV-export cells, not by the
partitioning function itself. CFOSAT ambiguity and plotting keys are likewise
consumed by notebooks rather than injected by package configuration.

## Plotting

The top-level `wasp.plot_directional_spectrum()` accepts an `(NF, ND)` spectrum,
frequency in Hz, and direction in degrees. Radius is peak period in seconds;
angle uses north at zero and clockwise rotation. Optional total or per-system
parameters are shown in a side panel.

`wasp.plotting_geo.plot_directional_spectrum()` has a similar signature but
interpolates and smooths the displayed field. It is not exported at package top
level and does not alter partitioning data.
