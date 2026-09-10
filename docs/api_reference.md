# Source-derived API inventory

This inventory records the current callable surface of `src/wasp/`. It is not a
promise to export every non-underscored helper forever. Exact signatures and
behavior come from source; detailed scientific rules live in the linked guides.
No classes are defined in the current package modules.

## Top-level facade

`import wasp` exposes:

- `partition_spectrum`
- `calculate_wave_parameters`
- `plot_directional_spectrum`
- `load_config`
- `__version__` (`"2.0.0"` in the current source tree)

Readers and PCSPM are imported from their modules, not from the facade.

## `wasp.partition`

| Callable | Current role |
|---|---|
| `partition_spectrum(E, frequencies, directions_rad, energy_threshold=None, max_partitions=3, threshold_mode='adaptive', threshold_percentile=98.0, merge_factor=0.315, spr_min_peak_prominence=0.4, spr_min_relative_energy=0.1, spr_diagnostic_filename=None)` | Full partition/SPR workflow; returns a result dictionary or `None` when no primary peak is usable |
| `identify_spectral_peaks(E, NF, ND, energy_threshold, max_partitions)` | Returns direction codes, seed mask, one-based primary peaks, primary counter, and one-based secondary peaks |
| `generate_mask(ICOD, MASK, NF, ND)` | Propagates/fills primary labels across the grid |
| `calculate_peak_distances(peaks, frequencies, directions_rad, nmask)` | Squared Cartesian peak-distance matrix |
| `calculate_peak_spreading(E, MASK, frequencies, directions_rad, NF, ND, nmask, Etot, delf, ddir)` | Primary `Eip` array |
| `merge_overlapping_systems(MASK, dist, Eip, peaks, nmask, merge_factor=0.5)` | Ordered primary-label merge; does not update `nmask` |
| `calculate_partitioned_energy(E, M, delf, ddir, NF, ND, nmask)` | Label-indexed energy and `Hs` arrays |
| `calculate_peak_parameters(E, mask, frequencies, directions_rad, NF, ND, nmask, delf, ddir)` | Label-indexed `Tp`/`Dp` arrays; empty labels are `NaN` |
| `secondary_peak_reassessment(..., merge_factor=None, alpha=0.05, ICOD=None, primary_peaks=None, min_peak_prominence=0.4, min_relative_energy=0.1)` | SPR; returns updated mask/energy/peak properties/counter/log |
| `renumber_partitions_by_energy(mask, Hs, e=None)` | Compacts existing positive labels by descending `Hs` |
| `calculate_spectral_moments(E, mask, freq, dirs_rad, delf, ddir, partition_idx=None)` | `(m0, m1, m2)` using angular-frequency moments |
| `build_partition_descriptors(...)` | Current `tp`/`dp` result descriptors |
| `build_descriptors_from_partition_results(...)` | Compatibility alias for the preceding builder |
| `calculate_directional_spreading(...)` | Circular width in degrees for a labelled partition |
| `calculate_partition_spectral_spreading(...)` | Dimensionless moment width |
| `plot_spr_diagnostic(..., filename=None)` | Returns a figure/axis and optionally writes PNG |
| `plot_directional_spectrum(E2d, freq, dirs, selected_time, hs, tp, dp)` | Legacy plotting function; not the top-level exported implementation |

See [partitioning](partitioning.md) and [SPR](WASP_SPR_specification.md).
The orchestrator's `moments` result field contains `total` as an
`(m0, m1, m2)` tuple plus label-indexed `m0`, `m1`, and `m2` arrays.
For direct SPR calls containing candidates, `merge_factor=None` has no fallback
and will fail when the compatibility calculation is reached; the public
orchestrator always supplies a numeric merge factor.

## `wasp.wave_params`

| Callable | Inputs and return |
|---|---|
| `calculate_wave_parameters(E2d, freq, dirs_rad)` | `(NF, ND)` spectrum; returns `(hs, tp, dp, m0, delf, ddir, i_peak, j_peak)` |
| `spectrum1d_from_2d(E2d, dirs_rad)` | Cleans invalid/negative values and returns direction-integrated spectrum plus `ddir` |
| `convert_meteorological_to_oceanographic(met_dir)` | Adds 180 degrees modulo 360 |
| `convert_spectrum_units(E2d, freq, dirs, from_unit, to_unit)` | Implements only named conversion branches; unsupported pairs return an unchanged copy |

The routines assume a uniform complete direction grid. `calculate_wave_parameters`
requires at least two frequency bins in practice because it repeats the last
forward frequency increment.

## `wasp.matching`

| Callable | Current role |
|---|---|
| `compute_partition_descriptors(E2d, frequencies, directions_rad, mask, partition_labels=None)` | Standalone descriptors using peak frequency/direction, Hz bandwidth, radian direction width, and Hz² Cartesian spread |
| `get_tp_tolerance(obs_tp, tp_break_1=10.0, tp_break_2=14.0, tp_rel_low=0.20, tp_rel_mid=0.25, tp_rel_high=0.30)` | Piecewise observed-period tolerance |
| `match_partition_properties(..., method='probabilistic', probability_threshold=0.05, physical_filter=False, max_direction_deg=60.0, ..., max_direction_cost=2.0)` | Augmented Hungarian PCSPM assignment and diagnostics |
| `match_spectral_partitions(..., method='weighted')` | Compatibility wrapper for already-computed descriptors |
| `haversine_distance(...)` | Great-circle distance in km |
| `check_ndbc_has_spectral_data(ndbc_file)` | Required-variable check; returns `False` on read errors |
| `scan_ndbc_stations(ndbc_base_dir, year_range=None)` | DataFrame inventory of station files and coordinates |
| `find_sar_ndbc_matches(..., max_distance_km=50, max_time_diff_hours=3, year_range=None, limit_files=None)` | Filesystem/time/space colocation DataFrame |
| `add_ww3_info(matches_df, ww3_dir)` | Adds WW3 availability within a hard-coded three-hour window |

See [PCSPM matching](matching.md) for required dictionaries, active equations,
assignment, ambiguity diagnostics, and empty-input behavior.

## Source adapters

| Module/callable | Input | Return and convention |
|---|---|---|
| `io_sar.load_sar_spectrum(ds, date_time=None, index=0)` | Open Xarray-like dataset | `(E2d, freq, dirs_deg, dirs_rad, time)`; `(NF, ND)`; going-to coordinate preserved |
| `io_sar.convert_sar_energy_units(E_sar, k, phi)` | Wavenumber spectrum, rad/m, degree directions | Converted `(NF, ND)` density, frequency, degrees, radians |
| `io_ndbc.load_ndbc_spectrum(ds, time_index, direction_resolution=15)` | Open Xarray dataset | Tuple `(E2d, freq, dirs_deg, dirs_rad, lon, lat)` or `None`; coming-from grid |
| `io_ndbc.load_ndbc_at_time(..., max_time_diff_hours=3.0)` | Directory/station/time | Metadata dictionary or `None` |
| `io_ndbc.find_station_file(...)` | Directory, numeric-like station identifier, year | Matching `Path` or `None` |
| `io_ndbc.find_closest_time(ds, target_time_dt)` | Open dataset/time | `(index, time, hours)` or three `None` values on error |
| `io_ww3.load_ww3_spectrum(file_path, time_index)` | NetCDF path/index | `(E2d, freq, dirs_deg, dirs_rad, lon, lat, wnd, wnddir)`; stored direction +180° |
| `io_ww3.find_closest_time(file_path, target_time_dt)` | NetCDF path/time | Nearest `(index, time, hours)` with no acceptance threshold |
| `io_cfosat.load_cfosat_spectrum(filepath, box, posneg=0, beam_index=None, apply_wavelength_limit=True, min_wavelength=500, normalize_to_file_hs=True)` | CFOSAT NetCDF selection | Dictionary; spectrum shape `(direction, frequency)` and therefore requires transpose for the core |
| `io_cfosat.load_cfosat_variables(filepath)` | NetCDF path | Metadata plus an open `cdf` object that the caller must close |
| `io_cfosat.find_cfosat_boxes_in_region(...)` | Path and bounds | List of box/side/location dictionaries |
| `io_cfosat.find_closest_cfosat_box(...)` | Path and target position | Closest box dictionary or `None` |
| `io_cfosat.k_to_wavelength(k)` / `k_to_frequency(k, gravity=9.81)` | Wavenumber | Deep-water wavelength/frequency |
| `io_cfosat.convert_cfosat_slope_to_elevation(...)` | Direction-by-frequency slope spectrum | Direction-by-frequency elevation spectrum |

See [processing and I/O](processing_guide.md) for variable names, units, shape
handoffs, direction decisions, and CFOSAT ambiguity handling.

## Plotting

`wasp.plotting.plot_directional_spectrum(..., n_levels=50, partitions=None)` is
the top-level implementation. It accepts degrees, uses period as radius, north
at zero, clockwise direction, a 0–25 s radial range, and defaults `vmax` to the
99th percentile of positive values. `partitions` entries require `Hs`, `Tp`, and
`Dp` keys.

`wasp.plotting_geo.plot_directional_spectrum(..., n_levels=100,
partitions=None)` interpolates with `RectBivariateSpline`, smooths with a Gaussian
filter, and defaults `vmax` to the data maximum. It requires enough valid bins
for bicubic interpolation and is not top-level exported.

Both functions return `(fig, ax)` and do not call `show()`.
Each plotting module also exposes `create_wave_energy_colormap()`, returning its
custom Matplotlib colormap.

## Metrics

- `compute_mad(sim, obs)` returns the time-dimension mean absolute difference.
- `compute_madp(sim, obs)` compares 0–100 percentile curves using NumPy axis 0.
- `metrics(data)` requires Xarray variables `SWH_mod` and `SWH_sat`; it returns
  bias, count, RMSE, means, normalized bias/RMSE, MAD, MADP, and MADC.

These are legacy SWH field metrics, not partition-assignment verification.

## Utilities

- `load_config(config_path=None)` reads an explicit YAML file or searches
  `config.yaml` in the current directory and `examples/`; it raises
  `FileNotFoundError` when absent.
- `format_partition_label(threshold, merge_factor)` formats an integer
  percentile and one-decimal merge factor.
- `build_case_name_cfosat(config)` and `build_case_name(config)` require
  orchestration keys absent from the minimal example YAML.
- `build_output_dir(config)` returns a `Path` and does not create the directory.
- Wave-parameter and conversion helpers are re-exported from `utils.py` for
  compatibility.
