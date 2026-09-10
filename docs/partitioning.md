# Partitioning and wave-property reference

This is the detailed current reference for `src/wasp/partition.py` and its use of
`src/wasp/wave_params.py`. Source code is authoritative. The implementation is
Hanson–Phillips-style scientific software with compatibility behavior that must
not be “cleaned up” without numerical review.

## Core contract

```python
partition_spectrum(
    E,
    frequencies,
    directions_rad,
    energy_threshold=None,
    max_partitions=3,
    threshold_mode="adaptive",
    threshold_percentile=98.0,
    merge_factor=0.315,
    spr_min_peak_prominence=0.4,
    spr_min_relative_energy=0.1,
    spr_diagnostic_filename=None,
)
```

`E` is expected to have shape `(len(frequencies), len(directions_rad))`.
Frequencies are in Hz. Directions are radians on a uniformly spaced complete
circle. Intended density units are `m² s rad⁻¹`. The core treats direction as
periodic, frequency as bounded, and does not change the supplied directional
convention.

The implementation assumes compatible, non-empty arrays and does not perform a
complete input-shape/grid validation. Callers should supply finite non-negative
energy and at least two frequency bins. `calculate_wave_parameters()` replaces
negative and non-finite values with zero, but all partition helpers do not
consistently apply that cleaning.

## 1. Threshold and peaks

Adaptive mode calculates the requested percentile from strictly positive
spectral cells. A spectrum with no positive value returns `None`. Absolute mode
requires `energy_threshold`; omission raises `ValueError`. An unknown mode also
raises `ValueError`.

`identify_spectral_peaks()` examines a 3×3 neighborhood, with circular direction
and hard frequency edges. Cells below `1e-15` are skipped. A cell is a local
maximum when no inspected neighbor has greater energy. Maxima are sorted by
descending cell energy and divided into:

- primary peaks: energy greater than or equal to the threshold;
- secondary peaks: energy strictly below the threshold.

Only primary peaks are limited by `max_partitions`. Above-threshold maxima beyond
that limit are dropped; they are not reclassified as secondary peaks. Peak arrays
contain one-based `(frequency_index, direction_index)` pairs.

## 2. Watershed-style mask

Only primary peaks seed the initial mask. `generate_mask()` repeatedly propagates
labels along the stored steepest-ascent direction codes, using forward/backward
sweeps for at most 50 iterations. Remaining zeros are filled by repeated
eight-neighbor majority voting; a cell with no labelled neighbor falls back to
label `1`.

## 3. Primary-system merge

Peak coordinates are mapped to frequency-direction Cartesian space:

```text
x = f cos(theta)
y = f sin(theta)
```

`calculate_peak_distances()` returns squared Euclidean distances. The spreading
measure `Eip` uses Cartesian first/second moments weighted by
`E * delf * ddir` and normalized by total `m0`. `delf` is the forward
frequency difference with the last value repeated.

For each original primary pair `(i, j)` in index order,
`merge_overlapping_systems()` absorbs label `j + 1` into `i + 1` when both
conditions hold:

```text
distance² <= merge_factor * Eip[i]
distance² <= merge_factor * Eip[j]
```

The standalone helper default is `merge_factor=0.5`; `partition_spectrum()`
passes its own default `0.315`. The merge is a single ordered pair pass. It does
not recompute spreading, compact labels, or decrement `nmask` during this stage.

## 4. Energy and peak properties

Partition energy uses trapezoidal frequency weights and constant
`ddir = 2π / ND`. Array index `0` stores unclassified/discarded energy; positive
indices correspond to labels. Significant wave height is `Hs = 4 sqrt(energy)`.

For each positive label, `Tp` is the reciprocal of the maximum of the
direction-integrated partition spectrum. `Dp` is the direction bin with maximum
energy at that peak frequency, modulo 360 degrees. It is not an energy-weighted
circular mean. Empty labels return `NaN` for `Tp` and `Dp`.

## 5. Secondary Peak Reassessment

SPR runs after the primary merge and before final relabelling. The orchestrator
passes `merge_factor`, `alpha=0.02`, `spr_min_peak_prominence` as
`min_peak_prominence`, and `spr_min_relative_energy` as `min_relative_energy`.
The direct SPR defaults are `alpha=0.05`, `min_peak_prominence=0.4`, and
`min_relative_energy=0.1`. `alpha` is compatibility/diagnostic-only and does not
control promotion. See [the SPR reference](WASP_SPR_specification.md) for the
candidate, saddle, incorporation, promotion, and discard rules.

## 6. Relabelling, moments, and descriptors

After SPR, existing positive labels are remapped in descending `Hs` order;
label `0` is retained. Because primary merging does not decrement the carried
counter, `result["nmask"]` can be larger than the number of positive labels in
the final mask. Empty trailing positions can therefore remain in `Hs`, `energy`,
`Tp`, `Dp`, and moment arrays. To count actual final systems, use positive unique
mask labels or `len(result["partition_descriptors"])`.

`calculate_spectral_moments()` first integrates over direction, then sums
`spec1d[i] * delf[i]` using angular frequency `omega = 2πf` for `m1` and `m2`.
This uses the forward-width `delf` array rather than the trapezoidal weights used
for `energy`/`total_m0`; consumers should not assume those independently
calculated `m0` values are bitwise identical.

Final descriptor dictionaries contain:

```python
{
    "partition": int,
    "tp": float,                    # s
    "dp": float,                    # degrees
    "energy": float,                # m²
    "m0": float,
    "m1": float,
    "m2": float,
    "bandwidth": float,             # dimensionless moment width
    "directional_spreading": float, # degrees
    "spectral_spreading": float,    # same moment width as bandwidth
}
```

Directional spreading uses unweighted frequency-bin sums of `E * ddir` inside
the partition and a circular resultant-length width. These descriptors do not
match the schema/units of `matching.compute_partition_descriptors()`; see
[PCSPM matching](matching.md).

## Result and edge behavior

The result dictionary contains `mask`, `energy`, `Hs`, `Tp`, `Dp`, total bulk
properties, `nmask`, `nmask_before_spr`, primary and secondary peak arrays,
`spr_log`, and `partition_descriptors`. Its `moments` dictionary contains
`total` as an `(m0, m1, m2)` tuple and label-indexed arrays under `m0`, `m1`,
and `m2`. The `spr_diagnostic` value is currently always `None`; when a filename
is supplied, the effect is the PNG side effect described in the SPR reference.

Diagnostic progress and energy checks are printed to standard output. These
messages are not a structured logging API.
