# Secondary Peak Reassessment (SPR): current implementation

This document describes the behavior implemented by
`wasp.partition.secondary_peak_reassessment()` in the released v2.0.0 line. It
does not prescribe a future calibration.

## Position in the pipeline

```text
all local maxima
  → threshold split into primary and secondary peaks
  → primary-seed watershed
  → Hanson–Phillips merge of primary systems
  → preliminary energy/Tp/Dp
  → Secondary Peak Reassessment
  → energy renumbering and final parameters
```

Primary peaks have energy greater than or equal to the selected threshold.
Secondary peaks are local maxima strictly below it. Only primary peaks seed the
initial partition mask. `max_partitions` limits only primary peaks;
above-threshold maxima dropped by that limit do not become SPR candidates.

## Candidate regions

When `ICOD` and the primary peaks are supplied, SPR constructs a separate mask
seeded by all primary and secondary peaks. Each secondary basin in this
all-peak mask defines that candidate's region and integrated `Esystem`.

This candidate-region mask is used only inside SPR. The input mask remains the
primary-only watershed after the Hanson–Phillips merge.

## Compatibility with existing systems

For every secondary candidate, SPR compares its peak with every positive label
currently present in the working mask. Peak locations are mapped into Cartesian
frequency-direction space:

```text
x = f cos(theta)
y = f sin(theta)
```

The candidate is compatible with a partition only when both Hanson–Phillips
style conditions hold:

```text
distance² <= merge_factor × Eip(candidate)
distance² <= merge_factor × Eip(partition)
```

If several partitions are compatible, the partition with the smallest squared
distance receives the candidate region. Primary partitions are not merged with
one another during SPR.

## Independence and promotion

If no compatible partition exists, SPR finds interfaces between the candidate
region and labelled neighboring regions. For each neighbor it retains the
largest interface saddle estimate, calculated from the minimum energy on the
two sides of adjacent boundary cells.

Peak prominence is:

```text
prominence = (Epeak - Esaddle) / Epeak
```

The candidate is considered independent when:

```text
prominence >= min_peak_prominence
```

The active default is `min_peak_prominence=0.4`.

The candidate energy fraction is:

```text
relative_energy = Esystem / Etotal
```

Promotion currently requires both independence and:

```text
relative_energy >= min_relative_energy
```

The active SPR default is `min_relative_energy=0.1`.

## Decisions

Each secondary region receives exactly one logged decision:

| Decision | Current condition | Mask action |
|---|---|---|
| `INCORPORATED` | Hanson–Phillips compatible | Assign region to nearest compatible partition |
| `NEW_SYSTEM` | Incompatible, independent, and energetic | Assign a new positive label |
| `INCORPORATED` | Incompatible but not independent | Assign region to the nearest existing partition by normalized distance |
| `DISCARDED` | Independent but below the energy requirement | Assign region to label 0 |

After all candidates are processed, partition energy, Tp, and Dp are
recalculated. Final energy-based renumbering is performed by
`partition_spectrum()` after SPR returns.

## Parameters and compatibility notes

- `merge_factor` controls the compatibility comparison.
- `min_peak_prominence` controls saddle separation.
- `min_relative_energy` controls promotion in the current implementation.
- `alpha` is retained in the SPR signature and diagnostics for compatibility,
  but it is not used by the active promotion condition.
- `partition_spectrum()` passes `alpha=0.02` and forwards
  `spr_min_peak_prominence`. It also forwards `spr_min_relative_energy` to SPR
  as `min_relative_energy`; both parameters default to 0.1.
- Direct calls with secondary candidates must supply a numeric `merge_factor`;
  the signature default `None` has no implemented fallback. The orchestrator
  always supplies its numeric value.
- `alpha` must be non-negative and `min_peak_prominence` must satisfy
  `0 <= value < 1`. `min_relative_energy` is not range-validated.
- `spr_log["secondary_peak"]` uses zero-based indices, whereas the
  `secondary_peaks` array returned by peak identification is one-based.

These compatibility details describe the current v2.0.0 behavior.

## Diagnostic plot

If `spr_diagnostic_filename` is supplied to `partition_spectrum()`,
`plot_spr_diagnostic()` writes a PNG showing:

- the original energy spectrum;
- final SPR labels;
- primary peaks;
- incorporated, promoted, and discarded secondary peaks.

No diagnostic file is produced when the filename is omitted.
The result dictionary's `spr_diagnostic` field remains `None` even when the PNG
side effect is requested.

## Labels and partition count

SPR initializes its label counter from the incoming primary `nmask` and adds one
for every promoted candidate. The preceding primary merge does not decrement
that counter when it absorbs a label. Final energy-based relabelling compacts
labels that actually exist, but `result["nmask"]` can still exceed the number of
positive labels in the returned mask. Count non-empty final systems from the
mask or `partition_descriptors` when this distinction matters.
