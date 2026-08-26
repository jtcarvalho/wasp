# PCSPM matching: current implementation

This is the single documentation source for the matching behavior implemented
in `src/wasp/matching.py` for the v2.0.0 pre-release. Older methodology files in
this directory are retained as historical research material and do not define
the package API.

## Public interfaces

### `match_partition_properties()`

Matches two sequences of precomputed partition-property dictionaries. This is
the high-level API for exported CSV properties.

Each descriptor must contain finite numeric values for:

```python
{
    "partition": 1,                 # optional identifier
    "tp": 14.0,                     # seconds
    "dp": 220.0,                    # degrees
    "energy": 0.8,                  # positive integrated energy
    "bandwidth": 0.2,               # non-negative; caller-defined contract
    "directional_spreading": 0.31,  # radians; non-negative
    "spectral_spreading": 0.2,      # non-negative; currently validation-only
}
```

The active cost interprets `directional_spreading` as radians. The implementation
does not translate `peak_frequency`/`peak_direction` keys to `tp`/`dp`.

### `match_spectral_partitions()`

Compatibility wrapper which delegates two already-computed descriptor sequences
to `match_partition_properties()`. Despite the historical name, it does not
accept spectra and does not call `compute_partition_descriptors()`.

Its default `method="weighted"` selects legacy result/dummy behavior; it does
not activate the earlier four-weight equation.

### `compute_partition_descriptors()`

Standalone producer accepting `E2d`, frequencies, directions in radians, and a
label mask. It returns:

```python
{
    "partition": ...,
    "peak_frequency": ...,
    "peak_direction": ...,
    "energy": ...,
    "bandwidth": ...,
    "directional_spreading": ...,
    "spectral_spreading": ...,
}
```

Its bandwidth is an energy-weighted frequency standard deviation in Hz,
directional spreading is in radians, and spectral spreading is Cartesian
frequency-direction variance normalized by partition energy.

This schema is not directly consumable by `match_partition_properties()`.
Likewise, `partition.build_partition_descriptors()` supplies `tp`/`dp` but
reports directional spreading in degrees, so passing its output directly would
give the matching cost the wrong angular unit. The pre-release cleanup documents
these known contract gaps instead of changing them, because choosing canonical
peak and spreading definitions can change matching results.

## Active pair cost

For each observed/modelled descriptor pair, the matcher converts Tp to peak
frequency only for an internal Cartesian-distance calculation:

```text
f = 1 / Tp
x = f cos(Dp)
y = f sin(Dp)
```

It calculates a normalized distance using `spectral_spreading`, but that value
is not part of the current active cost.

The returned cost is:

```text
C = Ctp + Cdirection + Cenergy + Cbandwidth + Cspread
```

with:

```text
Ctp = 0.5 [ln(Tp_obs / Tp_mod) / 0.1]²

sigma_theta = clip(
    sqrt(spread_obs² + spread_mod²),
    5 degrees,
    35 degrees,
)

Cdirection = 0.5 delta_theta² / sigma_theta²
             × [1 + (|delta_theta| / 40 degrees)⁴]

Cenergy    = 0.5 ln(E_obs / E_mod)²
Cbandwidth = 0.5 ln(max(BW_obs, 1e-3) / max(BW_mod, 1e-3))²
Cspread    = 0.5 ln(max(S_obs, 1e-3) / max(S_mod, 1e-3))²
```

`delta_theta` is the signed shortest circular difference in radians. Energy
must be positive; bandwidth and spreading values must be non-negative.

The parameters `alpha`, `beta`, `gamma`, and `delta` remain in the public
signature for compatibility but do not weight the active terms. The historical
weighted distance equation is commented out in the source. Likewise,
`max_direction_cost` is currently unused.

## Candidate matrix and optional filter

The complete observed-by-modelled candidate matrix is calculated unless
`physical_filter=True`. With the filter enabled, a candidate receives sentinel
cost `1e12` when either:

- circular direction difference exceeds `max_direction_deg`; or
- absolute Tp difference exceeds the adaptive tolerance returned by
  `get_tp_tolerance()`.

The adaptive default tolerance is:

| Observed Tp | Tolerance |
|---|---|
| `< 10 s` | `0.20 × Tp` |
| `10–14 s` | `0.25 × Tp` |
| `>= 14 s` | `0.30 × Tp` |

The default is `physical_filter=False`.

## Hungarian assignment and dummy systems

The candidate matrix is augmented with dummy rows and columns before calling
`scipy.optimize.linear_sum_assignment()`.

For default `method="probabilistic"`, dummy cost is:

```text
max(-ln(probability_threshold), 95th percentile of finite candidate costs)
```

The default probability threshold is 0.05. In other methods, dummy cost is the
75th percentile of finite candidate costs.

A real pair is accepted only when its candidate cost is strictly less than the
dummy cost. Consequently, some systems can remain unmatched even when systems
exist on both sides.

## Return structure

For non-empty candidate matrices, the returned dictionary contains:

- `matched_pairs`: accepted pairs and diagnostics;
- `unmatched_observed`: original observed dictionaries not accepted;
- `unmatched_modeled`: original modelled dictionaries not accepted;
- `cost_matrix`: unaugmented candidate matrix;
- `candidate_costs`: diagnostics for every candidate;
- `matching_ranking`: best/second-best diagnostics per observed system.

In probabilistic mode, each pair also contains individual cost components and
`match_probability = exp(-cost)`.

When either input side is empty, the function returns the first four core
collections but omits `candidate_costs` and `matching_ranking`.

## Example

```python
from wasp.matching import match_partition_properties

result = match_partition_properties(observed, modeled)

for pair in result["matched_pairs"]:
    print(pair["observed"]["partition"], pair["modeled"]["partition"])

print(result["unmatched_observed"])
print(result["unmatched_modeled"])
```

## Colocation helpers

`matching.py` also retains public helpers unrelated to partition assignment:

- `haversine_distance()`;
- `check_ndbc_has_spectral_data()`;
- `scan_ndbc_stations()`;
- `find_sar_ndbc_matches()`;
- `add_ww3_info()`.

They support filesystem discovery and spatial/temporal association. Their
placement is historical and retained for compatibility.
