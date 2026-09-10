# Examples, notebooks, and scripts

The files under `examples/` are source-specific research examples, not an
automated test suite. They require external NetCDF inputs under a local
`../data/` layout and write CSVs there. Stored outputs may have been produced by
an older working tree and are not validation evidence for current code.

## Current examples

| File | Workflow | Audit status |
|---|---|---|
| `examples/ww3_spec.ipynb` | WW3 reader → bulk parameters → partition → plot/export | Current WASP imports and keywords; requires `../data/ww3_41001.nc` |
| `examples/ndbc_spec.ipynb` | NDBC reconstruction → partition → plot/export | Current imports and `direction_resolution`; requires NDBC NetCDF data |
| `examples/sentinel_spec.ipynb` | SAR reader/conversion → partition → plot/export | Current imports and reader keywords; requires local SAR NetCDF data |
| `examples/cfosat_spec.ipynb` | CFOSAT conversion/ambiguity helper → transpose → partition → plot/export | Current imports and keywords; ambiguity removal is notebook-local |
| `examples/ex1.ipynb` | Older WW3-oriented example | Duplicates much of `ww3_spec.ipynb`; its hard-coded sample filename is CFOSAT-named, and stored output contains an old `partition_descriptors` `KeyError` |
| `examples/ex2.ipynb` | Older CFOSAT-oriented example | Duplicates much of `cfosat_spec.ipynb`; stored output contains the same old `partition_descriptors` `KeyError` |

`examples/config.yaml` is a notebook configuration example. `load_config()` does
not inject it into readers or partitioning. `min_energy_fraction` is an export
filter implemented by notebook cells, not a `partition_spectrum()` parameter.
The plotting settings are also consumed selectively by notebooks.

## Research notebook

`notebooks/01_analyze_1to2_transitions.ipynb` is a downstream analysis notebook
for previously generated CSV directories. It does not call the WASP package and
contains repeated analysis sections and stored output. Its expected directory
naming/schema may represent earlier validation runs.

## Scripts and tests

The current `scripts/` directory contains no executable workflow. There is no
project `tests/` directory, automated regression suite, or benchmark suite.
Virtual-environment package tests are third-party files and are not WASP tests.

## Recommended notebook backlog

1. Choose `ww3_spec.ipynb` and `cfosat_spec.ipynb` as canonical examples; decide
   whether `ex1.ipynb` and `ex2.ipynb` should be archived.
2. Clear or rerun stale stored outputs only in a controlled, data-available
   validation task; do not treat historical exceptions as current failures.
3. Add a small generated-spectrum notebook that needs no external data.
4. Document or package a reviewed CFOSAT ambiguity policy before moving the
   notebook-local helper into the library.
5. Consolidate repeated sections in the transition notebook without changing
   its scientific calculations.
6. Create unit/regression tests for partitioning, SPR, descriptors, matching,
   readers with fixtures, metrics, and plotting smoke behavior.
