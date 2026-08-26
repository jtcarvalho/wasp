# WASP — Wave Spectra Partitioning

[![PyPI version](https://img.shields.io/pypi/v/wasp-ocean.svg)](https://pypi.org/project/wasp-ocean/)
[![Python versions](https://img.shields.io/pypi/pyversions/wasp-ocean.svg)](https://pypi.org/project/wasp-ocean/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19415744.svg)](https://doi.org/10.5281/zenodo.19415744)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

WASP separates a two-dimensional directional ocean-wave spectrum into
physically distinct wave systems. Its core is a watershed-style implementation
based on Hanson and Phillips (2001), followed by primary-system merging and
Secondary Peak Reassessment (SPR).

The package also provides source-specific spectrum readers, bulk wave-parameter
calculations, polar plotting, and Physics-Constrained Spectral Partition
Matching (PCSPM) through a Hungarian assignment interface.

## Capabilities

- Partition spectra shaped `(frequency, direction)` into labelled systems.
- Calculate energy, significant wave height (`Hs`), peak period (`Tp`), peak
  direction (`Dp`), and spectral moments.
- Reassess sub-threshold peaks with SPR after the primary Hanson–Phillips merge.
- Match observed and modelled partition-property dictionaries with PCSPM.
- Load or reconstruct spectra from Sentinel/CMEMS SAR, CFOSAT SWIM, NDBC, and
  WaveWatch III data.
- Plot directional spectra in polar period–direction coordinates.

## Installation

WASP requires Python 3.10 or newer.

```bash
pip install wasp-ocean
```

For development from a local clone:

```bash
git clone https://github.com/jtcarvalho/wasp.git
cd wasp
python -m pip install -e .
```

## Minimal partitioning example

```python
import numpy as np
from wasp import partition_spectrum

frequencies = np.linspace(0.05, 0.30, 32)
directions_rad = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)

# E2d must have shape (len(frequencies), len(directions_rad)).
E2d = np.zeros((frequencies.size, directions_rad.size))
E2d[8, 4] = 1.0
E2d[18, 15] = 0.6

result = partition_spectrum(
    E2d,
    frequencies,
    directions_rad,
    threshold_mode="adaptive",
    threshold_percentile=98.0,
    max_partitions=3,
    merge_factor=0.315,
)

if result is not None:
    print(result["nmask"])
    print(result["Hs"][1 : result["nmask"] + 1])
    print(result["Tp"][1 : result["nmask"] + 1])
    print(result["Dp"][1 : result["nmask"] + 1])
```

The current result dictionary contains the final mask and parameters, total
parameters, primary and secondary peaks, SPR decisions, moment arrays, and
partition descriptors. Directions are not converted by the partitioning core;
their physical convention is inherited from `directions_rad`.

## Partitioning workflow

```text
threshold selection
  → local peak identification
  → primary-seed watershed mask
  → Hanson–Phillips primary merge
  → Secondary Peak Reassessment
  → energy-based renumbering
  → final parameters, moments, and descriptors
```

SPR compares each sub-threshold candidate with existing systems using spectral
distance and spreading. A compatible candidate is incorporated. Otherwise,
the current implementation uses candidate energy and peak-to-saddle prominence
to promote, incorporate, or discard its watershed region. See
[`WASP_SPR_specification.md`](WASP_SPR_specification.md) for the implemented
decision flow.

## PCSPM matching

The descriptor-based matching API is available from `wasp.matching`:

```python
from wasp.matching import match_partition_properties

observed = [{
    "partition": 1,
    "tp": 14.0,
    "dp": 220.0,
    "energy": 0.80,
    "bandwidth": 0.20,
    "directional_spreading": 0.31,  # radians in the active matching cost
    "spectral_spreading": 0.20,
}]
modeled = [{
    "partition": 1,
    "tp": 13.5,
    "dp": 215.0,
    "energy": 0.75,
    "bandwidth": 0.18,
    "directional_spreading": 0.35,  # radians in the active matching cost
    "spectral_spreading": 0.18,
}]

matching = match_partition_properties(observed, modeled)
print(matching["matched_pairs"])
print(matching["unmatched_observed"], matching["unmatched_modeled"])
```

The active matching cost interprets `directional_spreading` as radians.
`partition.build_partition_descriptors()` currently reports that field in
degrees, so its output requires an explicit unit/schema adapter before use with
the matcher; the package does not silently convert it.

The active implementation builds a candidate-cost matrix from peak-period,
circular-direction, energy, bandwidth, and directional-spreading penalties,
then solves an augmented Hungarian assignment. Exact equations, descriptor
requirements, dummy assignments, and compatibility parameters are documented
in [`docs/matching.md`](docs/matching.md).

## Supported data sources

| Source | Module | Current adapter behavior |
|---|---|---|
| Sentinel/CMEMS SAR | `wasp.io_sar` | Reads current `wave_spec` names or legacy `oswPolSpec`; converts wavenumber spectra to frequency-direction spectra; preserves the supplied going-to direction coordinate. |
| CFOSAT SWIM | `wasp.io_cfosat` | Reads L2/L2PBOX products, converts slope/wavenumber spectra, expands directional ambiguity, optionally normalizes to processor Hs, and returns a dictionary with spectrum shape `(direction, frequency)`. |
| NDBC | `wasp.io_ndbc` | Reconstructs a 2-D spectrum from directional Fourier coefficients and returns a coming-from directional grid. |
| WaveWatch III | `wasp.io_ww3` | Reads `efth` and rotates stored propagation directions by 180° to the adapter's coming-from convention. |

The adapters do not yet return one common object. Check shape, units, and
direction convention before comparing different sources. The package core
expects spectral density in `m² s rad⁻¹` on a uniformly spaced directional grid.

## Plotting and analysis

```python
from wasp import plot_directional_spectrum

fig, ax = plot_directional_spectrum(
    E2d,
    frequencies,
    np.degrees(directions_rad),
)
```

`wasp.plotting_geo` provides an optional interpolated presentation variant.
`wasp.metrics` contains legacy Xarray significant-wave-height metrics. Research
notebooks and generated analysis artefacts are included for inspection, but the
current checkout does not contain executable sensor-validation scripts.

## Documentation

- [Architecture](README_archi.md)
- [Processing and I/O guide](docs/processing_guide.md)
- [Current PCSPM implementation](docs/matching.md)
- [SPR implementation](WASP_SPR_specification.md)

## Citation

Please cite WASP using [`CITATION.cff`](CITATION.cff). The archived software DOI
is [10.5281/zenodo.19415744](https://doi.org/10.5281/zenodo.19415744).

The partitioning methodology is based on:

> Hanson, J. L., & Phillips, O. M. (2001). Automated Analysis of Ocean Surface
> Directional Wave Spectra. *Journal of Atmospheric and Oceanic Technology*,
> 18(2), 277–293.

## Project links

- Repository: https://github.com/jtcarvalho/wasp
- PyPI: https://pypi.org/project/wasp-ocean/
- Issues: https://github.com/jtcarvalho/wasp/issues

WASP is distributed under the MIT License.
