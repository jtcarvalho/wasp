# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 2024-12-20

### Added

- Full support for NDBC buoy spectral data partitioning
- Adaptive threshold partitioning algorithm with source-specific parameters
- Source-specific partitioning configurations:
  - SAR (Sentinel-1): High threshold (99.5%), conservative merging (0.3)
  - WW3 (WaveWatch3): Moderate threshold (99.0%), balanced merging (0.5)
  - NDBC (Buoys): Lower threshold (98.0%), aggressive merging (0.7)
- NDBC I/O functions (`io_ndbc.py`) for reading buoy spectral data
- Convergence-based iteration in mask generation (up to 50 iterations)
- Comprehensive documentation for partitioning algorithm and implementation

### Changed

- `partition_spectrum()` function now accepts configurable parameters:
  - `threshold_mode`: 'adaptive' (default) or 'fixed'
  - `threshold_percentile`: Percentile for adaptive threshold (default: 99.0)
  - `merge_factor`: Factor for merging overlapping systems (default: 0.5)
  - `max_partitions`: Maximum number of partitions (default: 5)
- `merge_overlapping_systems()` now accepts configurable `merge_factor` parameter
- `generate_mask()` iterates until convergence instead of fixed iterations
- Updated all example scripts (01_partition_sar.py, 02_partition_ww3.py, 03_partition_ndbc.py) with adaptive parameters
- Updated all notebooks (sarspec.ipynb, ww3spec.ipynb, ndbcspec.ipynb) with source-specific configurations

### Fixed

- Over-partitioning issue in NDBC data due to low directional resolution
- Fixed parameter approach that caused inconsistent partitioning across data sources
- Energy conservation in partitioning algorithm

### Documentation

- Added `PARTITION_ALGORITHM_ANALYSIS.md` - Detailed analysis of the watershed algorithm
- Added `PARTITION_IMPROVEMENTS.md` - Documentation of adaptive partitioning improvements
- Added `GUIA_IMPLEMENTACAO_PARTICIONAMENTO.md` - Implementation guide (Portuguese)

## [0.1.1] - 2024-12-17

### Added

- Initial release with SAR and WW3 support
- Basic spectral partitioning using watershed algorithm (Hanson & Phillips, 2001)
- Wave parameter calculations
- Plotting utilities for directional spectra

### Features

- SAR (Sentinel-1) spectral data processing
- WW3 (WaveWatch3) model data processing
- 2D directional spectrum partitioning
- Wave parameter computation (Hs, Tp, Dp, moments)
- Visualization tools for spectra and partitions

## [0.1.0] - Initial Development

### Added

- Core package structure
- Basic partitioning algorithm implementation
- I/O modules for SAR and WW3 data
