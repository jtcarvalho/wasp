# **WASP** - **WA**ve **S**pectra **P**artitioning

Watershed Algorithm for partitioning ocean wave spectra from model  and observation

[![PyPI version](https://img.shields.io/pypi/v/wasp-ocean.svg)](https://pypi.org/project/wasp-ocean/)
[![Python Version](https://img.shields.io/pypi/pyversions/wasp-ocean.svg)](https://pypi.org/project/wasp-ocean/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19415744.svg)](https://doi.org/10.5281/zenodo.19415744)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## 📋 What is WASP?

WASP focuses exclusively on **spectral partitioning** - the process of separating ocean wave spectra into individual wave systems (partitions) via watershed algorithm. Each partition represents a distinct wave system characterized by significant wave height (Hs), peak period (Tp), and direction (Dp).

**WASP handles:**

- ✅ Spectral partitioning using watershed algorithm
- ✅ Processing SAR (Sentinel-1 and CFOSAT), NDBC and WW3 model spectra
- ✅ Extracting wave parameters (Hs, Tp, Dp) for each partition

## 🚀 Installation

> ⚠️ **IMPORTANT:** Python 3.10 or higher is required.

### Install from PyPI (Recommended)

```bash
pip install wasp-ocean
```

### Verify Installation

```bash
# Test the import
python -c "import wasp; print(f'WASP version: {wasp.__version__}')"

# Test main functions
python -c "from wasp import partition_spectrum, calculate_wave_parameters; print('✓ Installation successful!')"
```

### Development Installation

For development or local modifications:

```bash
# Clone the repository
git clone https://github.with/jtcarvalho/wasp.git
cd wasp

# Install in editable mode
pip install -e .
```

## 📦 Key Dependencies

- **Python >= 3.10** (required)
- **NumPy >= 2.1.0** (required for `np.trapezoid`)
- pandas >= 2.2.0
- xarray >= 2024.11.0
- matplotlib >= 3.8.0
- scipy >= 1.14.0
- scikit-image >= 0.22.0
- netCDF4 >= 1.5.4

> ⚠️ **Note:** NumPy < 2.1.0 will cause errors as `np.trapezoid` is not available.


## 🏗️ Project Structure

```
wasp/
├── src/
│   └── wasp/              # Main package
│       ├── partition.py   # Watershed partitioning algorithm
│       ├── wave_params.py # Wave parameter calculations
│       ├── io_sar.py      # SAR Sentinel data I/O
│       ├── io_cfosat.py   # SAR CFOSAT data I/O
│       ├── io_ww3.py      # WW3 data I/O
│       ├── io_ndbc.py     # NDBC data I/O
│       └── utils.py       # Utility functions
└── examples/              # Usage examples

```

## 📄 License

This project is licensed under the MIT License.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📧 Contact

For questions or issues, please open an issue on [GitHub](https://github.with/jtcarvalho/wasp/issues).
