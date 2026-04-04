"""
Utility functions for WASP — configuration loading and shared helpers.

Functions that were previously defined here (calculate_wave_parameters,
spectrum1d_from_2d, convert_meteorological_to_oceanographic,
convert_spectrum_units) are now canonical in wave_params.py and are
re-exported here for backward compatibility.
"""

import yaml
from pathlib import Path

# Re-export from canonical modules so existing code using
# `from wasp.utils import …` continues to work without modification.
from .wave_params import (
    calculate_wave_parameters,
    spectrum1d_from_2d,
    convert_meteorological_to_oceanographic,
    convert_spectrum_units,
)


def load_config(config_path=None):
    """
    Load configuration from a YAML file.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to configuration file. If provided, loads this file directly.
        If not provided (None), searches for config.yaml in default locations.

    Returns
    -------
    dict
        Configuration dictionary with paths and parameters.

    Raises
    ------
    FileNotFoundError
        If the config file cannot be found.
    """
    if config_path is not None:
        config_file = Path(config_path)
        if not config_file.exists():
            raise FileNotFoundError(
                f"Configuration file '{config_path}' not found."
            )
    else:
        possible_paths = [
            Path.cwd() / 'config.yaml',
            Path.cwd() / 'examples' / 'config.yaml',
            Path(__file__).parent.parent.parent / 'examples' / 'config.yaml',
        ]

        config_file = None
        for path in possible_paths:
            if path.exists():
                config_file = path
                break

        if config_file is None:
            raise FileNotFoundError(
                "Configuration file 'config.yaml' not found.\n"
                "Please provide the path explicitly or place config.yaml in the "
                "working directory or examples/ folder."
            )

    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


