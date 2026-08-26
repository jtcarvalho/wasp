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




def format_partition_label(threshold, merge_factor):
    """
    Examples
    --------
    98, 0.315 -> 98-0.3
    98, 0.500 -> 98-0.5
    95, 0.700 -> 95-0.7
    """

    threshold = int(threshold)

    merge = f"{float(merge_factor):.1f}".rstrip("0").rstrip(".")

    return f"{threshold}-{merge}"



def build_case_name_cfosat(config):
    """
    Build:

    partition_system_analysis_cfosat_obs98-03_ww398-03
    """

    obs_type = config["processing"]["obs_type"]

    obs_cfg = config["partitioning"]["cfosat"]
    ww3_cfg = config["partitioning"]["ww3"]

    obs_label = format_partition_label(
        obs_cfg["threshold_percentile"],
        obs_cfg["merge_factor"],
    )

    ww3_label = format_partition_label(
        ww3_cfg["threshold_percentile"],
        ww3_cfg["merge_factor"],
    )

    return (
        f"partition_system_analysis_"
        f"{obs_type}"
        f"_obs{obs_label}"
        f"_ww3{ww3_label}"
    )

def build_case_name(config):
    """Build an observation-versus-WW3 analysis name from a full config.

    The mapping must contain ``processing.obs_type`` and matching entries under
    ``partitioning``. These orchestration keys are not present in the minimal
    example configuration shipped with the package.
    """

    obs_type = config["processing"]["obs_type"]

    obs_cfg = config["partitioning"][obs_type]
    ww3_cfg = config["partitioning"]["ww3"]

    obs_label = format_partition_label(
        obs_cfg["threshold_percentile"],
        obs_cfg["merge_factor"],
    )

    ww3_label = format_partition_label(
        ww3_cfg["threshold_percentile"],
        ww3_cfg["merge_factor"],
    )

    return (
        f"partition_system_analysis_"
        f"{obs_type}"
        f"_obs{obs_label}"
        f"_ww3{ww3_label}"
    )



def build_output_dir(config):
    """Build the configured analysis output path without creating it.

    In addition to the keys required by :func:`build_case_name`, the mapping
    must contain ``paths.output_root``.
    """

    output_dir = (
        Path(config["paths"]["output_root"])
        / build_case_name(config)
    )

    station_filter = config["processing"].get("station_filter")

    if station_filter:
        if isinstance(station_filter, (list, tuple)):
            output_dir /= "_".join(map(str, station_filter))
        else:
            output_dir /= str(station_filter)

    return output_dir
