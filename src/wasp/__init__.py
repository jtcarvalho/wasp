"""Public package facade for WASP wave-spectrum partitioning.

The top-level API exposes the stable partitioning, bulk-parameter, plotting,
and configuration entry points. Source-specific loaders and PCSPM matching
remain available from their respective :mod:`wasp.io_*` and
:mod:`wasp.matching` modules.
"""

__version__ = "2.0.0"
__author__ = "J.T. Carvalho"
__email__ = "jtcarvalho@gmail.com"

from .partition import partition_spectrum
from .wave_params import calculate_wave_parameters
from .plotting import plot_directional_spectrum
from .utils import load_config

__all__ = [
    "partition_spectrum",
    "calculate_wave_parameters",
    "plot_directional_spectrum",
    "load_config",
    "__version__",
]
