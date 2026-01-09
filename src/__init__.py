# AXCSAS - Advanced XRD Crystallite Size Analysis System
"""
AXCSAS - Advanced XRD Crystallite Size Analysis System

A modular Python package for automated XRD pattern analysis,
crystallite size calculation, and texture analysis.
"""

from . import preprocessing
from . import fitting
from . import physics
from . import validation
from . import utils

__version__ = "0.1.0"
__author__ = "AXCSAS Team"

__all__ = [
    "preprocessing",
    "fitting",
    "physics",
    "validation",
    "utils",
]
