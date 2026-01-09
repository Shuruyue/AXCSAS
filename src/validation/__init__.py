"""
Validation Module
Contains validation functions for XRD analysis results.
"""

from .goodness_of_fit import (
    calculate_rwp,
    calculate_r_squared,
    calculate_chi_squared,
    FitQualityResult
)
from .error_analysis import (
    ErrorAnalyzer,
    validate_size_range,
    check_broadening_ratio
)

__all__ = [
    "calculate_rwp",
    "calculate_r_squared", 
    "calculate_chi_squared",
    "FitQualityResult",
    "ErrorAnalyzer",
    "validate_size_range",
    "check_broadening_ratio",
]
