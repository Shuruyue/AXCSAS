# Preprocessing Module
"""
Data preprocessing module for XRD analysis.
Includes data loading, validation, smoothing, background subtraction,
Kα2 stripping, and pipeline orchestration.
"""

from .data_loader import XRDDataLoader, load_xrd_data
from .smoothing import SavitzkyGolayFilter, smooth_xrd_data
from .background import BackgroundSubtractor, subtract_background
from .kalpha_strip import KalphaStripper, strip_kalpha2
from .validation import (
    XRDDataset,
    DataValidationResult,
    ValidationWarning,
    WarningLevel,
    validate_xrd_data,
    check_negative_values,
)
from .pipeline import (
    PreprocessingPipeline,
    PreprocessingResult,
    PreprocessingStep,
    should_apply_kalpha_stripping,
)

# Aliases for convenience
apply_smoothing = smooth_xrd_data

__all__ = [
    # Data loading
    "XRDDataLoader",
    "load_xrd_data",
    # Validation
    "XRDDataset",
    "DataValidationResult",
    "ValidationWarning",
    "WarningLevel",
    "validate_xrd_data",
    "check_negative_values",
    # Smoothing
    "SavitzkyGolayFilter",
    "smooth_xrd_data",
    "apply_smoothing",
    # Background
    "BackgroundSubtractor",
    "subtract_background",
    # Kα2 stripping
    "KalphaStripper",
    "strip_kalpha2",
    # Pipeline
    "PreprocessingPipeline",
    "PreprocessingResult",
    "PreprocessingStep",
    "should_apply_kalpha_stripping",
]
