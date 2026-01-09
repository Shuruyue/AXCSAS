# Fitting Module
"""
Peak fitting module for XRD analysis.
Includes peak detection, Pseudo-Voigt fitting, quality metrics, and hkl assignment.
"""

from .peak_detection import PeakDetector, PeakInfo, find_peaks
from .pseudo_voigt import PseudoVoigt, PseudoVoigtParams, pseudo_voigt_function
from .lm_optimizer import LMOptimizer, FitResult, fit_peaks
from .quality_metrics import (
    QualityLevel,
    FitQualityReport,
    calculate_r_wp,
    calculate_rss,
    calculate_r_squared,
    validate_fit_parameters,
    generate_quality_report,
)
from .hkl_assignment import (
    PeakAssignment,
    assign_hkl,
    assign_hkl_detailed,
    assign_all_peaks,
    format_hkl,
    JCPDS_COPPER_PEAKS,
)

__all__ = [
    # Peak detection
    "PeakDetector",
    "PeakInfo",
    "find_peaks",
    # Pseudo-Voigt
    "PseudoVoigt",
    "PseudoVoigtParams",
    "pseudo_voigt_function",
    # LM Optimizer
    "LMOptimizer",
    "FitResult",
    "fit_peaks",
    # Quality metrics
    "QualityLevel",
    "FitQualityReport",
    "calculate_r_wp",
    "calculate_rss",
    "calculate_r_squared",
    "validate_fit_parameters",
    "generate_quality_report",
    # HKL assignment
    "PeakAssignment",
    "assign_hkl",
    "assign_hkl_detailed",
    "assign_all_peaks",
    "format_hkl",
    "JCPDS_COPPER_PEAKS",
]
