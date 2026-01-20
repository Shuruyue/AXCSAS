"""
XRD Data Validation Module XRD 資料驗證模組
===========================================

Input validation and data quality checks for XRD preprocessing.
XRD 預處理的輸入驗證和資料品質檢查。
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional
from enum import Enum


class WarningLevel(Enum):
    """Severity level for validation warnings."""
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass
class ValidationWarning:
    """A single validation warning."""
    level: WarningLevel
    code: str
    message: str
    

@dataclass
class DataValidationResult:
    """
    Result of XRD data validation.
    
    Attributes:
        is_valid: True if data passes all critical checks
        warnings: List of validation warnings
        stats: Dictionary of computed statistics
    """
    is_valid: bool
    warnings: List[ValidationWarning] = field(default_factory=list)
    stats: dict = field(default_factory=dict)
    
    def has_errors(self) -> bool:
        """Check if any error-level warnings exist."""
        return any(w.level == WarningLevel.ERROR for w in self.warnings)
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            f"Validation: {'PASS' if self.is_valid else 'FAIL'}",
            f"Warnings: {len(self.warnings)}"
        ]
        for w in self.warnings:
            lines.append(f"  [{w.level.value.upper()}] {w.code}: {w.message}")
        return "\n".join(lines)


@dataclass
class XRDDataset:
    """
    Standardized XRD data container.
    
    Provides a consistent interface for XRD data with computed properties
    and validation support.
    
    Attributes:
        two_theta: 2θ angle array in degrees
        intensity: Intensity array (counts or arbitrary units)
        metadata: Optional metadata dictionary
    """
    two_theta: np.ndarray
    intensity: np.ndarray
    metadata: dict = field(default_factory=dict)
    
    def __post_init__(self):
        """Ensure arrays are numpy arrays."""
        self.two_theta = np.asarray(self.two_theta, dtype=float)
        self.intensity = np.asarray(self.intensity, dtype=float)
    
    @property
    def step_size(self) -> float:
        """Calculate average 2θ step size in degrees."""
        if len(self.two_theta) < 2:
            return 0.0
        return float(np.mean(np.diff(self.two_theta)))
    
    @property
    def theta_range(self) -> Tuple[float, float]:
        """Return (min, max) 2θ range."""
        return float(self.two_theta.min()), float(self.two_theta.max())
    
    @property
    def n_points(self) -> int:
        """Number of data points."""
        return len(self.two_theta)
    
    def validate(self) -> DataValidationResult:
        """Run validation on this dataset."""
        return validate_xrd_data(self.two_theta, self.intensity)


# =============================================================================
# Validation Functions
# =============================================================================

# Validation thresholds
TWO_THETA_MIN = 10.0   # degrees
TWO_THETA_MAX = 150.0  # degrees
MIN_DATA_POINTS = 100
STEP_UNIFORMITY_TOLERANCE = 0.01  # degrees


def validate_xrd_data(
    two_theta: np.ndarray,
    intensity: np.ndarray
) -> DataValidationResult:
    """
    Validate XRD data for physical reasonableness.
    
    Checks performed:
    1. 2θ range within [10°, 150°]
    2. Intensity values non-negative
    3. Sufficient data points (>100)
    4. Uniform step size
    5. Array length consistency
    
    Args:
        two_theta: 2θ angle array in degrees
        intensity: Intensity array
        
    Returns:
        DataValidationResult with validation status and warnings
        
    Physical Rationale:
        - 2θ < 10°: Typically dominated by direct beam/air scatter
        - 2θ > 150°: Rarely used, geometric limitations
        - Negative intensity: Indicates background subtraction error
        - Non-uniform steps: May cause issues in FFT-based algorithms
    """
    warnings = []
    stats = {}
    is_valid = True
    
    # Check array lengths match
    if len(two_theta) != len(intensity):
        warnings.append(ValidationWarning(
            level=WarningLevel.ERROR,
            code="ARRAY_MISMATCH",
            message=f"Array lengths differ: 2θ={len(two_theta)}, I={len(intensity)}"
        ))
        return DataValidationResult(is_valid=False, warnings=warnings, stats=stats)
    
    # Check data points count
    stats['n_points'] = len(two_theta)
    if len(two_theta) < MIN_DATA_POINTS:
        warnings.append(ValidationWarning(
            level=WarningLevel.WARNING,
            code="INSUFFICIENT_POINTS",
            message=f"Only {len(two_theta)} points (recommend ≥{MIN_DATA_POINTS})"
        ))
    
    # Check 2θ range
    theta_min, theta_max = float(two_theta.min()), float(two_theta.max())
    stats['theta_range'] = (theta_min, theta_max)
    
    if theta_min < TWO_THETA_MIN:
        warnings.append(ValidationWarning(
            level=WarningLevel.WARNING,
            code="LOW_ANGLE",
            message=f"Data starts at {theta_min:.2f}° (below recommended {TWO_THETA_MIN}°)"
        ))
    
    if theta_max > TWO_THETA_MAX:
        warnings.append(ValidationWarning(
            level=WarningLevel.WARNING,
            code="HIGH_ANGLE",
            message=f"Data ends at {theta_max:.2f}° (above recommended {TWO_THETA_MAX}°)"
        ))
    
    # Check for negative intensity
    n_negative = np.sum(intensity < 0)
    stats['n_negative'] = int(n_negative)
    
    if n_negative > 0:
        pct_negative = 100 * n_negative / len(intensity)
        warnings.append(ValidationWarning(
            level=WarningLevel.WARNING if pct_negative < 5 else WarningLevel.ERROR,
            code="NEGATIVE_INTENSITY",
            message=f"{n_negative} negative values ({pct_negative:.1f}%) - likely background subtraction issue"
        ))
        if pct_negative >= 5:
            is_valid = False
    
    # Check step uniformity
    if len(two_theta) >= 2:
        steps = np.diff(two_theta)
        step_mean = np.mean(steps)
        step_std = np.std(steps)
        stats['step_size'] = float(step_mean)
        stats['step_std'] = float(step_std)
        
        if step_std > STEP_UNIFORMITY_TOLERANCE:
            warnings.append(ValidationWarning(
                level=WarningLevel.WARNING,
                code="NON_UNIFORM_STEP",
                message=f"Step size varies: mean={step_mean:.4f}°, std={step_std:.4f}°"
            ))
    
    # Check for monotonically increasing 2θ
    if not np.all(np.diff(two_theta) > 0):
        warnings.append(ValidationWarning(
            level=WarningLevel.ERROR,
            code="NON_MONOTONIC",
            message="2θ values are not monotonically increasing"
        ))
        is_valid = False
    
    return DataValidationResult(
        is_valid=is_valid,
        warnings=warnings,
        stats=stats
    )


def check_negative_values(
    intensity: np.ndarray,
    auto_correct: bool = False
) -> Tuple[np.ndarray, List[int]]:
    """
    Check for and optionally correct negative intensity values.
    
    Physical Context:
        Negative values typically result from over-aggressive background
        subtraction and should be flagged for review.
    
    Args:
        intensity: Intensity array
        auto_correct: If True, set negative values to zero
        
    Returns:
        Tuple of (corrected_intensity, indices_of_negative_values)
    """
    negative_indices = np.where(intensity < 0)[0].tolist()
    
    if auto_correct and negative_indices:
        corrected = intensity.copy()
        corrected[corrected < 0] = 0
        return corrected, negative_indices
    
    return intensity, negative_indices
