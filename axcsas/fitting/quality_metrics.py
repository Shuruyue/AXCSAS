"""
Fitting Quality Metrics Module
==============================

Quality assessment for XRD peak fitting results.

Reference: 計劃書/03_峰值擬合與參數提取.md §7
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from enum import Enum


class QualityLevel(Enum):
    """Fit quality classification based on R_wp."""
    EXCELLENT = "excellent"     # R_wp ≤ 5%
    GOOD = "good"               # 5% < R_wp ≤ 10%
    LOW_CONFIDENCE = "low_confidence"  # R_wp > 10%
    FAILED = "failed"           # Convergence failure


@dataclass
class FitQualityReport:
    """
    Comprehensive quality report for peak fitting.
    
    Attributes:
        r_wp: Weighted profile R-factor (%)
        r_squared: Coefficient of determination
        rss: Residual sum of squares
        is_valid: True if fit parameters are physically valid
        quality_level: Quality classification
        warnings: List of warning messages
    """
    r_wp: float
    r_squared: float
    rss: float
    is_valid: bool
    quality_level: QualityLevel
    warnings: List[str] = field(default_factory=list)
    
    def summary(self) -> str:
        """Generate human-readable quality summary."""
        status = "✓" if self.is_valid else "✗"
        lines = [
            f"Fit Quality: {self.quality_level.value.upper()} {status}",
            f"  R_wp: {self.r_wp:.2f}%",
            f"  R²: {self.r_squared:.4f}",
            f"  RSS: {self.rss:.2e}",
        ]
        if self.warnings:
            lines.append("  Warnings:")
            for w in self.warnings:
                lines.append(f"    - {w}")
        return "\n".join(lines)


# =============================================================================
# R_wp Calculation (文件 03 §7.1)
# =============================================================================

def calculate_r_wp(
    observed: np.ndarray,
    calculated: np.ndarray,
    weights: Optional[np.ndarray] = None
) -> float:
    """
    Calculate Weighted Profile R-factor (R_wp).
    
    Reference: 文件 03 §7.1
    
    R_wp = sqrt(Σ w_i (I_obs - I_calc)² / Σ w_i I_obs²) × 100%
    
    Quality thresholds:
      - ≤ 5%: Excellent
      - 5-10%: Good (acceptable)
      - > 10%: Low confidence (requires review)
    
    Args:
        observed: Observed intensity array
        calculated: Calculated (fitted) intensity array
        weights: Optional weights (default: Poisson weights 1/I_obs)
        
    Returns:
        R_wp value in percent (%)
    """
    if weights is None:
        # Poisson weights: w = 1/I for counting statistics
        weights = 1.0 / np.maximum(observed, 1.0)
    
    residuals = observed - calculated
    numerator = np.sum(weights * residuals**2)
    denominator = np.sum(weights * observed**2)
    
    if denominator == 0:
        return float('inf')
    
    return 100.0 * np.sqrt(numerator / denominator)


def calculate_rss(
    observed: np.ndarray,
    calculated: np.ndarray
) -> float:
    """
    Calculate Residual Sum of Squares (RSS).
    
    RSS = Σ (I_obs - I_calc)²
    
    Args:
        observed: Observed intensity array
        calculated: Calculated (fitted) intensity array
        
    Returns:
        RSS value
    """
    return float(np.sum((observed - calculated)**2))


def calculate_r_squared(
    observed: np.ndarray,
    calculated: np.ndarray
) -> float:
    """
    Calculate coefficient of determination (R²).
    
    R² = 1 - SS_res / SS_tot
    
    Args:
        observed: Observed intensity array
        calculated: Calculated (fitted) intensity array
        
    Returns:
        R² value (0-1, higher is better)
    """
    ss_res = np.sum((observed - calculated)**2)
    ss_tot = np.sum((observed - np.mean(observed))**2)
    
    if ss_tot == 0:
        return 0.0
    
    return 1.0 - (ss_res / ss_tot)


# =============================================================================
# Fit Validation (文件 03 §7.2)
# =============================================================================

def validate_fit_parameters(
    center: float,
    amplitude: float,
    fwhm: float,
    eta: float,
    two_theta_range: Tuple[float, float] = (10, 150)
) -> Tuple[bool, List[str]]:
    """
    Validate fitted parameters for physical reasonableness.
    
    Reference: 文件 03 §7.2 / §2.2
    
    Physical constraints:
      - η must be in [0, 1]
      - FWHM must be > 0
      - amplitude must be > 0
      - center must be in valid 2θ range
    
    Args:
        center: Peak center position (degrees)
        amplitude: Peak amplitude
        fwhm: Full width at half maximum (degrees)
        eta: Mixing parameter
        two_theta_range: Valid 2θ range
        
    Returns:
        Tuple of (is_valid, list_of_warnings)
    """
    warnings = []
    is_valid = True
    
    # Check η bounds
    if eta < 0 or eta > 1:
        warnings.append(f"η = {eta:.3f} out of bounds [0,1] - convergence failure")
        is_valid = False
    
    # Check FWHM positive
    if fwhm <= 0:
        warnings.append(f"FWHM = {fwhm:.4f}° ≤ 0 - physically invalid")
        is_valid = False
    elif fwhm > 5.0:
        warnings.append(f"FWHM = {fwhm:.3f}° unusually large - check data")
    
    # Check amplitude positive
    if amplitude <= 0:
        warnings.append(f"Amplitude = {amplitude:.1f} ≤ 0 - invalid")
        is_valid = False
    
    # Check center in range
    if center < two_theta_range[0] or center > two_theta_range[1]:
        warnings.append(
            f"Center = {center:.2f}° outside expected range "
            f"[{two_theta_range[0]}, {two_theta_range[1]}]"
        )
    
    return is_valid, warnings


def generate_quality_report(
    observed: np.ndarray,
    calculated: np.ndarray,
    center: float,
    amplitude: float,
    fwhm: float,
    eta: float
) -> FitQualityReport:
    """
    Generate comprehensive quality report for a fitted peak.
    
    Args:
        observed: Observed intensity array
        calculated: Fitted intensity array
        center, amplitude, fwhm, eta: Fitted parameters
        
    Returns:
        FitQualityReport with all quality metrics
    """
    # Calculate quality metrics
    r_wp = calculate_r_wp(observed, calculated)
    r_squared = calculate_r_squared(observed, calculated)
    rss = calculate_rss(observed, calculated)
    
    # Validate parameters
    is_valid, param_warnings = validate_fit_parameters(
        center, amplitude, fwhm, eta
    )
    
    # Determine quality level
    if not is_valid:
        quality_level = QualityLevel.FAILED
    elif r_wp <= 5.0:
        quality_level = QualityLevel.EXCELLENT
    elif r_wp <= 10.0:
        quality_level = QualityLevel.GOOD
    else:
        quality_level = QualityLevel.LOW_CONFIDENCE
    
    # Collect all warnings
    warnings = param_warnings.copy()
    
    if r_wp > 10.0 and is_valid:
        warnings.append(f"R_wp = {r_wp:.1f}% exceeds 10% threshold")
    
    return FitQualityReport(
        r_wp=r_wp,
        r_squared=r_squared,
        rss=rss,
        is_valid=is_valid,
        quality_level=quality_level,
        warnings=warnings
    )


def get_quality_threshold() -> dict:
    """Return quality threshold constants."""
    return {
        "r_wp_excellent": 5.0,
        "r_wp_good": 10.0,
        "fwhm_max": 5.0,
        "eta_min": 0.0,
        "eta_max": 1.0,
    }
