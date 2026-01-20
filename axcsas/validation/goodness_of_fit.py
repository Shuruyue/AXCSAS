"""
Goodness of Fit Module 擬合優度模組
====================================
Implements quality metrics for XRD peak fitting.
實現 XRD 峰擬合的品質指標。
"""

import numpy as np
from typing import Optional
from dataclasses import dataclass


@dataclass
class FitQualityResult:
    """Result of fit quality assessment."""
    rwp: float              # Weighted Profile R-factor (%)
    rp: float               # Profile R-factor (%)
    r_squared: float        # Coefficient of determination
    chi_squared: float      # Chi-squared statistic
    reduced_chi_squared: float  # Reduced chi-squared
    is_acceptable: bool     # Whether fit meets quality threshold
    
    def __repr__(self):
        status = "✓" if self.is_acceptable else "✗"
        return f"FitQuality({status} Rwp={self.rwp:.2f}%, R²={self.r_squared:.4f})"


def calculate_rwp(
    observed: np.ndarray,
    calculated: np.ndarray,
    weights: Optional[np.ndarray] = None
) -> float:
    """
    Calculate weighted profile R-factor (Rwp).
    
    Rwp = √[Σ w_i (y_obs - y_calc)² / Σ w_i y_obs²] × 100%
    
    Args:
        observed: Observed intensity array
        calculated: Calculated (fitted) intensity array
        weights: Optional weight array. If None, uses 1/y_obs
        
    Returns:
        Rwp value in percent
    """
    observed = np.asarray(observed)
    calculated = np.asarray(calculated)
    
    if len(observed) != len(calculated):
        raise ValueError("observed and calculated must have same length")
    
    if weights is None:
        # Standard weighting: w = 1/y_obs (avoid division by zero)
        weights = 1.0 / np.maximum(observed, 1.0)
    
    residual_sq = (observed - calculated) ** 2
    
    numerator = np.sum(weights * residual_sq)
    denominator = np.sum(weights * observed ** 2)
    
    if denominator < 1e-10:
        return 100.0  # Maximum error
    
    rwp = np.sqrt(numerator / denominator) * 100
    return rwp


def calculate_rp(
    observed: np.ndarray,
    calculated: np.ndarray
) -> float:
    """
    Calculate profile R-factor (Rp).
    
    Rp = Σ |y_obs - y_calc| / Σ y_obs × 100%
    
    Args:
        observed: Observed intensity array
        calculated: Calculated intensity array
        
    Returns:
        Rp value in percent
    """
    observed = np.asarray(observed)
    calculated = np.asarray(calculated)
    
    numerator = np.sum(np.abs(observed - calculated))
    denominator = np.sum(observed)
    
    if denominator < 1e-10:
        return 100.0
    
    return (numerator / denominator) * 100


def calculate_r_squared(
    observed: np.ndarray,
    calculated: np.ndarray
) -> float:
    """
    Calculate coefficient of determination (R²).
    
    R² = 1 - SS_res / SS_tot
    
    where:
    - SS_res = Σ (y_obs - y_calc)²
    - SS_tot = Σ (y_obs - y_mean)²
    
    Args:
        observed: Observed intensity array
        calculated: Calculated intensity array
        
    Returns:
        R² value (0 to 1, higher is better)
    """
    observed = np.asarray(observed)
    calculated = np.asarray(calculated)
    
    y_mean = np.mean(observed)
    
    ss_res = np.sum((observed - calculated) ** 2)
    ss_tot = np.sum((observed - y_mean) ** 2)
    
    if ss_tot < 1e-10:
        return 0.0
    
    return 1 - (ss_res / ss_tot)


def calculate_chi_squared(
    observed: np.ndarray,
    calculated: np.ndarray,
    weights: Optional[np.ndarray] = None,
    n_params: int = 0
) -> tuple:
    """
    Calculate chi-squared statistic.
    
    χ² = Σ w_i (y_obs - y_calc)²
    χ²_reduced = χ² / (N - n_params)
    
    Args:
        observed: Observed intensity array
        calculated: Calculated intensity array
        weights: Optional weights (default: 1/y_obs)
        n_params: Number of fitted parameters
        
    Returns:
        Tuple of (chi_squared, reduced_chi_squared)
    """
    observed = np.asarray(observed)
    calculated = np.asarray(calculated)
    
    if weights is None:
        weights = 1.0 / np.maximum(observed, 1.0)
    
    chi_sq = np.sum(weights * (observed - calculated) ** 2)
    
    dof = len(observed) - n_params
    if dof <= 0:
        dof = 1
    
    reduced_chi_sq = chi_sq / dof
    
    return chi_sq, reduced_chi_sq


def assess_fit_quality(
    observed: np.ndarray,
    calculated: np.ndarray,
    weights: Optional[np.ndarray] = None,
    n_params: int = 0,
    rwp_threshold: float = 10.0
) -> FitQualityResult:
    """
    Comprehensive fit quality assessment.
    
    Args:
        observed: Observed intensity array
        calculated: Calculated intensity array
        weights: Optional weight array
        n_params: Number of fitted parameters
        rwp_threshold: Maximum acceptable Rwp (%)
        
    Returns:
        FitQualityResult with all metrics
    """
    rwp = calculate_rwp(observed, calculated, weights)
    rp = calculate_rp(observed, calculated)
    r_squared = calculate_r_squared(observed, calculated)
    chi_sq, reduced_chi_sq = calculate_chi_squared(
        observed, calculated, weights, n_params
    )
    
    is_acceptable = rwp <= rwp_threshold
    
    return FitQualityResult(
        rwp=rwp,
        rp=rp,
        r_squared=r_squared,
        chi_squared=chi_sq,
        reduced_chi_squared=reduced_chi_sq,
        is_acceptable=is_acceptable
    )
