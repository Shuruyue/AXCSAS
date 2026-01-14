"""
AXCSAS Fitting API
==================

Centralized interface for peak fitting operations across the project.
Ensures consistency between:
1. Analysis Pipeline (Scherrer calculation)
2. Diagnostic Plots (Visual verification)
3. Trend Plots (FWHM evolution)

This module enforces the use of **Kα Doublet Fitting** as the standard for 
crystallite size analysis, correcting the previous "Single Peak" error 
which overestimated FWHM.
"""

from dataclasses import dataclass
from typing import Dict, Any, Optional, Tuple, NamedTuple
import numpy as np
from pathlib import Path

# Import core fitting logic
# Note: We wrap the lower-level classes to provide a simplified API
from axcsas.fitting.pseudo_voigt import PseudoVoigt, PseudoVoigtParams
from axcsas.fitting.ka_doublet import DoubletFitter, calculate_ka2_position
from axcsas.fitting.lm_optimizer import LMOptimizer

@dataclass
class FittingResult:
    """Standardized result from any fitting method."""
    success: bool
    method: str  # 'doublet', 'single-pv', 'enhanced-pv'
    
    # Primary peak (Kα1) parameters
    center: float
    fwhm: float
    amplitude: float
    eta: float  # Gaussian-Lorentzian mixing (0=Gaussian, 1=Lorentzian)
    
    # Goodness of fit
    r_squared: float
    chi2_red: float = 1.0  # Reduced Chi-squared
    
    # Uncertainty (standard errors)
    center_err: float = 0.0
    fwhm_err: float = 0.0
    eta_err: float = 0.0
    
    # For plotting
    fitted_curve: Optional[np.ndarray] = None
    theta_range: Optional[np.ndarray] = None
    intensity_range: Optional[np.ndarray] = None
    
    # Doublet specific
    center_ka2: Optional[float] = None


def fit_peak_doublet(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    expected_center: float,
    window: float = 2.5
) -> FittingResult:
    """
    Fit a peak using the Kα Doublet Model (PHYSICALLY CORRECT).
    
    This splits the peak into Kα1 and Kα2 components based on Bragg's law 
    and fixed intensity ratio (approx 2:1).
    
    Args:
        two_theta: X-axis data
        intensity: Y-axis data
        expected_center: Approximate position of Kα1 peak
        window: Data range to use around center (+/- degrees)
        
    Returns:
        FittingResult with Kα1 parameters (FWHM is for Kα1 only).
    """
    # 1. Select data range
    mask = (two_theta >= expected_center - window) & (two_theta <= expected_center + window)
    if not np.any(mask):
        return _make_failed_result("doublet")
    
    theta_range = two_theta[mask]
    int_range = intensity[mask]
    
    # 2. Check signal quality
    if np.max(int_range) < 50:
        return _make_failed_result("doublet", theta_range, int_range)

    # 3. Perform Fitting
    try:
        fitter = DoubletFitter(max_iterations=5000)
        
        # Estimate initial FWHM
        peak_idx = np.argmax(int_range)
        peak_int = int_range[peak_idx]
        half_max = peak_int / 2
        initial_fwhm = 0.3 # Default guess
        
        # Simple width estimation
        l_idx = peak_idx
        while l_idx > 0 and int_range[l_idx] > half_max: l_idx -= 1
        r_idx = peak_idx
        while r_idx < len(int_range)-1 and int_range[r_idx] > half_max: r_idx += 1
        width_guess = theta_range[r_idx] - theta_range[l_idx]
        if width_guess > 0.05:
            initial_fwhm = width_guess

        res = fitter.fit(theta_range, int_range, expected_center, initial_fwhm)
        
        if res.success and res.r_squared > 0.6: # Relaxed threshold for difficult peaks
            return FittingResult(
                success=True,
                method='doublet',
                center=res.center_ka1,
                fwhm=res.fwhm, # This is Kα1 FWHM
                amplitude=res.amplitude_ka1,
                eta=res.eta,
                r_squared=res.r_squared,
                chi2_red=1.0, # Not currently calculated by DoubletFitter, placeholder
                center_ka2=res.center_ka2,
                theta_range=theta_range,
                intensity_range=int_range,
                fitted_curve=res.fitted_curve
            )
            
    except Exception as e:
        # print(f"Doublet fitting error for {expected_center}: {e}")
        pass
        
    return _make_failed_result("doublet", theta_range, int_range)


def fit_peak_simple(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    expected_center: float,
    window: float = 2.5
) -> FittingResult:
    """
    Fit a peak using Single Pseudo-Voigt (LEGACY/INACCURATE).
    
    Included only for comparison. Treating a doublet as a single peak 
    results in virtually inflated FWHM.
    """
    mask = (two_theta >= expected_center - window) & (two_theta <= expected_center + window)
    if not np.any(mask):
        return _make_failed_result("single-pv")
        
    theta_range = two_theta[mask]
    int_range = intensity[mask]
    
    try:
        optimizer = LMOptimizer()
        res = optimizer.fit_single_peak(theta_range, int_range)
        
        if res.success:
            fitted_curve = PseudoVoigt.profile(
                theta_range, res.params.center, res.params.amplitude, 
                res.params.fwhm, res.params.eta
            ) + 0 # Simple optimizer usually assumes constant background removed or handled
            
            return FittingResult(
                success=True,
                method='single-pv',
                center=res.params.center,
                fwhm=res.params.fwhm,
                amplitude=res.params.amplitude,
                eta=res.params.eta,
                r_squared=res.r_squared,
                theta_range=theta_range,
                intensity_range=int_range,
                fitted_curve=fitted_curve
            )
    except Exception:
        pass
        
    return _make_failed_result("single-pv", theta_range, int_range)


def _make_failed_result(method: str, theta=None, intensity=None):
    """Helper to create failed result."""
    return FittingResult(
        success=False,
        method=method,
        center=np.nan,
        fwhm=np.nan,
        amplitude=np.nan,
        eta=np.nan,
        r_squared=0.0,
        theta_range=theta,
        intensity_range=intensity
    )
