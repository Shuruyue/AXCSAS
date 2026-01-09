"""
Levenberg-Marquardt Optimizer Module
Implements non-linear least squares fitting for XRD peaks.
"""

import numpy as np
from scipy.optimize import curve_fit, least_squares
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass, field

from .pseudo_voigt import PseudoVoigt, PseudoVoigtParams


@dataclass
class FitResult:
    """
    Result of peak fitting.
    
    Enhanced with R_wp quality metric, integrated area, and hkl assignment.
    Reference: 文件 03 §6-7
    """
    params: PseudoVoigtParams
    covariance: Optional[np.ndarray]
    residuals: np.ndarray
    chi_squared: float
    r_squared: float
    success: bool
    message: str
    # Phase 03 enhancements
    r_wp: float = 0.0              # Weighted Profile R-factor (%)
    area: float = 0.0              # Integrated intensity (peak area)
    hkl: Optional[Tuple[int, int, int]] = None  # Miller indices
    
    def is_valid(self) -> bool:
        """Check if fit result is physically valid."""
        if not self.success:
            return False
        p = self.params
        return 0 <= p.eta <= 1 and p.fwhm > 0 and p.amplitude > 0


class LMOptimizer:
    """
    Levenberg-Marquardt optimizer for Pseudo-Voigt peak fitting.
    
    Minimizes the residual sum of squares (RSS):
    RSS = Σ(I_obs - I_calc)²
    """
    
    def __init__(
        self,
        max_iterations: int = 1000,
        tolerance: float = 1e-8
    ):
        """
        Initialize optimizer.
        
        Args:
            max_iterations: Maximum number of iterations
            tolerance: Convergence tolerance
        """
        self.max_iterations = max_iterations
        self.tolerance = tolerance
    
    def fit_single_peak(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        initial_guess: Optional[PseudoVoigtParams] = None,
        peak_idx: Optional[int] = None
    ) -> FitResult:
        """
        Fit a single Pseudo-Voigt peak.
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            initial_guess: Initial parameter guess
            peak_idx: Peak index for automatic initial guess
            
        Returns:
            FitResult object
        """
        # Generate initial guess if not provided
        if initial_guess is None:
            if peak_idx is not None:
                initial_guess = self._auto_initial_guess(
                    two_theta, intensity, peak_idx
                )
            else:
                # Use maximum as peak
                peak_idx = np.argmax(intensity)
                initial_guess = self._auto_initial_guess(
                    two_theta, intensity, peak_idx
                )
        
        p0 = initial_guess.to_array()
        
        # Parameter bounds
        bounds = (
            [two_theta.min(), 0, 0.01, 0],      # Lower bounds
            [two_theta.max(), np.inf, 5.0, 1]   # Upper bounds
        )
        
        try:
            popt, pcov = curve_fit(
                lambda x, c, a, w, e: PseudoVoigt.profile(x, c, a, w, e),
                two_theta,
                intensity,
                p0=p0,
                bounds=bounds,
                maxfev=self.max_iterations,
                ftol=self.tolerance
            )
            
            # Calculate fit quality metrics
            fitted = PseudoVoigt.profile(two_theta, *popt)
            residuals = intensity - fitted
            
            chi_sq = np.sum(residuals ** 2)
            ss_tot = np.sum((intensity - np.mean(intensity)) ** 2)
            r_sq = 1 - (chi_sq / ss_tot) if ss_tot > 0 else 0
            
            return FitResult(
                params=PseudoVoigtParams.from_array(popt),
                covariance=pcov,
                residuals=residuals,
                chi_squared=chi_sq,
                r_squared=r_sq,
                success=True,
                message="Fit converged"
            )
            
        except Exception as e:
            return FitResult(
                params=initial_guess,
                covariance=None,
                residuals=np.zeros_like(intensity),
                chi_squared=float('inf'),
                r_squared=0,
                success=False,
                message=str(e)
            )
    
    def fit_multi_peak(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        n_peaks: int,
        initial_guesses: Optional[List[PseudoVoigtParams]] = None
    ) -> List[FitResult]:
        """
        Fit multiple peaks simultaneously.
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            n_peaks: Number of peaks to fit
            initial_guesses: List of initial guesses
            
        Returns:
            List of FitResult objects
        """
        # For now, fit peaks sequentially
        # TODO: Implement simultaneous multi-peak fitting
        results = []
        
        if initial_guesses:
            for guess in initial_guesses:
                # Define region around peak
                center = guess.center
                margin = 2.0  # degrees
                mask = (two_theta >= center - margin) & (two_theta <= center + margin)
                
                result = self.fit_single_peak(
                    two_theta[mask],
                    intensity[mask],
                    initial_guess=guess
                )
                results.append(result)
        
        return results
    
    def _auto_initial_guess(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        peak_idx: int
    ) -> PseudoVoigtParams:
        """
        Generate automatic initial guess for peak fitting.
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            peak_idx: Peak index
            
        Returns:
            PseudoVoigtParams with initial guess
        """
        center = two_theta[peak_idx]
        amplitude = intensity[peak_idx]
        
        # Estimate FWHM
        half_max = amplitude / 2
        left_idx = peak_idx
        while left_idx > 0 and intensity[left_idx] > half_max:
            left_idx -= 1
        right_idx = peak_idx
        while right_idx < len(intensity) - 1 and intensity[right_idx] > half_max:
            right_idx += 1
        fwhm = max(two_theta[right_idx] - two_theta[left_idx], 0.1)
        
        # Default eta = 0.5 (equal Gaussian/Lorentzian mixing)
        eta = 0.5
        
        return PseudoVoigtParams(center, amplitude, fwhm, eta)


def fit_peaks(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    peak_positions: Optional[List[float]] = None
) -> List[FitResult]:
    """
    Convenience function for peak fitting.
    
    Args:
        two_theta: 2θ array
        intensity: Intensity array
        peak_positions: Optional list of peak positions
        
    Returns:
        List of FitResult objects
    """
    optimizer = LMOptimizer()
    
    if peak_positions is None:
        # Fit single dominant peak
        return [optimizer.fit_single_peak(two_theta, intensity)]
    
    # Fit each specified peak
    results = []
    for pos in peak_positions:
        # Find closest index
        idx = np.argmin(np.abs(two_theta - pos))
        
        # Define local region
        margin = 2.0
        mask = (two_theta >= pos - margin) & (two_theta <= pos + margin)
        
        result = optimizer.fit_single_peak(
            two_theta[mask],
            intensity[mask],
            peak_idx=np.argmin(np.abs(two_theta[mask] - pos))
        )
        results.append(result)
    
    return results
