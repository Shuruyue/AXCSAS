"""
Levenberg-Marquardt Optimizer Module
Implements non-linear least squares fitting for XRD peaks.
"""

import numpy as np
from scipy.optimize import curve_fit, least_squares
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass, field

from .pseudo_voigt import PseudoVoigt, PseudoVoigtParams, TrueVoigt, VoigtParams


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
        peak_idx: Optional[int] = None,
        use_linear_background: bool = True
    ) -> FitResult:
        """
        Fit a single Pseudo-Voigt peak with linear background.
        
        Enhanced model: I(2θ) = PV(2θ) + (slope × 2θ) + intercept
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            initial_guess: Initial parameter guess
            peak_idx: Peak index for automatic initial guess
            use_linear_background: If True, use linear bg; else constant
            
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
        
        # Estimate background from edges of data
        n_edge = max(5, len(intensity) // 8)
        left_bg = np.mean(intensity[:n_edge])
        right_bg = np.mean(intensity[-n_edge:])
        left_x = np.mean(two_theta[:n_edge])
        right_x = np.mean(two_theta[-n_edge:])
        
        # Linear background: y = slope * x + intercept
        if right_x != left_x:
            slope_estimate = (right_bg - left_bg) / (right_x - left_x)
        else:
            slope_estimate = 0.0
        intercept_estimate = left_bg - slope_estimate * left_x
        
        # Subtract background estimate from amplitude
        center_x = initial_guess.center
        bg_at_center = slope_estimate * center_x + intercept_estimate
        amplitude_above_bg = max(initial_guess.amplitude - bg_at_center, 100)
        
        if use_linear_background:
            # Parameters: [center, amplitude, fwhm, eta, slope, intercept]
            p0 = [
                initial_guess.center,
                amplitude_above_bg,
                initial_guess.fwhm,
                initial_guess.eta,
                slope_estimate,
                intercept_estimate
            ]
            
            # Parameter bounds
            bounds = (
                [two_theta.min(), 10, 0.02, 0, -np.inf, 0],        # Lower bounds
                [two_theta.max(), np.inf, 5.0, 1, np.inf, np.inf]  # Upper bounds
            )
            
            # Fitting function with linear background
            def pv_with_linear_bg(x, c, a, w, e, slope, intercept):
                return PseudoVoigt.profile(x, c, a, w, e) + slope * x + intercept
            
            fit_func = pv_with_linear_bg
        else:
            # Parameters: [center, amplitude, fwhm, eta, background]
            p0 = [
                initial_guess.center,
                amplitude_above_bg,
                initial_guess.fwhm,
                initial_guess.eta,
                max(intercept_estimate, 0)
            ]
            
            bounds = (
                [two_theta.min(), 10, 0.02, 0, 0],
                [two_theta.max(), np.inf, 5.0, 1, np.inf]
            )
            
            def pv_with_const_bg(x, c, a, w, e, bg):
                return PseudoVoigt.profile(x, c, a, w, e) + bg
            
            fit_func = pv_with_const_bg
        
        try:
            popt, pcov = curve_fit(
                fit_func,
                two_theta,
                intensity,
                p0=p0,
                bounds=bounds,
                maxfev=self.max_iterations * 2,  # More iterations
                ftol=self.tolerance,
                method='trf'  # Trust Region Reflective for better convergence
            )
            
            # Calculate fit quality metrics
            fitted = fit_func(two_theta, *popt)
            residuals = intensity - fitted
            
            chi_sq = np.sum(residuals ** 2)
            ss_tot = np.sum((intensity - np.mean(intensity)) ** 2)
            r_sq = 1 - (chi_sq / ss_tot) if ss_tot > 0 else 0
            
            # Extract Pseudo-Voigt parameters
            params = PseudoVoigtParams(
                center=popt[0],
                amplitude=popt[1],
                fwhm=popt[2],
                eta=popt[3]
            )
            
            if use_linear_background:
                bg_info = f"linear bg: {popt[4]:.2f}x + {popt[5]:.1f}"
            else:
                bg_info = f"const bg: {popt[4]:.1f}"
            
            return FitResult(
                params=params,
                covariance=pcov,
                residuals=residuals,
                chi_squared=chi_sq,
                r_squared=r_sq,
                success=True,
                message=f"Fit converged ({bg_info})"
            )
            
        except Exception as e:
            # Try fallback with simpler model
            return self._fallback_fit(two_theta, intensity, initial_guess, str(e))
    
    def _fallback_fit(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        initial_guess: PseudoVoigtParams,
        original_error: str
    ) -> FitResult:
        """Fallback fitting with simpler constant background model."""
        try:
            # Simple constant background
            bg = np.min(intensity)
            
            def simple_pv(x, c, a, w, e):
                return PseudoVoigt.profile(x, c, a, w, e) + bg
            
            p0 = [initial_guess.center, initial_guess.amplitude - bg, 
                  initial_guess.fwhm, initial_guess.eta]
            bounds = ([two_theta.min(), 0, 0.01, 0], 
                     [two_theta.max(), np.inf, 5.0, 1])
            
            popt, pcov = curve_fit(simple_pv, two_theta, intensity, p0=p0, 
                                   bounds=bounds, maxfev=500)
            
            fitted = simple_pv(two_theta, *popt)
            residuals = intensity - fitted
            chi_sq = np.sum(residuals ** 2)
            ss_tot = np.sum((intensity - np.mean(intensity)) ** 2)
            r_sq = 1 - (chi_sq / ss_tot) if ss_tot > 0 else 0
            
            return FitResult(
                params=PseudoVoigtParams.from_array(popt),
                covariance=pcov, residuals=residuals, chi_squared=chi_sq,
                r_squared=r_sq, success=True,
                message=f"Fallback fit (const bg={bg:.1f})"
            )
        except Exception as e2:
            return FitResult(
                params=initial_guess, covariance=None,
                residuals=np.zeros_like(intensity), chi_squared=float('inf'),
                r_squared=0, success=False,
                message=f"Original: {original_error}; Fallback: {str(e2)}"
            )
    
    def fit_voigt(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        initial_fwhm: float = 0.3,
        initial_eta: float = 0.5
    ) -> FitResult:
        """
        Fit using true Voigt profile (convolution of Gaussian and Lorentzian).
        
        This is more physically accurate than Pseudo-Voigt.
        Parameters: [center, amplitude, sigma, gamma, slope, intercept]
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            initial_fwhm: Initial FWHM estimate
            initial_eta: Initial Lorentzian fraction estimate
            
        Returns:
            FitResult with FWHM calculated from sigma and gamma
        """
        # Find peak
        peak_idx = np.argmax(intensity)
        center_init = two_theta[peak_idx]
        amplitude_init = intensity[peak_idx]
        
        # Estimate background
        n_edge = max(5, len(intensity) // 8)
        left_bg = np.mean(intensity[:n_edge])
        right_bg = np.mean(intensity[-n_edge:])
        left_x = np.mean(two_theta[:n_edge])
        right_x = np.mean(two_theta[-n_edge:])
        
        if right_x != left_x:
            slope_init = (right_bg - left_bg) / (right_x - left_x)
        else:
            slope_init = 0.0
        intercept_init = left_bg - slope_init * left_x
        
        # Convert FWHM and eta to sigma and gamma
        sigma_init, gamma_init = TrueVoigt.params_from_fwhm(initial_fwhm, initial_eta)
        
        # Parameters: [center, amplitude, sigma, gamma, slope, intercept]
        p0 = [center_init, amplitude_init - (slope_init * center_init + intercept_init),
              sigma_init, gamma_init, slope_init, intercept_init]
        
        bounds = (
            [two_theta.min(), 10, 0.001, 0.001, -np.inf, 0],
            [two_theta.max(), np.inf, 2.0, 2.0, np.inf, np.inf]
        )
        
        def voigt_with_bg(x, c, a, sigma, gamma, slope, intercept):
            return TrueVoigt.profile(x, c, a, sigma, gamma) + slope * x + intercept
        
        try:
            popt, pcov = curve_fit(
                voigt_with_bg,
                two_theta,
                intensity,
                p0=p0,
                bounds=bounds,
                maxfev=self.max_iterations * 3,
                ftol=self.tolerance,
                method='trf'
            )
            
            fitted = voigt_with_bg(two_theta, *popt)
            residuals = intensity - fitted
            chi_sq = np.sum(residuals ** 2)
            ss_tot = np.sum((intensity - np.mean(intensity)) ** 2)
            r_sq = 1 - (chi_sq / ss_tot) if ss_tot > 0 else 0
            
            # Calculate total FWHM from sigma and gamma
            fwhm_total = TrueVoigt.fwhm_from_params(popt[2], popt[3])
            
            # Calculate eta equivalent (Lorentzian fraction)
            fG = 2.0 * popt[2] * np.sqrt(2.0 * np.log(2.0))
            fL = 2.0 * popt[3]
            eta_equiv = fL / (fG + fL) if (fG + fL) > 0 else 0.5
            
            # Store as PseudoVoigtParams for compatibility
            params = PseudoVoigtParams(
                center=popt[0],
                amplitude=popt[1],
                fwhm=fwhm_total,
                eta=eta_equiv
            )
            
            return FitResult(
                params=params,
                covariance=pcov,
                residuals=residuals,
                chi_squared=chi_sq,
                r_squared=r_sq,
                success=True,
                message=f"Voigt fit (σ={popt[2]:.4f}, γ={popt[3]:.4f})"
            )
            
        except Exception as e:
            # Fallback to Pseudo-Voigt
            return self.fit_single_peak(two_theta, intensity)
    
    def fit_multi_peak(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        n_peaks: int,
        initial_guesses: Optional[List[PseudoVoigtParams]] = None
    ) -> List[FitResult]:
        """
        Fit multiple peaks sequentially.
        依序擬合多個峰。
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            n_peaks: Number of peaks to fit
            initial_guesses: List of initial guesses
            
        Returns:
            List of FitResult objects
            
        Note / 注意:
            目前使用依序擬合（每個峰獨立擬合），而非同時擬合。
            對於嚴重重疊的峰（雖然在正常銅樣品中少見），
            同時擬合的準確度可能更高。
            
            Currently uses sequential fitting (each peak fitted independently),
            not simultaneous optimization. For severely overlapping peaks
            (rare in normal copper samples), simultaneous fitting may be more accurate.
            
            Future improvement: implement global multi-peak optimization.
        """
        # 演算法限制：依序擬合 - 列入未來優化項目
        # Algorithm limitation: sequential fitting - future optimization target
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
