"""
Kα Doublet Handling Module Kα 雙峰處理模組
==========================================

Two approaches for handling Cu Kα₁/Kα₂ doublet in XRD data.
處理 XRD 數據中 Cu Kα₁/Kα₂ 雙峰的兩種方法。

1. Ka2Stripper: Remove Kα₂ contribution from spectrum / 從光譜中移除 Kα₂ 貢獨
2. DoubletFitter: Fit both Kα₁ and Kα₂ peaks simultaneously / 同時擬合 Kα₁ 和 Kα₂ 峰

Cu Kα wavelengths 波長 (Bearden 1967, Rev. Mod. Phys. 39, 78):
- Kα₁ = 1.540562 Å (stronger, 2x intensity / 較強)
- Kα₂ = 1.544390 Å (weaker, 1x intensity / 較弱)
- Kα₂/Kα₁ intensity ratio 強度比 ≈ 0.5 (Burger-Dorgelo rule)
"""

import numpy as np
from scipy.optimize import curve_fit, least_squares, minimize
from scipy.ndimage import gaussian_filter1d
from typing import Tuple, Optional, Dict, Any
from dataclasses import dataclass

from .pseudo_voigt import PseudoVoigt, PseudoVoigtParams, TrueVoigt
from axcsas.core.constants import CU_KA1, CU_KA2, KA2_KA1_RATIO


# 波長常數直接使用 CU_KA1/CU_KA2 / Wavelength: use CU_KA1/CU_KA2 directly


@dataclass
class DoubletFitResult:
    """Result of Kα doublet fitting."""
    center_ka1: float       # 2θ position of Kα₁
    center_ka2: float       # 2θ position of Kα₂
    amplitude_ka1: float    # Amplitude of Kα₁
    amplitude_ka2: float    # Amplitude of Kα₂
    fwhm: float            # Shared FWHM
    eta: float             # Shared mixing parameter
    r_squared: float       # Goodness of fit
    success: bool
    message: str
    fwhm_error: float = 0.0    # Standard error of FWHM
    center_error: float = 0.0  # Standard error of center
    eta_error: float = 0.0     # Standard error of eta
    fitted_curve: Optional[np.ndarray] = None
    opt_result: Optional[Any] = None  # Optimization result with Jacobian
    
    @property
    def fwhm_ka1(self) -> float:
        """FWHM of Kα₁ peak (primary result)."""
        return self.fwhm





def theta2_from_wavelength_shift(theta1: float, lambda1: float, lambda2: float) -> float:
    """
    Calculate 2θ shift due to wavelength difference.
    
    Uses Bragg's law: n·λ = 2d·sin(θ)
    For same d-spacing: sin(θ₂)/sin(θ₁) = λ₂/λ₁
    
    Args:
        theta1: 2θ position for wavelength λ₁ (degrees)
        lambda1: First wavelength (Å)
        lambda2: Second wavelength (Å)
        
    Returns:
        2θ position for wavelength λ₂ (degrees)
    """
    theta1_rad = np.radians(theta1 / 2)  # Convert to θ (half of 2θ)
    sin_theta2 = np.sin(theta1_rad) * (lambda2 / lambda1)
    
    # Check if sin value is valid
    if abs(sin_theta2) > 1:
        return theta1  # Return original if invalid
    
    theta2_rad = np.arcsin(sin_theta2)
    return np.degrees(theta2_rad) * 2  # Convert back to 2θ


def calculate_ka2_position(two_theta_ka1: float) -> float:
    """
    Calculate Kα₂ peak position from Kα₁ position.
    
    Args:
        two_theta_ka1: 2θ position of Kα₁ peak (degrees)
        
    Returns:
        2θ position of Kα₂ peak (degrees)
    """
    return theta2_from_wavelength_shift(two_theta_ka1, CU_KA1, CU_KA2)


class Ka2Stripper:
    """
    Remove Kα₂ contribution from XRD spectrum.
    
    The Kα₂ peak is shifted to higher angles and has ~50% intensity of Kα₁.
    This class estimates and subtracts the Kα₂ contribution.
    """
    
    def __init__(self, intensity_ratio: float = KA2_KA1_RATIO):
        """
        Initialize stripper.
        
        Args:
            intensity_ratio: Kα₂/Kα₁ intensity ratio (default: 0.5)
        """
        self.intensity_ratio = intensity_ratio
    
    def strip(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        smooth_sigma: float = 0.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Strip Kα₂ from spectrum using Rachinger correction.
        
        Algorithm:
        1. For each point at 2θ, find the Kα₁ position that would produce Kα₂ here
        2. Subtract intensity_ratio × I(Kα₁ position)
        
        Args:
            two_theta: 2θ array (degrees)
            intensity: Intensity array
            smooth_sigma: Gaussian smoothing sigma (0 = no smoothing)
            
        Returns:
            (two_theta, stripped_intensity) tuple
        """
        stripped = np.copy(intensity).astype(float)
        
        # Apply optional smoothing first
        if smooth_sigma > 0:
            intensity_smooth = gaussian_filter1d(intensity, smooth_sigma)
        else:
            intensity_smooth = intensity
        
        # Process from high to low angles (Kα₂ appears at higher angles)
        for i in range(len(two_theta) - 1, -1, -1):
            theta_ka2 = two_theta[i]
            
            # Calculate where the corresponding Kα₁ would be
            # Kα₂ is at higher angle than Kα₁
            theta_ka1 = theta2_from_wavelength_shift(theta_ka2, CU_KA2, CU_KA1)
            
            # Find the index of Kα₁ position
            idx_ka1 = np.searchsorted(two_theta, theta_ka1)
            
            if 0 <= idx_ka1 < len(two_theta):
                # Interpolate intensity at Kα₁ position
                if idx_ka1 > 0 and idx_ka1 < len(two_theta):
                    # Linear interpolation
                    x1, x2 = two_theta[idx_ka1 - 1], two_theta[idx_ka1]
                    y1, y2 = intensity_smooth[idx_ka1 - 1], intensity_smooth[idx_ka1]
                    if x2 != x1:
                        ka1_intensity = y1 + (y2 - y1) * (theta_ka1 - x1) / (x2 - x1)
                    else:
                        ka1_intensity = y1
                else:
                    ka1_intensity = intensity_smooth[max(0, min(idx_ka1, len(intensity_smooth) - 1))]
                
                # Subtract Kα₂ contribution
                stripped[i] -= self.intensity_ratio * ka1_intensity
        
        # Ensure non-negative
        stripped = np.maximum(stripped, 0)
        
        return two_theta, stripped
    
    def strip_peak_region(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        peak_center: float,
        window: float = 3.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Strip Kα₂ from a specific peak region.
        
        Args:
            two_theta: Full 2θ array
            intensity: Full intensity array
            peak_center: Approximate Kα₁ peak center
            window: Half-width of region to process
            
        Returns:
            (theta_region, stripped_region) tuple
        """
        mask = (two_theta >= peak_center - window) & (two_theta <= peak_center + window + 0.5)
        return self.strip(two_theta[mask], intensity[mask])


class DoubletFitter:
    """
    Fit Kα₁/Kα₂ doublet simultaneously.
    
    Constraints:
    - Kα₂ position calculated from Kα₁ via Bragg's law
    - Kα₂ amplitude = Kα₁ amplitude × 0.5
    - Both peaks share same FWHM and eta
    """
    
    def __init__(
        self,
        intensity_ratio: float = KA2_KA1_RATIO,
        max_iterations: int = 100000,  # Increased for better convergence
        tolerance: float = 1e-8
    ):
        """
        Initialize fitter.
        
        Args:
            intensity_ratio: Fixed Kα₂/Kα₁ intensity ratio
            max_iterations: Maximum fitting iterations
            tolerance: Convergence tolerance
        """
        self.intensity_ratio = intensity_ratio
        self.max_iterations = max_iterations
        self.tolerance = tolerance
    
    def _doublet_profile(
        self,
        x: np.ndarray,
        center_ka1: float,
        amplitude_ka1: float,
        sigma: float,
        gamma: float,
        slope: float,
        intercept: float
    ) -> np.ndarray:
        """
        Calculate Kα₁ + Kα₂ doublet profile using True Voigt.
        
        Args:
            x: 2θ array
            center_ka1: Kα₁ peak center
            amplitude_ka1: Kα₁ amplitude
            sigma: Gaussian width parameter
            gamma: Lorentzian width parameter
            slope, intercept: Linear background
            
        Returns:
            Combined intensity profile
        """
        # Calculate Kα₂ position from Kα₁
        center_ka2 = calculate_ka2_position(center_ka1)
        amplitude_ka2 = amplitude_ka1 * self.intensity_ratio
        
        # Sum of both peaks + background using True Voigt
        ka1 = TrueVoigt.profile(x, center_ka1, amplitude_ka1, sigma, gamma)
        ka2 = TrueVoigt.profile(x, center_ka2, amplitude_ka2, sigma, gamma)
        background = slope * x + intercept
        
        return ka1 + ka2 + background
    
    def fit(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        initial_center: Optional[float] = None,
        initial_fwhm: float = 0.3
    ) -> DoubletFitResult:
        """
        Fit Kα₁/Kα₂ doublet to data using True Voigt model.
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            initial_center: Initial Kα₁ center guess (auto if None)
            initial_fwhm: Initial FWHM guess
            
        Returns:
            DoubletFitResult with fit parameters
        """
        # Auto-detect peak center
        if initial_center is None:
            initial_center = two_theta[np.argmax(intensity)]
        
        peak_amp = np.max(intensity)
        
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
        
        # Improved initial FWHM estimation
        bg_level = min(left_bg, right_bg)
        half_max = (peak_amp + bg_level) / 2
        above_half = intensity > half_max
        if np.any(above_half):
            left_idx = np.argmax(above_half)
            right_idx = len(intensity) - 1 - np.argmax(above_half[::-1])
            if right_idx > left_idx:
                initial_fwhm = max(two_theta[right_idx] - two_theta[left_idx], 0.15)
        
        # Initial guesses for Sigma and Gamma from FWHM
        # Assume initial mix is 50/50 Gaussian/Lorentzian (eta=0.5)
        # fL = 0.5 * fwhm, fG = 0.5 * fwhm
        # sigma = fG / (2*sqrt(2*ln2)) ≈ fG / 2.355
        # gamma = fL / 2
        sigma_init = (0.5 * initial_fwhm) / 2.35482
        gamma_init = (0.5 * initial_fwhm) / 2.0
        
        # Multi-start optimization for robustness
        best_result = None
        best_r_sq = -1
        
        # Define residual function ONCE outside loop
        profile_func = self._doublet_profile
        
        def make_residual(x_data, y_data):
            def residual(params):
                return y_data - profile_func(x_data, *params)
            return residual
        
        residual_func = make_residual(two_theta, intensity)
        
        # Define bounds for [center, amp, sigma, gamma, slope, intercept]
        # sigma, gamma > 0
        lb = np.array([two_theta.min(), 10, 1e-4, 1e-4, -1e6, -1e6])
        ub = np.array([two_theta.max(), peak_amp * 3, 2.0, 2.0, 1e6, 1e6])
        
        # Try multiple initial ratios
        for ratio in [0.3, 0.5, 0.7]:
            # Adjust init params based on ratio (ratio ~ gamma/(sigma+gamma) roughly)
            this_gamma = initial_fwhm * ratio / 2.0
            this_sigma = initial_fwhm * (1-ratio) / 2.35482
            
            x0 = np.array([
                initial_center, 
                peak_amp * 0.8, 
                this_sigma, 
                this_gamma, 
                slope_init, 
                intercept_init
            ])
            
            try:
                # Use least_squares with strict iteration limit
                result = least_squares(
                    residual_func,
                    x0,
                    bounds=(lb, ub),
                    method='trf',
                    ftol=self.tolerance, # Inherited from class init
                    xtol=1e-12,         # High precision
                    gtol=1e-12,
                    max_nfev=self.max_iterations,
                    verbose=0
                )
                
                if result.success or result.status >= 1:
                    popt = result.x
                    fitted = self._doublet_profile(two_theta, *popt)
                    residuals = intensity - fitted
                    chi_sq = np.sum(residuals ** 2)
                    ss_tot = np.sum((intensity - np.mean(intensity)) ** 2)
                    r_sq = 1 - (chi_sq / ss_tot) if ss_tot > 0 else 0
                    
                    if r_sq > best_r_sq:
                        best_r_sq = r_sq
                        
                        # Extract params
                        center_val = popt[0]
                        amp_val = popt[1]
                        sigma_val = popt[2]
                        gamma_val = popt[3]
                        
                        center_ka2 = calculate_ka2_position(center_val)
                        
                        # Convert back to standard reporting metrics
                        # FWHM total
                        fwhm_total = TrueVoigt.fwhm_from_params(sigma_val, gamma_val)
                        
                        # Effective Eta: fL / FWHM_total (Standard approximation)
                        # fL = 2 * gamma
                        fwhm_lorentzian = 2.0 * gamma_val
                        eta_eff = fwhm_lorentzian / fwhm_total if fwhm_total > 0 else 0.5
                        
                        # Calculate standard errors from Jacobian
                        fwhm_err = 0.0
                        center_err = 0.0
                        eta_err = 0.0
                        
                        try:
                            jac = result.jac
                            s_sq = 2 * result.cost / (len(intensity) - len(popt))
                            pcov = s_sq * np.linalg.inv(jac.T @ jac)
                            perr = np.sqrt(np.diag(pcov))
                            
                            # error[0] is center error
                            center_err = perr[0]
                            
                            # Error propagation for FWHM (approximate)
                            sigma_err = perr[2]
                            gamma_err = perr[3]
                            fwhm_err = sigma_err * 2.35 + gamma_err * 2.0
                            
                            # Error propagation for Eta (approximate)
                            # Eta ~ 2*Gamma / FWHM
                            if fwhm_total > 0 and gamma_val > 0:
                                rel_gamma = gamma_err / gamma_val
                                rel_fwhm = fwhm_err / fwhm_total
                                eta_err = eta_eff * (rel_gamma + rel_fwhm)
                            else:
                                eta_err = 0.01
                                
                        except Exception:
                            fwhm_err = fwhm_total * 0.05
                            center_err = 0.001
                            eta_err = 0.05
                            
                        best_result = DoubletFitResult(
                            center_ka1=center_val,
                            center_ka2=center_ka2,
                            amplitude_ka1=amp_val,
                            amplitude_ka2=amp_val * self.intensity_ratio,
                            fwhm=fwhm_total,
                            eta=eta_eff,
                            r_squared=r_sq,
                            success=True,
                            message=f"True Voigt fit converged (R²={r_sq:.4f})",
                            fitted_curve=fitted,
                            fwhm_error=fwhm_err,
                            center_error=center_err,
                            eta_error=eta_err,
                            opt_result=result
                        )
            except Exception:
                continue
        
        if best_result is not None:
            return best_result
        
        # Fallback
        return DoubletFitResult(
            center_ka1=initial_center,
            center_ka2=calculate_ka2_position(initial_center),
            amplitude_ka1=peak_amp,
            amplitude_ka2=peak_amp * self.intensity_ratio,
            fwhm=initial_fwhm,
            eta=0.5,
            r_squared=0,
            success=False,
            message="Multi-start True Voigt failed"
        )


def compare_fitting_methods(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    peak_center: float,
    window: float = 2.5
) -> Dict:
    """
    Compare single peak, Ka2-stripped, and doublet fitting methods.
    
    Args:
        two_theta: Full 2θ array
        intensity: Full intensity array
        peak_center: Expected peak center
        window: Half-width of fitting region
        
    Returns:
        Dict with results from all three methods
    """
    from .lm_optimizer import LMOptimizer
    
    # Select region
    mask = (two_theta >= peak_center - window) & (two_theta <= peak_center + window)
    theta_region = two_theta[mask]
    int_region = intensity[mask]
    
    results = {}
    
    # Method 1: Simple Pseudo-Voigt fit
    optimizer = LMOptimizer()
    pv_result = optimizer.fit_single_peak(theta_region, int_region)
    results['pseudo_voigt'] = {
        'fwhm': pv_result.params.fwhm,
        'eta': pv_result.params.eta,
        'r_squared': pv_result.r_squared,
        'center': pv_result.params.center
    }
    
    # Method 2: Strip Kα₂ then fit
    stripper = Ka2Stripper()
    _, stripped = stripper.strip_peak_region(two_theta, intensity, peak_center, window)
    stripped_result = optimizer.fit_single_peak(theta_region, stripped)
    results['ka2_stripped'] = {
        'fwhm': stripped_result.params.fwhm,
        'eta': stripped_result.params.eta,
        'r_squared': stripped_result.r_squared,
        'center': stripped_result.params.center
    }
    
    # Method 3: Fit doublet directly
    doublet_fitter = DoubletFitter()
    doublet_result = doublet_fitter.fit(theta_region, int_region, peak_center)
    results['doublet'] = {
        'fwhm': doublet_result.fwhm,
        'eta': doublet_result.eta,
        'r_squared': doublet_result.r_squared,
        'center_ka1': doublet_result.center_ka1,
        'center_ka2': doublet_result.center_ka2
    }
    
    return results
