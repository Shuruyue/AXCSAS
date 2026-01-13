"""
Kα Doublet Handling Module
=========================

Two approaches for handling Cu Kα₁/Kα₂ doublet in XRD data:

1. Ka2Stripper: Remove Kα₂ contribution from spectrum
2. DoubletFitter: Fit both Kα₁ and Kα₂ peaks simultaneously

Cu Kα wavelengths:
- Kα₁ = 1.54056 Å (stronger, 2x intensity)
- Kα₂ = 1.54439 Å (weaker, 1x intensity)
- Kα₂/Kα₁ intensity ratio ≈ 0.5
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from typing import Tuple, Optional, Dict
from dataclasses import dataclass

from .pseudo_voigt import PseudoVoigt, PseudoVoigtParams, TrueVoigt


# Cu Kα wavelengths in Angstrom
LAMBDA_KA1 = 1.54056
LAMBDA_KA2 = 1.54439
KA2_KA1_RATIO = 0.5  # Intensity ratio Kα₂/Kα₁


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
    fitted_curve: Optional[np.ndarray] = None
    
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
    return theta2_from_wavelength_shift(two_theta_ka1, LAMBDA_KA1, LAMBDA_KA2)


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
            theta_ka1 = theta2_from_wavelength_shift(theta_ka2, LAMBDA_KA2, LAMBDA_KA1)
            
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
        max_iterations: int = 2000,
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
        fwhm: float,
        eta: float,
        slope: float,
        intercept: float
    ) -> np.ndarray:
        """
        Calculate Kα₁ + Kα₂ doublet profile.
        
        Args:
            x: 2θ array
            center_ka1: Kα₁ peak center
            amplitude_ka1: Kα₁ amplitude
            fwhm: Shared FWHM
            eta: Shared mixing parameter
            slope, intercept: Linear background
            
        Returns:
            Combined intensity profile
        """
        # Calculate Kα₂ position from Kα₁
        center_ka2 = calculate_ka2_position(center_ka1)
        amplitude_ka2 = amplitude_ka1 * self.intensity_ratio
        
        # Sum of both peaks + background
        ka1 = PseudoVoigt.profile(x, center_ka1, amplitude_ka1, fwhm, eta)
        ka2 = PseudoVoigt.profile(x, center_ka2, amplitude_ka2, fwhm, eta)
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
        Fit Kα₁/Kα₂ doublet to data.
        
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
        
        # Initial parameters: [center_ka1, amplitude_ka1, fwhm, eta, slope, intercept]
        p0 = [initial_center, peak_amp * 0.8, initial_fwhm, 0.5, slope_init, intercept_init]
        
        bounds = (
            [two_theta.min(), 10, 0.02, 0, -np.inf, 0],
            [two_theta.max(), np.inf, 3.0, 1, np.inf, np.inf]
        )
        
        try:
            popt, pcov = curve_fit(
                self._doublet_profile,
                two_theta,
                intensity,
                p0=p0,
                bounds=bounds,
                maxfev=self.max_iterations,
                ftol=self.tolerance,
                method='trf'
            )
            
            fitted = self._doublet_profile(two_theta, *popt)
            residuals = intensity - fitted
            chi_sq = np.sum(residuals ** 2)
            ss_tot = np.sum((intensity - np.mean(intensity)) ** 2)
            r_sq = 1 - (chi_sq / ss_tot) if ss_tot > 0 else 0
            
            center_ka2 = calculate_ka2_position(popt[0])
            
            return DoubletFitResult(
                center_ka1=popt[0],
                center_ka2=center_ka2,
                amplitude_ka1=popt[1],
                amplitude_ka2=popt[1] * self.intensity_ratio,
                fwhm=popt[2],
                eta=popt[3],
                r_squared=r_sq,
                success=True,
                message=f"Doublet fit converged (Δ2θ = {center_ka2 - popt[0]:.4f}°)",
                fitted_curve=fitted
            )
            
        except Exception as e:
            return DoubletFitResult(
                center_ka1=initial_center,
                center_ka2=calculate_ka2_position(initial_center),
                amplitude_ka1=peak_amp,
                amplitude_ka2=peak_amp * self.intensity_ratio,
                fwhm=initial_fwhm,
                eta=0.5,
                r_squared=0,
                success=False,
                message=str(e)
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
