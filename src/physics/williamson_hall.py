"""
Williamson-Hall Analysis Module
Separates crystallite size and microstrain contributions to peak broadening.

Reference: Williamson, G. K., & Hall, W. H. (1953). 
           X-ray line broadening from filed aluminium and wolfram.
           Acta Metallurgica, 1(1), 22-31.
"""

import numpy as np
from typing import Tuple, Optional, List, NamedTuple
from dataclasses import dataclass


@dataclass
class WHResult:
    """Result of Williamson-Hall analysis."""
    crystallite_size_nm: float    # Crystallite size in nm
    microstrain: float            # Dimensionless microstrain ε
    r_squared: float              # Goodness of fit R²
    intercept: float              # y-intercept (Kλ/D)
    slope: float                  # Slope (4ε)
    is_reliable: bool             # Whether fit is reliable
    warning: Optional[str]        # Warning message if any
    
    def __repr__(self):
        status = "✓" if self.is_reliable else "⚠"
        return (f"WHResult({status} D={self.crystallite_size_nm:.1f} nm, "
                f"ε={self.microstrain:.2e}, R²={self.r_squared:.3f})")


class WilliamsonHallAnalyzer:
    """
    Williamson-Hall analysis for separating size and strain broadening.
    
    The W-H equation:
    
    β cos θ = (K λ / D) + 4 ε sin θ
    
    where:
    - β: Total peak broadening (FWHM in radians)
    - θ: Bragg angle
    - K: Scherrer constant
    - λ: X-ray wavelength
    - D: Crystallite size
    - ε: Microstrain
    
    By plotting β cos θ vs sin θ, we get:
    - Intercept = K λ / D → D = K λ / intercept
    - Slope = 4 ε → ε = slope / 4
    """
    
    # Default parameters
    DEFAULT_WAVELENGTH = 1.54056  # Cu Kα1 in Å
    DEFAULT_K = 0.89              # Spherical grains
    
    # Reliability thresholds
    MIN_R_SQUARED = 0.7
    MIN_PEAKS = 3
    
    def __init__(
        self,
        wavelength: float = DEFAULT_WAVELENGTH,
        k_factor: float = DEFAULT_K
    ):
        """
        Initialize Williamson-Hall analyzer.
        
        Args:
            wavelength: X-ray wavelength in Ångströms
            k_factor: Scherrer constant K
        """
        self.wavelength = wavelength
        self.k_factor = k_factor
    
    def analyze(
        self,
        two_theta: np.ndarray,
        fwhm: np.ndarray,
        fwhm_in_radians: bool = False
    ) -> WHResult:
        """
        Perform Williamson-Hall analysis.
        
        Args:
            two_theta: Array of 2θ peak positions (degrees)
            fwhm: Array of FWHM values (degrees or radians)
            fwhm_in_radians: If True, FWHM is in radians
            
        Returns:
            WHResult with crystallite size and strain
        """
        two_theta = np.asarray(two_theta)
        fwhm = np.asarray(fwhm)
        
        # Validate input
        if len(two_theta) != len(fwhm):
            raise ValueError("two_theta and fwhm must have same length")
        
        n_peaks = len(two_theta)
        warning = None
        is_reliable = True
        
        if n_peaks < self.MIN_PEAKS:
            warning = f"Only {n_peaks} peaks, need at least {self.MIN_PEAKS} for reliable W-H analysis"
            is_reliable = False
        
        # Convert FWHM to radians if needed
        if fwhm_in_radians:
            beta = fwhm
        else:
            beta = np.radians(fwhm)
        
        # Calculate θ in radians
        theta = np.radians(two_theta / 2)
        
        # W-H plot coordinates
        # x-axis: sin θ
        # y-axis: β cos θ
        x = np.sin(theta)
        y = beta * np.cos(theta)
        
        # Linear regression: y = slope * x + intercept
        slope, intercept = self._linear_fit(x, y)
        
        # Calculate R²
        r_squared = self._calculate_r_squared(x, y, slope, intercept)
        
        # Check R² reliability
        if r_squared < self.MIN_R_SQUARED:
            if warning:
                warning += f"; Low R² ({r_squared:.3f})"
            else:
                warning = f"Low R² ({r_squared:.3f}), fit may be unreliable"
            is_reliable = False
        
        # Extract physical parameters
        # intercept = K λ / D → D = K λ / intercept
        # ε = slope / 4
        
        if intercept > 0:
            size_angstrom = (self.k_factor * self.wavelength) / intercept
            size_nm = size_angstrom / 10.0
        else:
            size_nm = float('inf')
            if warning:
                warning += "; Negative intercept"
            else:
                warning = "Negative intercept, size estimation unreliable"
            is_reliable = False
        
        microstrain = slope / 4.0
        
        # Check for negative strain (physically unreasonable for uniform strain)
        if microstrain < 0:
            if warning:
                warning += "; Negative strain"
            else:
                warning = "Negative strain detected, may indicate non-uniform strain"
        
        return WHResult(
            crystallite_size_nm=size_nm,
            microstrain=abs(microstrain),
            r_squared=r_squared,
            intercept=intercept,
            slope=slope,
            is_reliable=is_reliable,
            warning=warning
        )
    
    def analyze_with_correction(
        self,
        two_theta: np.ndarray,
        fwhm_observed: np.ndarray,
        fwhm_instrumental: np.ndarray
    ) -> WHResult:
        """
        Perform W-H analysis with instrumental broadening correction.
        
        Uses geometric approximation for Pseudo-Voigt:
        β_sample = β_obs - β²_inst / β_obs
        
        Args:
            two_theta: Peak positions (degrees)
            fwhm_observed: Observed FWHM (degrees)
            fwhm_instrumental: Instrumental FWHM (degrees)
            
        Returns:
            WHResult
        """
        fwhm_obs = np.asarray(fwhm_observed)
        fwhm_inst = np.asarray(fwhm_instrumental)
        
        # Geometric correction
        fwhm_corrected = fwhm_obs - (fwhm_inst**2 / fwhm_obs)
        
        # Ensure positive values
        fwhm_corrected = np.maximum(fwhm_corrected, 0.001)
        
        return self.analyze(two_theta, fwhm_corrected)
    
    def get_plot_data(
        self,
        two_theta: np.ndarray,
        fwhm: np.ndarray,
        fwhm_in_radians: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get data for W-H plot.
        
        Args:
            two_theta: Peak positions (degrees)
            fwhm: FWHM values
            fwhm_in_radians: If True, FWHM is in radians
            
        Returns:
            Tuple of (x_data, y_data, x_fit, y_fit)
        """
        two_theta = np.asarray(two_theta)
        fwhm = np.asarray(fwhm)
        
        if not fwhm_in_radians:
            beta = np.radians(fwhm)
        else:
            beta = fwhm
        
        theta = np.radians(two_theta / 2)
        
        x_data = np.sin(theta)
        y_data = beta * np.cos(theta)
        
        # Fit line
        slope, intercept = self._linear_fit(x_data, y_data)
        
        x_fit = np.array([0, x_data.max() * 1.1])
        y_fit = slope * x_fit + intercept
        
        return x_data, y_data, x_fit, y_fit
    
    @staticmethod
    def _linear_fit(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
        """Simple linear least squares fit."""
        # y = slope * x + intercept
        n = len(x)
        sum_x = np.sum(x)
        sum_y = np.sum(y)
        sum_xy = np.sum(x * y)
        sum_xx = np.sum(x * x)
        
        denom = n * sum_xx - sum_x**2
        
        if abs(denom) < 1e-10:
            return (0.0, np.mean(y))
        
        slope = (n * sum_xy - sum_x * sum_y) / denom
        intercept = (sum_y - slope * sum_x) / n
        
        return slope, intercept
    
    @staticmethod
    def _calculate_r_squared(
        x: np.ndarray,
        y: np.ndarray,
        slope: float,
        intercept: float
    ) -> float:
        """Calculate coefficient of determination R²."""
        y_mean = np.mean(y)
        y_pred = slope * x + intercept
        
        ss_tot = np.sum((y - y_mean)**2)
        ss_res = np.sum((y - y_pred)**2)
        
        if ss_tot < 1e-10:
            return 0.0
        
        return 1 - (ss_res / ss_tot)


def williamson_hall_analysis(
    two_theta: np.ndarray,
    fwhm: np.ndarray,
    wavelength: float = 1.54056,
    k_factor: float = 0.89
) -> Tuple[float, float, float]:
    """
    Convenience function for Williamson-Hall analysis.
    
    Args:
        two_theta: Peak positions (degrees)
        fwhm: FWHM values (degrees)
        wavelength: X-ray wavelength (Å)
        k_factor: Scherrer constant
        
    Returns:
        Tuple of (crystallite_size_nm, microstrain, r_squared)
    """
    analyzer = WilliamsonHallAnalyzer(wavelength, k_factor)
    result = analyzer.analyze(two_theta, fwhm)
    return (result.crystallite_size_nm, result.microstrain, result.r_squared)
