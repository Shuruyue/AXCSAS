"""
Pseudo-Voigt Function Module
Implements the Pseudo-Voigt profile function for XRD peak fitting.
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass


@dataclass
class PseudoVoigtParams:
    """Parameters for a Pseudo-Voigt peak."""
    center: float       # 2θ center position (degrees)
    amplitude: float    # Peak amplitude (intensity)
    fwhm: float        # Full Width at Half Maximum (degrees)
    eta: float         # Mixing parameter (0 = Gaussian, 1 = Lorentzian)
    
    def to_array(self) -> np.ndarray:
        """Convert to numpy array [center, amplitude, fwhm, eta]."""
        return np.array([self.center, self.amplitude, self.fwhm, self.eta])
    
    @classmethod
    def from_array(cls, arr: np.ndarray) -> 'PseudoVoigtParams':
        """Create from numpy array."""
        return cls(
            center=arr[0],
            amplitude=arr[1],
            fwhm=arr[2],
            eta=arr[3]
        )


class PseudoVoigt:
    """
    Pseudo-Voigt profile function.
    
    The Pseudo-Voigt is a linear combination of Gaussian and Lorentzian:
    
    I(2θ) = I₀ × [η × L(2θ) + (1-η) × G(2θ)]
    
    where:
    - L(2θ): Lorentzian component (size broadening)
    - G(2θ): Gaussian component (strain + instrumental)
    - η: Mixing parameter (0 < η < 1)
    """
    
    @staticmethod
    def gaussian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
        """
        Normalized Gaussian function.
        
        G(x) = exp(-4 ln(2) × ((x - center) / fwhm)²)
        """
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        return np.exp(-0.5 * ((x - center) / sigma) ** 2)
    
    @staticmethod
    def lorentzian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
        """
        Normalized Lorentzian function.
        
        L(x) = 1 / (1 + 4 × ((x - center) / fwhm)²)
        """
        gamma = fwhm / 2
        return 1 / (1 + ((x - center) / gamma) ** 2)
    
    @classmethod
    def profile(
        cls,
        x: np.ndarray,
        center: float,
        amplitude: float,
        fwhm: float,
        eta: float
    ) -> np.ndarray:
        """
        Calculate Pseudo-Voigt profile.
        
        Args:
            x: 2θ array
            center: Peak center position
            amplitude: Peak amplitude
            fwhm: Full Width at Half Maximum
            eta: Mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian)
            
        Returns:
            Intensity array
        """
        # Ensure eta is bounded
        eta = np.clip(eta, 0, 1)
        
        gaussian = cls.gaussian(x, center, fwhm)
        lorentzian = cls.lorentzian(x, center, fwhm)
        
        return amplitude * (eta * lorentzian + (1 - eta) * gaussian)
    
    @classmethod
    def multi_peak(
        cls,
        x: np.ndarray,
        params_list: list,
        background: float = 0
    ) -> np.ndarray:
        """
        Calculate sum of multiple Pseudo-Voigt peaks.
        
        Args:
            x: 2θ array
            params_list: List of PseudoVoigtParams or parameter arrays
            background: Constant background level
            
        Returns:
            Total intensity array
        """
        result = np.full_like(x, background, dtype=float)
        
        for params in params_list:
            if isinstance(params, PseudoVoigtParams):
                center, amplitude, fwhm, eta = (
                    params.center, params.amplitude, params.fwhm, params.eta
                )
            else:
                center, amplitude, fwhm, eta = params[:4]
            
            result += cls.profile(x, center, amplitude, fwhm, eta)
        
        return result


def pseudo_voigt_function(
    x: np.ndarray,
    center: float,
    amplitude: float,
    fwhm: float,
    eta: float
) -> np.ndarray:
    """
    Convenience function for Pseudo-Voigt calculation.
    
    Args:
        x: 2θ array
        center: Peak center (degrees)
        amplitude: Peak amplitude
        fwhm: Full Width at Half Maximum (degrees)
        eta: Mixing parameter (0-1)
        
    Returns:
        Intensity array
    """
    return PseudoVoigt.profile(x, center, amplitude, fwhm, eta)
