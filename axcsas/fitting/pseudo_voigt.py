"""
Pseudo-Voigt Function Module 偽Voigt函數模組
=============================================
Implements the Pseudo-Voigt profile function for XRD peak fitting.
實現用於 XRD 峰擬合的偽Voigt剖面函數。

Also includes True Voigt profile using scipy.special.voigt_profile.
另包含使用 scipy.special.voigt_profile 的真Voigt剖面。
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass
from scipy.special import voigt_profile


@dataclass
class VoigtParams:
    """Parameters for a true Voigt peak."""
    center: float       # 2θ center position (degrees)
    amplitude: float    # Peak amplitude (intensity)
    sigma: float        # Gaussian width parameter
    gamma: float        # Lorentzian width parameter (half-width)
    
    @property
    def fwhm_gaussian(self) -> float:
        """FWHM of Gaussian component."""
        return 2.0 * self.sigma * np.sqrt(2.0 * np.log(2.0))
    
    @property
    def fwhm_lorentzian(self) -> float:
        """FWHM of Lorentzian component."""
        return 2.0 * self.gamma
    
    @property
    def fwhm_total(self) -> float:
        """Approximate total FWHM using Olivero-Longbothum approximation."""
        fG = self.fwhm_gaussian
        fL = self.fwhm_lorentzian
        return 0.5346 * fL + np.sqrt(0.2166 * fL**2 + fG**2)
    
    def to_array(self) -> np.ndarray:
        """Convert to numpy array [center, amplitude, sigma, gamma]."""
        return np.array([self.center, self.amplitude, self.sigma, self.gamma])
    
    @classmethod
    def from_array(cls, arr: np.ndarray) -> 'VoigtParams':
        """Create from numpy array."""
        return cls(center=arr[0], amplitude=arr[1], sigma=arr[2], gamma=arr[3])


class TrueVoigt:
    """
    True Voigt profile function (convolution of Gaussian and Lorentzian).
    
    Mathematical Definition (Faddeeva Function):
    --------------------------------------------
    V(x; σ, γ) = Re[w(z)] / (σ√(2π))
    
    where w(z) is the Faddeeva function (scaled complex error function):
        w(z) = exp(-z²) * erfc(-iz)
        
    and z is the complex argument:
        z = (x + iγ) / (σ√2)
        
    Physical Interpretation in XRD:
    -------------------------------
    - x (Real part): Deviation from Bragg angle (2θ - 2θ₀).
      Analogous to "Frequency Drift" in spectroscopy.
    - γ (Imaginary part): Lorentzian HWHM (Size Broadening + Lifetime).
      Analogous to "Damping Coefficient" in spectroscopy.
    - σ (Normalization): Gaussian Width (Strain + Instrument).
    
    Note on Instrument Parameters:
    ------------------------------
    This function fits the *Total* profile (Sample ⊗ Instrument).
    You do NOT need instrument parameters to perform this fit. The optimizer
    will find the effective σ_total and γ_total from the raw data.
    Instrumental correction is performed *after* fitting during the
    analysis stage (e.g., in Williamson-Hall or Scherrer methods).
    
    References:
    -----------
    1. Armstrong, B. H. (1967). 
       "Spectrum Line Profiles: The Voigt Function". 
       J. Quant. Spectrosc. Radiat. Transfer, 7, 61-88.
       
    2. Poppe, G. P. M., & Wijers, C. M. J. (1990).
       "More efficient computation of the complex error function".
       ACM Trans. Math. Softw., 16, 38-46.
       (Standard algorithm used by scipy.special.wofz)
    """
    
    @staticmethod
    def profile(
        x: np.ndarray,
        center: float,
        amplitude: float,
        sigma: float,
        gamma: float
    ) -> np.ndarray:
        """
        Calculate true Voigt profile.
        
        Args:
            x: 2θ array
            center: Peak center position
            amplitude: Peak amplitude
            sigma: Gaussian width parameter (σ)
            gamma: Lorentzian half-width parameter (γ)
            
        Returns:
            Intensity array
        """
        # Shift x to be centered at 0
        x_shifted = x - center
        
        # Use scipy's voigt_profile (normalized to unit area)
        # voigt_profile(x, sigma, gamma) computes the Voigt function
        v = voigt_profile(x_shifted, sigma, gamma)
        
        # Normalize and scale by amplitude
        # The voigt_profile is normalized, so we scale to match peak amplitude
        v_max = voigt_profile(0, sigma, gamma)
        if v_max > 0:
            return amplitude * (v / v_max)
        else:
            return np.zeros_like(x)
    
    @staticmethod
    def fwhm_from_params(sigma: float, gamma: float) -> float:
        """
        Calculate total FWHM from sigma and gamma.
        Uses Olivero-Longbothum approximation.
        """
        fG = 2.0 * sigma * np.sqrt(2.0 * np.log(2.0))
        fL = 2.0 * gamma
        return 0.5346 * fL + np.sqrt(0.2166 * fL**2 + fG**2)
    
    @staticmethod
    def params_from_fwhm(fwhm: float, eta: float = 0.5) -> Tuple[float, float]:
        """
        Estimate sigma and gamma from total FWHM and mixing ratio.
        
        Args:
            fwhm: Total FWHM
            eta: Lorentzian fraction (0 = Gaussian, 1 = Lorentzian)
            
        Returns:
            (sigma, gamma) tuple
        """
        # Approximate split between Gaussian and Lorentzian FWHM
        fL = eta * fwhm
        fG = (1 - eta) * fwhm
        
        sigma = fG / (2.0 * np.sqrt(2.0 * np.log(2.0))) if fG > 0 else 0.01
        gamma = fL / 2.0 if fL > 0 else 0.01
        
        return sigma, gamma

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
