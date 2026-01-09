"""
Scherrer Equation Module
Implements crystallite size calculation from XRD peak broadening.

Reference: Langford, J. I., & Wilson, A. J. C. (1978). 
           Scherrer after sixty years. J. Appl. Cryst., 11, 102-113.
"""

import numpy as np
from typing import Optional, Tuple, NamedTuple
from dataclasses import dataclass
from enum import Enum


class GrainShape(Enum):
    """Grain shape for Scherrer constant selection."""
    SPHERICAL = "spherical"
    CUBIC = "cubic"
    CUSTOM = "custom"


# Scherrer constants from Langford & Wilson (1978)
SCHERRER_CONSTANTS = {
    GrainShape.SPHERICAL: 0.89,
    GrainShape.CUBIC: 0.94,
}


@dataclass
class ScherrerResult:
    """Result of Scherrer crystallite size calculation."""
    size_nm: float           # Crystallite size in nanometers
    size_angstrom: float     # Crystallite size in Ångströms
    peak_position: float     # 2θ position (degrees)
    fwhm_corrected: float    # Corrected FWHM (degrees)
    is_reliable: bool        # Whether result is within reliable range
    warning: Optional[str]   # Warning message if any
    
    def __repr__(self):
        status = "✓" if self.is_reliable else "⚠"
        return f"ScherrerResult({status} D={self.size_nm:.1f} nm @ 2θ={self.peak_position:.2f}°)"


class ScherrerCalculator:
    """
    Scherrer equation for crystallite size calculation.
    
    D = K × λ / (β × cos θ)
    
    where:
    - D: Crystallite size
    - K: Scherrer constant (shape factor)
    - λ: X-ray wavelength
    - β: Peak broadening (FWHM in radians, corrected for instrumental)
    - θ: Bragg angle (half of 2θ)
    
    Attributes:
        wavelength: X-ray wavelength in Ångströms
        k_factor: Scherrer constant K
        grain_shape: Assumed grain shape
    """
    
    # Default Cu Kα1 wavelength
    DEFAULT_WAVELENGTH = 1.54056  # Å
    
    # Reliable size range (nm)
    MIN_RELIABLE_SIZE = 2.0
    MAX_RELIABLE_SIZE = 200.0
    
    def __init__(
        self,
        wavelength: float = DEFAULT_WAVELENGTH,
        grain_shape: GrainShape = GrainShape.SPHERICAL,
        k_factor: Optional[float] = None
    ):
        """
        Initialize Scherrer calculator.
        
        Args:
            wavelength: X-ray wavelength in Ångströms (default: Cu Kα1)
            grain_shape: Grain shape assumption for K selection
            k_factor: Custom K value (overrides grain_shape if provided)
        """
        self.wavelength = wavelength
        self.grain_shape = grain_shape
        
        if k_factor is not None:
            self.k_factor = k_factor
            self.grain_shape = GrainShape.CUSTOM
        else:
            self.k_factor = SCHERRER_CONSTANTS.get(grain_shape, 0.89)
    
    def calculate(
        self,
        two_theta: float,
        fwhm: float,
        fwhm_is_radians: bool = False
    ) -> ScherrerResult:
        """
        Calculate crystallite size using Scherrer equation.
        
        Args:
            two_theta: Peak position in degrees (2θ)
            fwhm: Full Width at Half Maximum (degrees or radians)
            fwhm_is_radians: If True, fwhm is in radians; otherwise degrees
            
        Returns:
            ScherrerResult with calculated size and metadata
        """
        # Convert FWHM to radians if needed
        if fwhm_is_radians:
            beta_rad = fwhm
            fwhm_deg = np.degrees(fwhm)
        else:
            beta_rad = np.radians(fwhm)
            fwhm_deg = fwhm
        
        # Validate inputs
        if beta_rad <= 0:
            raise ValueError(f"FWHM must be positive, got {fwhm}")
        if two_theta <= 0 or two_theta >= 180:
            raise ValueError(f"2θ must be between 0 and 180°, got {two_theta}")
        
        # Calculate Bragg angle θ in radians
        theta_rad = np.radians(two_theta / 2)
        cos_theta = np.cos(theta_rad)
        
        # Scherrer equation: D = K × λ / (β × cos θ)
        # Result in Ångströms (since λ is in Å)
        size_angstrom = (self.k_factor * self.wavelength) / (beta_rad * cos_theta)
        
        # Convert to nanometers
        size_nm = size_angstrom / 10.0
        
        # Check reliability
        warning = None
        is_reliable = True
        
        if size_nm < self.MIN_RELIABLE_SIZE:
            warning = f"Size {size_nm:.1f} nm below precision limit ({self.MIN_RELIABLE_SIZE} nm)"
            is_reliable = False
        elif size_nm > self.MAX_RELIABLE_SIZE:
            warning = f"Size {size_nm:.1f} nm exceeds detection limit ({self.MAX_RELIABLE_SIZE} nm)"
            is_reliable = False
        
        return ScherrerResult(
            size_nm=size_nm,
            size_angstrom=size_angstrom,
            peak_position=two_theta,
            fwhm_corrected=fwhm_deg,
            is_reliable=is_reliable,
            warning=warning
        )
    
    def calculate_with_correction(
        self,
        two_theta: float,
        fwhm_observed: float,
        fwhm_instrumental: float
    ) -> ScherrerResult:
        """
        Calculate crystallite size with instrumental broadening correction.
        
        Uses geometric approximation for Pseudo-Voigt profiles:
        β_sample = β_obs - β²_inst / β_obs
        
        Args:
            two_theta: Peak position (degrees)
            fwhm_observed: Observed FWHM (degrees)
            fwhm_instrumental: Instrumental FWHM at this angle (degrees)
            
        Returns:
            ScherrerResult with calculated size
        """
        # Geometric correction for Pseudo-Voigt
        if fwhm_observed <= fwhm_instrumental:
            # Can't correct - instrumental dominates
            return ScherrerResult(
                size_nm=float('inf'),
                size_angstrom=float('inf'),
                peak_position=two_theta,
                fwhm_corrected=0.0,
                is_reliable=False,
                warning="Observed FWHM smaller than instrumental broadening"
            )
        
        fwhm_corrected = fwhm_observed - (fwhm_instrumental**2 / fwhm_observed)
        
        # Ensure positive FWHM
        fwhm_corrected = max(fwhm_corrected, 0.001)
        
        return self.calculate(two_theta, fwhm_corrected)
    
    def batch_calculate(
        self,
        peaks: list,
        fwhm_is_radians: bool = False
    ) -> list:
        """
        Calculate crystallite sizes for multiple peaks.
        
        Args:
            peaks: List of (two_theta, fwhm) tuples
            fwhm_is_radians: If True, FWHM values are in radians
            
        Returns:
            List of ScherrerResult objects
        """
        results = []
        for two_theta, fwhm in peaks:
            result = self.calculate(two_theta, fwhm, fwhm_is_radians)
            results.append(result)
        return results
    
    def average_size(self, results: list) -> Tuple[float, float]:
        """
        Calculate average crystallite size from multiple peaks.
        
        Only includes reliable results.
        
        Args:
            results: List of ScherrerResult objects
            
        Returns:
            Tuple of (average_size_nm, std_dev_nm)
        """
        reliable = [r.size_nm for r in results if r.is_reliable]
        
        if not reliable:
            return (float('nan'), float('nan'))
        
        avg = np.mean(reliable)
        std = np.std(reliable) if len(reliable) > 1 else 0.0
        
        return (avg, std)


def calculate_crystallite_size(
    two_theta: float,
    fwhm: float,
    wavelength: float = 1.54056,
    k_factor: float = 0.89
) -> float:
    """
    Convenience function for quick crystallite size calculation.
    
    Args:
        two_theta: Peak position (degrees)
        fwhm: FWHM in degrees
        wavelength: X-ray wavelength (Å), default Cu Kα1
        k_factor: Scherrer constant, default 0.89 (spherical)
        
    Returns:
        Crystallite size in nanometers
    """
    calc = ScherrerCalculator(wavelength=wavelength, k_factor=k_factor)
    result = calc.calculate(two_theta, fwhm)
    return result.size_nm
