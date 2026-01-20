"""
Unit Conversion Module 單位轉換模組
===================================
Provides unit conversion utilities for XRD analysis.
提供 XRD 分析所需的單位轉換工具。

Standard definitions 標準定義:
- 1 Å (Ångström) = 10⁻¹⁰ m = 0.1 nm
- 2θ = 2 × θ (Bragg angle)

All functions are pure mathematical transformations with no physical constants.
所有函數皆為純數學轉換，不含物理常數。
"""

import numpy as np
from typing import Union

from axcsas.core.constants import CU_KA1


# =============================================================================
# Angle Conversions 角度轉換
# =============================================================================

def deg_to_rad(degrees: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert degrees to radians.
    
    Args:
        degrees: Angle in degrees
        
    Returns:
        Angle in radians
    """
    return np.radians(degrees)


def rad_to_deg(radians: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert radians to degrees.
    
    Args:
        radians: Angle in radians
        
    Returns:
        Angle in degrees
    """
    return np.degrees(radians)


def two_theta_to_theta(two_theta: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert 2θ to θ (Bragg angle).
    
    Args:
        two_theta: Diffraction angle 2θ (degrees)
        
    Returns:
        Bragg angle θ (degrees)
    """
    return two_theta / 2.0


def theta_to_two_theta(theta: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert θ (Bragg angle) to 2θ.
    
    Args:
        theta: Bragg angle θ (degrees)
        
    Returns:
        Diffraction angle 2θ (degrees)
    """
    return theta * 2.0


# =============================================================================
# Length Conversions
# =============================================================================

def angstrom_to_nm(angstrom: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert Ångströms to nanometers.
    
    Args:
        angstrom: Length in Ångströms
        
    Returns:
        Length in nanometers
    """
    return angstrom / 10.0


def nm_to_angstrom(nm: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert nanometers to Ångströms.
    
    Args:
        nm: Length in nanometers
        
    Returns:
        Length in Ångströms
    """
    return nm * 10.0


def nm_to_meter(nm: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert nanometers to meters.
    
    Args:
        nm: Length in nanometers
        
    Returns:
        Length in meters
    """
    return nm * 1e-9


def meter_to_nm(meter: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert meters to nanometers.
    
    Args:
        meter: Length in meters
        
    Returns:
        Length in nanometers
    """
    return meter * 1e9


# =============================================================================
# FWHM/Beta Conversions
# =============================================================================

def fwhm_deg_to_rad(fwhm_deg: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert FWHM from degrees to radians.
    
    This is specifically for peak widths (FWHM/β).
    
    Args:
        fwhm_deg: FWHM in degrees
        
    Returns:
        FWHM in radians
    """
    return np.radians(fwhm_deg)


def fwhm_rad_to_deg(fwhm_rad: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert FWHM from radians to degrees.
    
    Args:
        fwhm_rad: FWHM in radians
        
    Returns:
        FWHM in degrees
    """
    return np.degrees(fwhm_rad)


# =============================================================================
# Crystallographic Conversions
# =============================================================================

def d_spacing_to_two_theta(
    d_spacing: float,
    wavelength: float = CU_KA1  # Default: Cu Kα1 from constants
) -> float:
    """
    Calculate 2θ from d-spacing using Bragg's law.
    
    Bragg's Law: nλ = 2d sin θ
    For n=1: sin θ = λ / (2d)
    
    Args:
        d_spacing: Interplanar spacing in Ångströms
        wavelength: X-ray wavelength in Ångströms
        
    Returns:
        2θ angle in degrees
    """
    if d_spacing <= 0:
        raise ValueError("d-spacing must be positive")
    
    sin_theta = wavelength / (2 * d_spacing)
    
    if sin_theta > 1:
        raise ValueError(f"d-spacing {d_spacing} Å too small for wavelength {wavelength} Å")
    
    theta_rad = np.arcsin(sin_theta)
    two_theta = 2 * np.degrees(theta_rad)
    
    return two_theta


def two_theta_to_d_spacing(
    two_theta: float,
    wavelength: float = CU_KA1  # Default: Cu Kα1 from constants
) -> float:
    """
    Calculate d-spacing from 2θ using Bragg's law.
    
    d = λ / (2 sin θ)
    
    Args:
        two_theta: 2θ angle in degrees
        wavelength: X-ray wavelength in Ångströms
        
    Returns:
        d-spacing in Ångströms
    """
    theta_rad = np.radians(two_theta / 2)
    sin_theta = np.sin(theta_rad)
    
    if sin_theta <= 0:
        raise ValueError(f"Invalid 2θ angle: {two_theta}")
    
    d_spacing = wavelength / (2 * sin_theta)
    
    return d_spacing
