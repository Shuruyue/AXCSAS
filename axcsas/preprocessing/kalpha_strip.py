"""
Kα2 Stripping Module Kα2 剥離模組
==================================
Implements Rachinger Correction for Cu Kα2 removal.
實現用於移除 Cu Kα2 的 Rachinger 校正。

Reference 出處:
    Rachinger, W. A. (1948). J. Sci. Instr., 25(7), 254-255.
"""

import numpy as np
from typing import Tuple

from axcsas.core.constants import CU_KA1, CU_KA2, KA2_KA1_RATIO


class KalphaStripper:
    """
    Kα2 stripping using Rachinger Correction.
    使用 Rachinger 校正的 Kα2 剥離。
    
    For Cu Kα wavelengths (Bearden 1967, Rev. Mod. Phys. 39, 78):
    - Kα1: 1.540562 Å
    - Kα2: 1.544390 Å
    - Intensity ratio: Kα2/Kα1 ≈ 0.5 (Burger-Dorgelo rule)
    """

    def __init__(
        self,
        ka1_lambda: float = CU_KA1,
        ka2_lambda: float = CU_KA2,
        ka_ratio: float = KA2_KA1_RATIO
    ):
        """
        Initialize Kα2 stripper.
        
        Args:
            wavelength_ka2: Kα2 wavelength in Å (default: Cu Kα2)
            ka2_ratio: I(Kα2) / I(Kα1) ratio (default: 0.5)
        """
        self.wavelength_ka1 = ka1_lambda
        self.wavelength_ka2 = ka2_lambda
        self.ka2_ratio = ka_ratio
    
    def strip(
        self, 
        two_theta: np.ndarray, 
        intensity: np.ndarray
    ) -> np.ndarray:
        """
        Remove Kα2 contribution using Rachinger correction.
        
        Algorithm:
        1. For each point, calculate the 2θ position where Kα2 
           would produce its peak (shifted from Kα1)
        2. Recursively subtract Kα2 contribution
        
        Args:
            two_theta: 2θ angle array (degrees)
            intensity: Intensity array
            
        Returns:
            Kα1-only intensity array
        """
        # Calculate angular shift between Kα2 and Kα1
        # Using Bragg's law: sin(θ₂)/sin(θ₁) = λ₂/λ₁
        
        corrected = intensity.copy()
        step = np.mean(np.diff(two_theta))
        
        for i in range(len(corrected)):
            theta_rad = np.radians(two_theta[i] / 2)
            sin_theta = np.sin(theta_rad)
            
            # Calculate the 2θ position for Kα2
            sin_theta_ka2 = sin_theta * (self.wavelength_ka2 / self.wavelength_ka1)
            
            if sin_theta_ka2 > 1:
                continue
            
            theta_ka2_rad = np.arcsin(sin_theta_ka2)
            two_theta_ka2 = 2 * np.degrees(theta_ka2_rad)
            
            # Find index offset
            delta_theta = two_theta_ka2 - two_theta[i]
            idx_offset = int(round(delta_theta / step))
            
            # Subtract Kα2 contribution
            if 0 <= i - idx_offset < len(corrected):
                ka2_contribution = self.ka2_ratio * corrected[i - idx_offset]
                corrected[i] = max(0, corrected[i] - ka2_contribution)
        
        return corrected
    
    def estimate_angular_shift(self, two_theta: float) -> float:
        """
        Estimate the angular shift between Kα1 and Kα2 peaks.
        
        Args:
            two_theta: 2θ angle in degrees
            
        Returns:
            Angular shift in degrees
        """
        # Calculate 2θ shift: Δ(2θ) = 2 * arcsin(λ2/λ1 * sin(θ1)) - 2θ1
        theta_rad = np.radians(two_theta / 2.0)
        sin_theta1 = np.sin(theta_rad)
        
        # Check domain for arcsin
        term = (self.wavelength_ka2 / self.wavelength_ka1) * sin_theta1       
        if term > 1:
            return float('nan')
        
        theta_ka2_rad = np.arcsin(term)
        return 2 * np.degrees(theta_ka2_rad) - two_theta


def strip_kalpha2(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    ka2_ratio: float = 0.5
) -> np.ndarray:
    """
    Convenience function for Kα2 stripping.
    
    Args:
        two_theta: 2θ angle array (degrees)
        intensity: Intensity array
        ka2_ratio: I(Kα2) / I(Kα1) ratio
        
    Returns:
        Kα1-only intensity array
    """
    stripper = KalphaStripper(ka2_ratio=ka2_ratio)
    return stripper.strip(two_theta, intensity)
