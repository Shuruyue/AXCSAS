"""
Smoothing Module
Implements Savitzky-Golay filter for XRD data smoothing.
"""

import numpy as np
from scipy.signal import savgol_filter
from typing import Optional


class SavitzkyGolayFilter:
    """
    Savitzky-Golay filter for XRD data smoothing.
    
    This filter preserves peak characteristics (position, height, width)
    while removing high-frequency noise.
    
    Recommended parameters:
    - window_size: 5-15 points (must be odd)
    - poly_order: 2-3
    """
    
    def __init__(self, window_size: int = 11, poly_order: int = 3):
        """
        Initialize the Savitzky-Golay filter.
        
        Args:
            window_size: Filter window size (must be odd, default: 11)
            poly_order: Polynomial order for fitting (default: 3)
        """
        if window_size % 2 == 0:
            raise ValueError("window_size must be odd")
        if poly_order >= window_size:
            raise ValueError("poly_order must be less than window_size")
        
        self.window_size = window_size
        self.poly_order = poly_order
    
    def apply(self, intensity: np.ndarray) -> np.ndarray:
        """
        Apply Savitzky-Golay smoothing to intensity data.
        
        Args:
            intensity: Raw intensity array
            
        Returns:
            Smoothed intensity array
            
        Raises:
            ValueError: If window_size >= len(intensity)
        """
        if self.window_size >= len(intensity):
            raise ValueError(
                f"Window size ({self.window_size}) must be less than "
                f"data length ({len(intensity)})"
            )
        return savgol_filter(
            intensity, 
            window_length=self.window_size,
            polyorder=self.poly_order
        )
    
    def apply_derivative(
        self, 
        intensity: np.ndarray, 
        deriv: int = 1
    ) -> np.ndarray:
        """
        Apply Savitzky-Golay filter and compute derivative.
        
        Useful for peak detection (2nd derivative) and 
        background identification.
        
        Args:
            intensity: Raw intensity array
            deriv: Derivative order (1 or 2)
            
        Returns:
            Derivative of smoothed intensity
        """
        return savgol_filter(
            intensity,
            window_length=self.window_size,
            polyorder=self.poly_order,
            deriv=deriv
        )


def smooth_xrd_data(
    intensity: np.ndarray,
    window_size: int = 11,
    poly_order: int = 3
) -> np.ndarray:
    """
    Convenience function for XRD data smoothing.
    
    Args:
        intensity: Raw intensity array
        window_size: Filter window size (must be odd)
        poly_order: Polynomial order
        
    Returns:
        Smoothed intensity array
    """
    sg_filter = SavitzkyGolayFilter(window_size, poly_order)
    return sg_filter.apply(intensity)
