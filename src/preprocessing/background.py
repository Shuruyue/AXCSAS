"""
Background Subtraction Module
Implements Sonneveld-Visser and Chebyshev polynomial methods.
"""

import numpy as np
from scipy.signal import savgol_filter
from typing import Tuple, Optional


class BackgroundSubtractor:
    """
    Background subtraction for XRD data.
    
    Two methods available:
    1. Sonneveld-Visser: Second derivative method
    2. Chebyshev: Polynomial fitting to peak-free regions
    """
    
    def __init__(self, method: str = "chebyshev", **kwargs):
        """
        Initialize background subtractor.
        
        Args:
            method: "chebyshev" or "sonneveld_visser"
            **kwargs: Method-specific parameters
                - poly_degree: For Chebyshev method (default: 5)
                - iterations: For Sonneveld-Visser (default: 100)
        """
        self.method = method.lower()
        self.poly_degree = kwargs.get('poly_degree', 5)
        self.iterations = kwargs.get('iterations', 100)
        
        if self.method not in ['chebyshev', 'sonneveld_visser']:
            raise ValueError(f"Unknown method: {method}")
    
    def subtract(
        self, 
        two_theta: np.ndarray, 
        intensity: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Subtract background from XRD data.
        
        Args:
            two_theta: 2θ angle array
            intensity: Intensity array
            
        Returns:
            Tuple of (corrected_intensity, background)
        """
        if self.method == 'chebyshev':
            return self._chebyshev_background(two_theta, intensity)
        else:
            return self._sonneveld_visser_background(two_theta, intensity)
    
    def _chebyshev_background(
        self, 
        two_theta: np.ndarray, 
        intensity: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Chebyshev polynomial background fitting.
        
        Iteratively fits polynomial, excluding points above the fit.
        """
        # Normalize x for numerical stability
        x_norm = 2 * (two_theta - two_theta.min()) / (two_theta.max() - two_theta.min()) - 1
        
        # Initial fit
        mask = np.ones(len(intensity), dtype=bool)
        
        for _ in range(self.iterations):
            # Fit Chebyshev polynomial to masked data
            coeffs = np.polynomial.chebyshev.chebfit(
                x_norm[mask], 
                intensity[mask], 
                self.poly_degree
            )
            background = np.polynomial.chebyshev.chebval(x_norm, coeffs)
            
            # Update mask: keep points below or near background
            new_mask = intensity <= background * 1.05
            
            if np.array_equal(mask, new_mask):
                break
            mask = new_mask
        
        corrected = intensity - background
        corrected[corrected < 0] = 0  # Ensure non-negative
        
        return corrected, background
    
    def _sonneveld_visser_background(
        self, 
        two_theta: np.ndarray, 
        intensity: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sonneveld-Visser second derivative method.
        
        Uses second derivative to identify background regions.
        """
        # Compute smoothed second derivative
        window = min(21, len(intensity) // 10)
        if window % 2 == 0:
            window += 1
        
        deriv2 = savgol_filter(intensity, window, 3, deriv=2)
        
        # Background points have near-zero second derivative
        threshold = np.std(deriv2) * 0.1
        bg_mask = np.abs(deriv2) < threshold
        
        if np.sum(bg_mask) < 10:
            # Fallback: use minimum envelope
            bg_mask = intensity < np.percentile(intensity, 20)
        
        # Interpolate background from identified points
        bg_indices = np.where(bg_mask)[0]
        bg_values = intensity[bg_mask]
        
        background = np.interp(
            np.arange(len(intensity)),
            bg_indices,
            bg_values
        )
        
        corrected = intensity - background
        corrected[corrected < 0] = 0
        
        return corrected, background


def subtract_background(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    method: str = "chebyshev",
    **kwargs
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convenience function for background subtraction.
    
    Args:
        two_theta: 2θ angle array
        intensity: Intensity array
        method: "chebyshev" or "sonneveld_visser"
        **kwargs: Method-specific parameters
        
    Returns:
        Tuple of (corrected_intensity, background)
    """
    subtractor = BackgroundSubtractor(method, **kwargs)
    return subtractor.subtract(two_theta, intensity)
