"""
Peak Detection Module
Implements second derivative method for automatic peak finding.
"""

import numpy as np
from scipy.signal import savgol_filter, find_peaks as scipy_find_peaks
from typing import List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class PeakInfo:
    """Data class for peak information."""
    index: int
    two_theta: float
    intensity: float
    estimated_fwhm: float


class PeakDetector:
    """
    Automatic peak detection using second derivative method.
    
    The second derivative identifies inflection points where
    peaks have local minima (negative curvature).
    """
    
    def __init__(
        self,
        min_height: float = 100,
        min_distance: float = 0.5,
        window_size: int = 11,
        poly_order: int = 3
    ):
        """
        Initialize peak detector.
        
        Args:
            min_height: Minimum peak intensity threshold
            min_distance: Minimum distance between peaks (degrees)
            window_size: Savitzky-Golay window size
            poly_order: Savitzky-Golay polynomial order
        """
        self.min_height = min_height
        self.min_distance = min_distance
        self.window_size = window_size
        self.poly_order = poly_order
    
    def detect(
        self, 
        two_theta: np.ndarray, 
        intensity: np.ndarray
    ) -> List[PeakInfo]:
        """
        Detect peaks in XRD data.
        
        Args:
            two_theta: 2θ angle array
            intensity: Intensity array
            
        Returns:
            List of PeakInfo objects
        """
        # Calculate step size
        step = np.mean(np.diff(two_theta))
        min_distance_points = max(1, int(self.min_distance / step))
        
        # Compute second derivative
        deriv2 = savgol_filter(
            intensity, 
            self.window_size, 
            self.poly_order, 
            deriv=2
        )
        
        # Find peaks in original data
        peaks_idx, properties = scipy_find_peaks(
            intensity,
            height=self.min_height,
            distance=min_distance_points,
            prominence=self.min_height * 0.1
        )
        
        # Estimate FWHM for each peak
        peaks = []
        for idx in peaks_idx:
            fwhm = self._estimate_fwhm(two_theta, intensity, idx)
            peaks.append(PeakInfo(
                index=idx,
                two_theta=two_theta[idx],
                intensity=intensity[idx],
                estimated_fwhm=fwhm
            ))
        
        return peaks
    
    def _estimate_fwhm(
        self, 
        two_theta: np.ndarray, 
        intensity: np.ndarray, 
        peak_idx: int
    ) -> float:
        """
        Estimate FWHM for a single peak.
        
        Args:
            two_theta: 2θ array
            intensity: Intensity array
            peak_idx: Peak index
            
        Returns:
            Estimated FWHM in degrees
        """
        peak_height = intensity[peak_idx]
        half_max = peak_height / 2
        
        # Find left half-max point
        left_idx = peak_idx
        while left_idx > 0 and intensity[left_idx] > half_max:
            left_idx -= 1
        
        # Find right half-max point
        right_idx = peak_idx
        while right_idx < len(intensity) - 1 and intensity[right_idx] > half_max:
            right_idx += 1
        
        # Calculate FWHM
        fwhm = two_theta[right_idx] - two_theta[left_idx]
        return max(fwhm, 0.05)  # Minimum 0.05 degrees


def find_peaks(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    min_height: float = 100,
    min_distance: float = 0.5
) -> List[PeakInfo]:
    """
    Convenience function for peak detection.
    
    Args:
        two_theta: 2θ angle array
        intensity: Intensity array
        min_height: Minimum peak intensity
        min_distance: Minimum distance between peaks (degrees)
        
    Returns:
        List of PeakInfo objects
    """
    detector = PeakDetector(min_height=min_height, min_distance=min_distance)
    return detector.detect(two_theta, intensity)
