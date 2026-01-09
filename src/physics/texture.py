"""
Texture Analysis Module
Implements Harris Texture Coefficient for preferred orientation analysis.

Reference: Harris, G. B. (1952). Quantitative measurement of preferred 
           orientation in rolled uranium bars. The London, Edinburgh, 
           and Dublin Philosophical Magazine and Journal of Science, 
           43(336), 113-123.
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


# JCPDS standard data for Cu (PDF 04-0836)
# Format: (hkl, 2θ position, relative intensity %)
CU_JCPDS_STANDARD = {
    (1, 1, 1): {"two_theta": 43.298, "intensity": 100.0},
    (2, 0, 0): {"two_theta": 50.434, "intensity": 46.0},
    (2, 2, 0): {"two_theta": 74.130, "intensity": 20.0},
    (3, 1, 1): {"two_theta": 89.931, "intensity": 17.0},
    (2, 2, 2): {"two_theta": 95.139, "intensity": 5.0},
}


@dataclass
class TextureResult:
    """Result for a single reflection's texture coefficient."""
    hkl: Tuple[int, int, int]  # Miller indices
    tc_value: float             # Texture Coefficient value
    intensity_observed: float   # Observed intensity
    intensity_standard: float   # Standard intensity
    is_preferred: bool          # Whether this is a preferred orientation
    
    def __repr__(self):
        hkl_str = f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]})"
        status = "★" if self.is_preferred else ""
        return f"TC{hkl_str}={self.tc_value:.2f}{status}"


@dataclass 
class TextureAnalysisResult:
    """Complete texture analysis result."""
    tc_values: Dict[Tuple[int, int, int], float]  # TC for each (hkl)
    preferred_orientation: Optional[Tuple[int, int, int]]  # Highest TC
    degree_of_texture: float  # Degree of preferred orientation
    is_random: bool           # Whether texture is random (TC ≈ 1 for all)
    details: List[TextureResult]  # Detailed results
    
    def __repr__(self):
        if self.is_random:
            return "TextureAnalysisResult(Random texture)"
        else:
            hkl = self.preferred_orientation
            hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})" if hkl else "None"
            return f"TextureAnalysisResult(Preferred: {hkl_str}, σ={self.degree_of_texture:.2f})"


class TextureAnalyzer:
    """
    Harris Texture Coefficient analyzer.
    
    The Texture Coefficient (TC) is defined as:
    
    TC(hkl) = [I(hkl) / I₀(hkl)] / [1/n × Σ I(hkl) / I₀(hkl)]
    
    where:
    - I(hkl): Observed intensity of (hkl) reflection
    - I₀(hkl): Standard intensity from JCPDS
    - n: Number of reflections
    
    Interpretation:
    - TC = 1 for all peaks: Random orientation (powder average)
    - TC > 1: Preferred orientation along this direction
    - TC < 1: Under-represented direction
    - Sum of all TC = n (number of reflections)
    
    The degree of texture σ is:
    σ = √[Σ(TC - 1)² / n]
    
    σ ≈ 0: Random texture
    σ >> 0: Strong preferred orientation
    """
    
    # Threshold for considering preferred orientation
    TC_PREFERRED_THRESHOLD = 1.5
    
    # Threshold for random texture
    RANDOM_TEXTURE_SIGMA = 0.3
    
    def __init__(
        self,
        standard_data: Optional[Dict] = None
    ):
        """
        Initialize texture analyzer.
        
        Args:
            standard_data: Custom JCPDS data, defaults to Cu standard
        """
        self.standard_data = standard_data or CU_JCPDS_STANDARD
    
    def analyze(
        self,
        intensities: Dict[Tuple[int, int, int], float],
        normalize: bool = True
    ) -> TextureAnalysisResult:
        """
        Perform texture coefficient analysis.
        
        Args:
            intensities: Dict mapping (hkl) tuples to observed intensities
            normalize: If True, normalize intensities to max=100
            
        Returns:
            TextureAnalysisResult with TC values and interpretation
        """
        # Filter to only include peaks present in both observed and standard
        common_hkl = set(intensities.keys()) & set(self.standard_data.keys())
        
        if len(common_hkl) < 2:
            raise ValueError(
                f"Need at least 2 common peaks for TC analysis. "
                f"Found: {len(common_hkl)}"
            )
        
        n = len(common_hkl)
        
        # Calculate I(hkl) / I₀(hkl) ratios
        ratios = {}
        for hkl in common_hkl:
            i_obs = intensities[hkl]
            i_std = self.standard_data[hkl]["intensity"]
            ratios[hkl] = i_obs / i_std
        
        # Calculate average ratio
        ratio_sum = sum(ratios.values())
        avg_ratio = ratio_sum / n
        
        # Calculate TC for each peak
        tc_values = {}
        details = []
        
        for hkl in common_hkl:
            tc = ratios[hkl] / avg_ratio
            tc_values[hkl] = tc
            
            is_preferred = tc > self.TC_PREFERRED_THRESHOLD
            
            details.append(TextureResult(
                hkl=hkl,
                tc_value=tc,
                intensity_observed=intensities[hkl],
                intensity_standard=self.standard_data[hkl]["intensity"],
                is_preferred=is_preferred
            ))
        
        # Sort details by TC value (descending)
        details.sort(key=lambda x: x.tc_value, reverse=True)
        
        # Find preferred orientation (highest TC)
        preferred_hkl = max(tc_values.keys(), key=lambda k: tc_values[k])
        max_tc = tc_values[preferred_hkl]
        
        # Calculate degree of texture σ
        sigma_sq = sum((tc - 1)**2 for tc in tc_values.values()) / n
        sigma = np.sqrt(sigma_sq)
        
        # Determine if texture is random
        is_random = sigma < self.RANDOM_TEXTURE_SIGMA
        
        return TextureAnalysisResult(
            tc_values=tc_values,
            preferred_orientation=preferred_hkl if not is_random else None,
            degree_of_texture=sigma,
            is_random=is_random,
            details=details
        )
    
    def analyze_from_peaks(
        self,
        peaks: List[Tuple[float, float]],
        hkl_assignments: List[Tuple[int, int, int]]
    ) -> TextureAnalysisResult:
        """
        Analyze texture from peak fitting results.
        
        Args:
            peaks: List of (two_theta, intensity) tuples
            hkl_assignments: List of (h, k, l) assignments for each peak
            
        Returns:
            TextureAnalysisResult
        """
        if len(peaks) != len(hkl_assignments):
            raise ValueError("peaks and hkl_assignments must have same length")
        
        intensities = {}
        for (two_theta, intensity), hkl in zip(peaks, hkl_assignments):
            intensities[hkl] = intensity
        
        return self.analyze(intensities)
    
    def get_hkl_for_angle(
        self,
        two_theta: float,
        tolerance: float = 1.0
    ) -> Optional[Tuple[int, int, int]]:
        """
        Find (hkl) assignment for a given 2θ angle.
        
        Args:
            two_theta: Observed 2θ position (degrees)
            tolerance: Maximum deviation from standard (degrees)
            
        Returns:
            (h, k, l) tuple or None if no match
        """
        best_match = None
        min_diff = tolerance
        
        for hkl, data in self.standard_data.items():
            diff = abs(two_theta - data["two_theta"])
            if diff < min_diff:
                min_diff = diff
                best_match = hkl
        
        return best_match


def get_standard_angles(material: str = "Cu") -> Dict[Tuple[int, int, int], float]:
    """
    Get standard 2θ angles for common materials.
    
    Args:
        material: Material name ("Cu" supported)
        
    Returns:
        Dict mapping (hkl) to 2θ angle (degrees)
    """
    if material.upper() == "CU":
        return {
            hkl: data["two_theta"] 
            for hkl, data in CU_JCPDS_STANDARD.items()
        }
    else:
        raise ValueError(f"Unknown material: {material}")


def calculate_texture_coefficient(
    observed_intensities: Dict[Tuple[int, int, int], float],
    standard_intensities: Optional[Dict[Tuple[int, int, int], float]] = None
) -> Dict[Tuple[int, int, int], float]:
    """
    Convenience function for texture coefficient calculation.
    
    Args:
        observed_intensities: Dict of {(hkl): intensity}
        standard_intensities: Dict of {(hkl): standard_intensity}
                             If None, uses Cu JCPDS data
        
    Returns:
        Dict of {(hkl): TC_value}
    """
    if standard_intensities is None:
        standard_data = CU_JCPDS_STANDARD
    else:
        standard_data = {
            hkl: {"intensity": i} 
            for hkl, i in standard_intensities.items()
        }
    
    analyzer = TextureAnalyzer(standard_data)
    result = analyzer.analyze(observed_intensities)
    return result.tc_values
