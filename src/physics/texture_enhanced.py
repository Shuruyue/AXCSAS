"""
Enhanced Texture Analysis Module
================================

Harris Texture Coefficient (TC) analysis for electroplated copper.

NOTE: This module provides DATA ONLY. Interpretation should be done by humans.

Reference: 計劃書/06_微結構分析_上篇_織構.md

Physical Basis:
    Harris TC = [I(hkl) / I₀(hkl)] / [1/N × Σ I(hkl) / I₀(hkl)]
    
    TC = 1.0: Random orientation (powder average)
    TC > 1.0: Preferred orientation (grain alignment)
    TC < 1.0: Suppressed orientation
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict
from enum import Enum


# =============================================================================
# JCPDS Standard Intensity Library (文件 06 §3.1)
# =============================================================================

JCPDS_STANDARD_INTENSITY: Dict[Tuple[int, int, int], float] = {
    (1, 1, 1): 100.0,
    (2, 0, 0): 46.0,
    (2, 2, 0): 20.0,
    (3, 1, 1): 17.0,
    (2, 2, 2): 5.0,
}

# Standard intensity ratios (文件 06 §3.2)
STANDARD_RATIO_200_111 = 0.46  # I(200)/I(111)
STANDARD_RATIO_220_111 = 0.20  # I(220)/I(111)

# TC interpretation thresholds (DATA MARKERS ONLY)
TC_RANDOM_MIN = 0.9
TC_RANDOM_MAX = 1.1
TC_PREFERRED_THRESHOLD = 1.0


# =============================================================================
# Enums
# =============================================================================

class OrientationType(Enum):
    """Classification of crystallographic orientation (DATA MARKER)."""
    PREFERRED = "preferred"    # TC > 1.0
    RANDOM = "random"          # TC ≈ 1.0 (0.9-1.1)
    SUPPRESSED = "suppressed"  # TC < 1.0


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class TCResult:
    """Single peak texture coefficient result."""
    hkl: Tuple[int, int, int]
    tc_value: float
    intensity_observed: float
    intensity_standard: float
    ratio: float
    orientation_type: OrientationType
    
    def __repr__(self) -> str:
        hkl_str = f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]})"
        return f"{hkl_str}: TC={self.tc_value:.2f} [{self.orientation_type.value}]"


@dataclass
class TextureResultEnhanced:
    """
    Complete texture analysis result.
    
    NOTE: Contains DATA ONLY. No automatic process diagnosis.
    
    Reference: 文件 06 §8
    """
    # TC values
    tc_values: Dict[Tuple[int, int, int], float]
    tc_details: List[TCResult] = field(default_factory=list)
    
    # Dominant orientation (numerical fact, not interpretation)
    dominant_hkl: Optional[Tuple[int, int, int]] = None
    dominant_tc: float = 1.0
    
    # Statistical assessment
    is_random: bool = True
    degree_of_texture: float = 0.0  # Standard deviation of TC values
    
    # Metadata
    n_peaks: int = 0
    intensity_type: str = "area"  # "area" or "height"
    
    def __repr__(self) -> str:
        if self.dominant_hkl:
            dom = f"({self.dominant_hkl[0]}{self.dominant_hkl[1]}{self.dominant_hkl[2]})"
        else:
            dom = "None"
        return f"TextureResult: Dominant={dom} (TC={self.dominant_tc:.2f})"


# =============================================================================
# Enhanced Texture Analyzer
# =============================================================================

class TextureAnalyzerEnhanced:
    """
    Harris Texture Coefficient analyzer.
    
    Provides DATA output only. Process interpretation should be done by humans.
    
    Physical Background (文件 06 §2):
        - TC(hkl) indicates relative grain orientation along [hkl]
        - TC = 1 for perfect random powder
        - TC > 1 indicates preferred orientation
        - TC < 1 indicates suppressed orientation
    
    Harris TC Formula:
        TC(hkl) = [I(hkl) / I₀(hkl)] / [1/N × Σ I(hkl) / I₀(hkl)]
    """
    
    def __init__(
        self,
        standard_data: Optional[Dict[Tuple[int, int, int], float]] = None,
        use_area: bool = True
    ):
        """
        Initialize texture analyzer.
        
        Args:
            standard_data: Custom JCPDS data, defaults to Cu standard
            use_area: If True, intensities are integrated areas (recommended)
        """
        self.standard_data = standard_data or JCPDS_STANDARD_INTENSITY
        self.use_area = use_area
        self.intensity_type = "area" if use_area else "height"
    
    def analyze(
        self,
        intensities: Dict[Tuple[int, int, int], float]
    ) -> TextureResultEnhanced:
        """
        Perform texture coefficient analysis.
        
        Args:
            intensities: Dict mapping (hkl) tuples to observed intensities
                        (should be integrated areas, not peak heights)
            
        Returns:
            TextureResultEnhanced with TC values (DATA ONLY)
        """
        # Step 1: Filter to peaks with standard data
        valid_hkls = [hkl for hkl in intensities if hkl in self.standard_data]
        n_peaks = len(valid_hkls)
        
        if n_peaks < 2:
            return self._create_failed_result("需要至少 2 個峰進行織構分析")
        
        # Step 2: Calculate intensity ratios I/I₀
        ratios: Dict[Tuple[int, int, int], float] = {}
        for hkl in valid_hkls:
            I_obs = intensities[hkl]
            I_std = self.standard_data[hkl]
            ratios[hkl] = I_obs / I_std
        
        # Step 3: Calculate average ratio
        avg_ratio = sum(ratios.values()) / n_peaks
        
        # Step 4: Calculate TC values
        tc_values: Dict[Tuple[int, int, int], float] = {}
        tc_details: List[TCResult] = []
        
        for hkl in valid_hkls:
            tc = ratios[hkl] / avg_ratio if avg_ratio > 0 else 1.0
            tc_values[hkl] = tc
            
            # Classify orientation type (DATA MARKER)
            if TC_RANDOM_MIN <= tc <= TC_RANDOM_MAX:
                otype = OrientationType.RANDOM
            elif tc > TC_PREFERRED_THRESHOLD:
                otype = OrientationType.PREFERRED
            else:
                otype = OrientationType.SUPPRESSED
            
            tc_details.append(TCResult(
                hkl=hkl,
                tc_value=tc,
                intensity_observed=intensities[hkl],
                intensity_standard=self.standard_data[hkl],
                ratio=ratios[hkl],
                orientation_type=otype
            ))
        
        # Step 5: Find dominant orientation (statistical fact)
        dominant_hkl = max(tc_values, key=tc_values.get)
        dominant_tc = tc_values[dominant_hkl]
        
        # Step 6: Check if random
        is_random = all(
            TC_RANDOM_MIN <= tc <= TC_RANDOM_MAX 
            for tc in tc_values.values()
        )
        
        # Step 7: Calculate degree of texture
        tc_list = list(tc_values.values())
        degree_of_texture = float(np.std(tc_list))
        
        return TextureResultEnhanced(
            tc_values=tc_values,
            tc_details=tc_details,
            dominant_hkl=dominant_hkl,
            dominant_tc=dominant_tc,
            is_random=is_random,
            degree_of_texture=degree_of_texture,
            n_peaks=n_peaks,
            intensity_type=self.intensity_type
        )
    
    def _create_failed_result(self, message: str) -> TextureResultEnhanced:
        """Create a failed result object."""
        return TextureResultEnhanced(
            tc_values={},
            n_peaks=0
        )


# =============================================================================
# Convenience Functions
# =============================================================================

def analyze_texture_enhanced(
    intensities: Dict[Tuple[int, int, int], float],
    use_area: bool = True
) -> TextureResultEnhanced:
    """
    Convenience function for texture analysis.
    
    Example (文件 06 §4):
        >>> intensities = {
        ...     (1,1,1): 15680,
        ...     (2,0,0): 5520,
        ...     (2,2,0): 4200
        ... }
        >>> result = analyze_texture_enhanced(intensities)
        >>> print(f"TC(111) = {result.tc_values[(1,1,1)]:.2f}")
    """
    analyzer = TextureAnalyzerEnhanced(use_area=use_area)
    return analyzer.analyze(intensities)


def generate_texture_report(
    result: TextureResultEnhanced,
    sample_name: str = "Unknown"
) -> str:
    """
    Generate formatted texture analysis DATA report.
    
    NOTE: No automatic interpretation. Data only.
    """
    lines = [
        "=== Texture Coefficient Analysis (Harris TC) ===",
        f"Sample: {sample_name}",
        f"Intensity Type: {result.intensity_type.upper()}",
        f"Number of Peaks: {result.n_peaks}",
        "",
        "Texture Coefficients:",
        "| Peak   | I_obs    | I_std | I/I₀    | TC    | Classification |",
        "|--------|----------|-------|---------|-------|----------------|",
    ]
    
    for detail in sorted(result.tc_details, key=lambda x: -x.tc_value):
        hkl_str = f"({detail.hkl[0]}{detail.hkl[1]}{detail.hkl[2]})"
        cls = detail.orientation_type.value.upper()
        lines.append(
            f"| {hkl_str:6s} | {detail.intensity_observed:8.0f} | "
            f"{detail.intensity_standard:5.0f} | {detail.ratio:7.1f} | "
            f"{detail.tc_value:.2f}  | {cls:14s} |"
        )
    
    # Statistical summary
    if result.dominant_hkl:
        dom = f"({result.dominant_hkl[0]}{result.dominant_hkl[1]}{result.dominant_hkl[2]})"
    else:
        dom = "N/A"
    
    lines.extend([
        "",
        "Statistical Summary:",
        f"  Highest TC: {dom} = {result.dominant_tc:.2f}",
        f"  Degree of Texture (σ): {result.degree_of_texture:.3f}",
        f"  All TC ≈ 1.0: {result.is_random}",
    ])
    
    return "\n".join(lines)


def get_standard_intensity(hkl: Tuple[int, int, int]) -> float:
    """
    Get JCPDS standard intensity for given hkl.
    
    Reference: JCPDS 04-0836
    """
    return JCPDS_STANDARD_INTENSITY.get(hkl, 0.0)
