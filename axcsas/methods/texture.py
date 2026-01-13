"""
Texture Analysis Module
=======================

Harris Texture Coefficient (TC) analysis for preferred orientation.
Harris 紋理係數（TC）分析用於優選取向。

NOTE: This module provides DATA ONLY. Interpretation should be done by humans.
注意：此模組僅提供資料，解讀應由人類進行。

Reference:
    Harris, G. B. (1952). Quantitative measurement of preferred orientation
    in rolled uranium bars. Phil. Mag., 43(336), 113-123.
"""

import numpy as np
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple

from axcsas.core.constants import CU_JCPDS


# =============================================================================
# Constants
# =============================================================================

# JCPDS standard data for Cu (PDF 04-0836)
JCPDS_STANDARD_INTENSITY: Dict[Tuple[int, int, int], float] = {
    (1, 1, 1): 100.0,
    (2, 0, 0): 46.0,
    (2, 2, 0): 20.0,
    (3, 1, 1): 17.0,
    (2, 2, 2): 5.0,
}

CU_JCPDS_STANDARD = {
    (1, 1, 1): {"two_theta": 43.298, "intensity": 100.0},
    (2, 0, 0): {"two_theta": 50.434, "intensity": 46.0},
    (2, 2, 0): {"two_theta": 74.130, "intensity": 20.0},
    (3, 1, 1): {"two_theta": 89.931, "intensity": 17.0},
    (2, 2, 2): {"two_theta": 95.139, "intensity": 5.0},
}

# TC interpretation thresholds (DATA MARKERS ONLY)
TC_RANDOM_MIN = 0.9
TC_RANDOM_MAX = 1.1
TC_PREFERRED_THRESHOLD = 1.5
RANDOM_TEXTURE_SIGMA = 0.3


# =============================================================================
# Enums
# =============================================================================

class OrientationType(Enum):
    """
    Classification of crystallographic orientation (DATA MARKER).
    結晶學取向分類（資料標記）。
    """
    PREFERRED = "preferred"    # TC > 1.0
    RANDOM = "random"          # TC ≈ 1.0 (0.9-1.1)
    SUPPRESSED = "suppressed"  # TC < 1.0


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class TCResult:
    """
    Single peak texture coefficient result.
    單一峰紋理係數結果。
    """
    hkl: Tuple[int, int, int]
    tc_value: float
    intensity_observed: float
    intensity_standard: float
    ratio: float = 0.0
    orientation_type: OrientationType = OrientationType.RANDOM

    def __repr__(self) -> str:
        hkl_str = f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]})"
        return f"{hkl_str}: TC={self.tc_value:.2f} [{self.orientation_type.value}]"


@dataclass
class TextureAnalysisResult:
    """
    Complete texture analysis result.
    完整紋理分析結果。

    NOTE: Contains DATA ONLY. No automatic process diagnosis.
    注意：僅包含資料，無自動過程診斷。

    Attributes:
        tc_values: TC for each (hkl). 每個 (hkl) 的 TC 值
        tc_details: Detailed results. 詳細結果
        dominant_hkl: Highest TC direction. 最高 TC 方向
        dominant_tc: Highest TC value. 最高 TC 值
        is_random: Whether texture is random. 是否為隨機紋理
        degree_of_texture: Standard deviation of TC values. TC 值標準差
        n_peaks: Number of peaks. 峰數
        intensity_type: Type of intensity used. 使用的強度類型
    """
    tc_values: Dict[Tuple[int, int, int], float] = field(default_factory=dict)
    tc_details: List[TCResult] = field(default_factory=list)
    dominant_hkl: Optional[Tuple[int, int, int]] = None
    dominant_tc: float = 1.0
    is_random: bool = True
    degree_of_texture: float = 0.0
    n_peaks: int = 0
    intensity_type: str = "area"

    def __repr__(self) -> str:
        if self.dominant_hkl:
            dom = f"({self.dominant_hkl[0]}{self.dominant_hkl[1]}{self.dominant_hkl[2]})"
        else:
            dom = "None"
        return f"TextureAnalysisResult(Dominant={dom}, TC={self.dominant_tc:.2f})"


# Backward compatibility alias
TextureResult = TCResult


# =============================================================================
# Texture Analyzer
# =============================================================================

class TextureAnalyzer:
    """
    Harris Texture Coefficient analyzer.
    Harris 紋理係數分析器。

    Provides DATA output only. Interpretation should be done by humans.
    僅提供資料輸出，解讀應由人類進行。

    Harris TC Formula:
        TC(hkl) = [I(hkl) / I₀(hkl)] / [1/N × Σ I(hkl) / I₀(hkl)]

    Interpretation:
        TC = 1 for all peaks: Random orientation (powder average)
        TC > 1: Preferred orientation along this direction
        TC < 1: Under-represented direction

    Args:
        standard_data: Custom JCPDS data (default: Cu)
        use_area: If True, intensities are integrated areas

    Example:
        >>> analyzer = TextureAnalyzer()
        >>> result = analyzer.analyze({(1,1,1): 15680, (2,0,0): 5520})
        >>> print(f"TC(111) = {result.tc_values[(1,1,1)]:.2f}")
    """

    def __init__(
        self,
        standard_data: Optional[Dict[Tuple[int, int, int], float]] = None,
        use_area: bool = True
    ) -> None:
        if standard_data is None:
            self.standard_data = JCPDS_STANDARD_INTENSITY
        else:
            self.standard_data = standard_data
        self.use_area = use_area
        self.intensity_type = "area" if use_area else "height"

    def analyze(
        self,
        intensities: Dict[Tuple[int, int, int], float],
        normalize: bool = True
    ) -> TextureAnalysisResult:
        """
        Perform texture coefficient analysis.
        執行紋理係數分析。

        Args:
            intensities: Dict mapping (hkl) to observed intensities.
                應為積分面積而非峰高。
            normalize: If True, normalize intensities (ignored in this version).

        Returns:
            TextureAnalysisResult with TC values (DATA ONLY).
        """
        # Filter to peaks with standard data
        valid_hkls = [hkl for hkl in intensities if hkl in self.standard_data]
        n_peaks = len(valid_hkls)

        if n_peaks < 2:
            return TextureAnalysisResult(
                tc_values={}, n_peaks=0,
                intensity_type=self.intensity_type
            )

        # Calculate intensity ratios I/I₀
        ratios: Dict[Tuple[int, int, int], float] = {}
        for hkl in valid_hkls:
            I_obs = intensities[hkl]
            I_std = self.standard_data[hkl]
            ratios[hkl] = I_obs / I_std

        # Calculate average ratio
        avg_ratio = sum(ratios.values()) / n_peaks

        # Calculate TC values
        tc_values: Dict[Tuple[int, int, int], float] = {}
        tc_details: List[TCResult] = []

        for hkl in valid_hkls:
            tc = ratios[hkl] / avg_ratio if avg_ratio > 0 else 1.0
            tc_values[hkl] = tc

            # Classify orientation type
            if TC_RANDOM_MIN <= tc <= TC_RANDOM_MAX:
                otype = OrientationType.RANDOM
            elif tc > 1.0:
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

        # Sort by TC value
        tc_details.sort(key=lambda x: x.tc_value, reverse=True)

        # Find dominant orientation
        dominant_hkl = max(tc_values, key=tc_values.get)
        dominant_tc = tc_values[dominant_hkl]

        # Check if random
        is_random = all(
            TC_RANDOM_MIN <= tc <= TC_RANDOM_MAX
            for tc in tc_values.values()
        )

        # Calculate degree of texture (σ)
        tc_list = list(tc_values.values())
        degree_of_texture = float(np.std(tc_list))

        return TextureAnalysisResult(
            tc_values=tc_values,
            tc_details=tc_details,
            dominant_hkl=dominant_hkl,
            dominant_tc=dominant_tc,
            is_random=is_random,
            degree_of_texture=degree_of_texture,
            n_peaks=n_peaks,
            intensity_type=self.intensity_type
        )

    def analyze_from_peaks(
        self,
        peaks: List[Tuple[float, float]],
        hkl_assignments: List[Tuple[int, int, int]]
    ) -> TextureAnalysisResult:
        """
        Analyze texture from peak fitting results.
        從峰擬合結果分析紋理。
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
        為給定的 2θ 角度查找 (hkl) 指派。
        """
        best_match = None
        min_diff = tolerance

        for hkl, data in CU_JCPDS_STANDARD.items():
            diff = abs(two_theta - data["two_theta"])
            if diff < min_diff:
                min_diff = diff
                best_match = hkl

        return best_match


# =============================================================================
# Convenience Functions
# =============================================================================

def analyze_texture(
    intensities: Dict[Tuple[int, int, int], float],
    use_area: bool = True
) -> TextureAnalysisResult:
    """
    Convenience function for texture analysis.
    紋理分析便利函式。

    Example:
        >>> intensities = {(1,1,1): 15680, (2,0,0): 5520, (2,2,0): 4200}
        >>> result = analyze_texture(intensities)
        >>> print(f"TC(111) = {result.tc_values[(1,1,1)]:.2f}")
    """
    analyzer = TextureAnalyzer(use_area=use_area)
    return analyzer.analyze(intensities)


def generate_texture_report(
    result: TextureAnalysisResult,
    sample_name: str = "Unknown"
) -> str:
    """
    Generate formatted texture analysis DATA report.
    產生格式化的紋理分析資料報告。

    NOTE: No automatic interpretation. Data only.
    """
    lines = [
        "=" * 70,
        "Texture Coefficient Analysis (Harris TC)",
        f"Sample: {sample_name}",
        f"Intensity Type: {result.intensity_type.upper()}",
        f"Number of Peaks: {result.n_peaks}",
        "=" * 70,
        "",
        f"{'Peak':^8} {'I_obs':>10} {'I_std':>8} {'I/I₀':>8} {'TC':>6} {'Type':^12}",
        "-" * 70,
    ]

    for detail in sorted(result.tc_details, key=lambda x: -x.tc_value):
        hkl_str = f"({detail.hkl[0]}{detail.hkl[1]}{detail.hkl[2]})"
        cls = detail.orientation_type.value.upper()
        lines.append(
            f"{hkl_str:^8} {detail.intensity_observed:>10.0f} "
            f"{detail.intensity_standard:>8.0f} {detail.ratio:>8.1f} "
            f"{detail.tc_value:>6.2f} {cls:^12}"
        )

    lines.append("-" * 70)

    if result.dominant_hkl:
        dom = f"({result.dominant_hkl[0]}{result.dominant_hkl[1]}{result.dominant_hkl[2]})"
    else:
        dom = "N/A"

    lines.extend([
        "",
        "Summary:",
        f"  Highest TC: {dom} = {result.dominant_tc:.2f}",
        f"  Degree of Texture (σ): {result.degree_of_texture:.3f}",
        f"  Random Texture: {result.is_random}",
        "=" * 70,
    ])

    return "\n".join(lines)


def get_standard_intensity(hkl: Tuple[int, int, int]) -> float:
    """
    Get JCPDS standard intensity for given hkl.
    取得指定 hkl 的 JCPDS 標準強度。

    Reference: JCPDS 04-0836
    """
    return JCPDS_STANDARD_INTENSITY.get(hkl, 0.0)


def get_standard_angles(material: str = "Cu") -> Dict[Tuple[int, int, int], float]:
    """
    Get standard 2θ angles for common materials.
    取得常見材料的標準 2θ 角度。
    """
    if material.upper() == "CU":
        return {hkl: data["two_theta"] for hkl, data in CU_JCPDS_STANDARD.items()}
    else:
        raise ValueError(f"Unknown material: {material}")


def calculate_texture_coefficient(
    observed_intensities: Dict[Tuple[int, int, int], float],
    standard_intensities: Optional[Dict[Tuple[int, int, int], float]] = None
) -> Dict[Tuple[int, int, int], float]:
    """
    Convenience function for texture coefficient calculation.
    紋理係數計算便利函式。
    """
    analyzer = TextureAnalyzer(standard_data=standard_intensities)
    result = analyzer.analyze(observed_intensities)
    return result.tc_values


# =============================================================================
# Backward Compatibility Aliases
# =============================================================================

TextureResultEnhanced = TextureAnalysisResult
TextureAnalyzerEnhanced = TextureAnalyzer
analyze_texture_enhanced = analyze_texture

