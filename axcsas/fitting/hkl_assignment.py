"""
HKL Peak Assignment Module 峰位指標指派模組
==========================================

Automatic assignment of XRD peaks to copper crystallographic indices.
XRD 峰自動指派到銅的結晶學指標。
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict

from axcsas.core.copper_crystal import CU_JCPDS_EXTENDED, is_fcc_allowed


@dataclass
class PeakAssignment:
    """
    Result of hkl peak assignment.
    
    Attributes:
        hkl: Miller indices (h, k, l) or None if unassigned
        measured_two_theta: Measured 2θ position
        standard_two_theta: JCPDS standard position
        deviation: Difference (measured - standard) in degrees
        confidence: Assignment confidence ("high", "medium", "low")
    """
    hkl: Optional[Tuple[int, int, int]]
    measured_two_theta: float
    standard_two_theta: Optional[float]
    deviation: float
    confidence: str
    
    def __repr__(self) -> str:
        if self.hkl:
            return f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]}) @ {self.measured_two_theta:.3f}°"
        return f"Unassigned @ {self.measured_two_theta:.3f}°"


# JCPDS 04-0836 standard peak positions for Cu / 銅的 JCPDS 標準峰位
# Generated from CU_JCPDS constants / 從 CU_JCPDS 常數生成
JCPDS_COPPER_PEAKS: Dict[Tuple[int, int, int], float] = {
    hkl: data["two_theta"] for hkl, data in CU_JCPDS_EXTENDED.items()
}


def assign_hkl(
    two_theta: float,
    tolerance: float = 0.5
) -> Optional[Tuple[int, int, int]]:
    """
    Assign (hkl) Miller indices to a peak position.
    將 (hkl) Miller 指標指派給峰位。
    
    The assignment uses JCPDS 04-0836 standard values with tolerance
    for residual stress shifts (|Δ2θ| < 0.5° is normal for ED-Cu).
    使用 JCPDS 04-0836 標準值，容許殘留應力偏移 (|Δ2θ| < 0.5° 對 ED-Cu 正常)。
    
    Args:
        two_theta: Measured 2θ position in degrees
        tolerance: Maximum allowed deviation from standard (default: 0.5°)
        
    Returns:
        (h, k, l) tuple if matched, None otherwise
        
    Examples:
        >>> assign_hkl(43.32)
        (1, 1, 1)
        >>> assign_hkl(50.5)
        (2, 0, 0)
        >>> assign_hkl(30.0)
        None
    """
    best_match = None
    best_deviation = float('inf')
    
    for hkl, standard_pos in JCPDS_COPPER_PEAKS.items():
        deviation = abs(two_theta - standard_pos)
        if deviation <= tolerance and deviation < best_deviation:
            best_match = hkl
            best_deviation = deviation
    
    return best_match


def assign_hkl_detailed(
    two_theta: float,
    tolerance: float = 0.5
) -> PeakAssignment:
    """
    Assign (hkl) with detailed information including deviation and confidence.
    
    Confidence levels:
      - "high": |Δ2θ| < 0.2°
      - "medium": 0.2° ≤ |Δ2θ| < 0.4°
      - "low": 0.4° ≤ |Δ2θ| ≤ tolerance
    
    Args:
        two_theta: Measured 2θ position in degrees
        tolerance: Maximum allowed deviation
        
    Returns:
        PeakAssignment with full details
    """
    hkl = assign_hkl(two_theta, tolerance)
    
    if hkl is None:
        return PeakAssignment(
            hkl=None,
            measured_two_theta=two_theta,
            standard_two_theta=None,
            deviation=float('nan'),
            confidence="none"
        )
    
    standard_pos = JCPDS_COPPER_PEAKS[hkl]
    deviation = two_theta - standard_pos
    abs_dev = abs(deviation)
    
    # Determine confidence
    if abs_dev < 0.2:
        confidence = "high"
    elif abs_dev < 0.4:
        confidence = "medium"
    else:
        confidence = "low"
    
    return PeakAssignment(
        hkl=hkl,
        measured_two_theta=two_theta,
        standard_two_theta=standard_pos,
        deviation=deviation,
        confidence=confidence
    )


def assign_all_peaks(
    peak_positions: List[float],
    tolerance: float = 0.5
) -> List[PeakAssignment]:
    """
    Assign hkl to multiple peaks.
    
    Args:
        peak_positions: List of 2θ positions
        tolerance: Maximum deviation from standard
        
    Returns:
        List of PeakAssignment objects
    """
    return [assign_hkl_detailed(pos, tolerance) for pos in peak_positions]


def get_expected_peak_range(hkl: Tuple[int, int, int]) -> Optional[Tuple[float, float]]:
    """
    Get expected 2θ range for a given hkl.
    獲取給定 hkl 的預期 2θ 範圍。
    
    Args:
        hkl: Miller indices
        
    Returns:
        (min_2theta, max_2theta) or None if hkl not in database
    """
    ranges = {
        (1, 1, 1): (43.0, 43.6),
        (2, 0, 0): (50.2, 50.8),
        (2, 2, 0): (73.8, 74.4),
        (3, 1, 1): (89.6, 90.3),
    }
    return ranges.get(hkl)


def format_hkl(hkl: Optional[Tuple[int, int, int]]) -> str:
    """Format hkl tuple as string like '(111)'."""
    if hkl is None:
        return "(?)"
    return f"({hkl[0]}{hkl[1]}{hkl[2]})"
