"""
Defect and Stress Analysis Module
=================================

Stacking fault detection, lattice constant monitoring, and residual stress analysis.

Reference: 計劃書/07_微結構分析_下篇_缺陷與應力.md
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict
from enum import Enum


# =============================================================================
# Constants (文件 07 §10)
# =============================================================================

# Standard peak separation (111)-(200)
STANDARD_PEAK_SEPARATION = 7.136  # degrees

# Warren geometric coefficient for FCC (111)-(200) stacking fault
# α = (Δ2θ_exp - Δ2θ_std) / G
# Based on design document (Doc 07 §3.3): 0.136° deviation → α ≈ 0.68%
# Empirical relation: every 0.2° deviation corresponds to α = 1%
# Therefore: α = deviation / G → 0.01 = -0.2 / G → G = -20
WARREN_G_COEFFICIENT = -20.0

# Standard lattice constant for Cu
STANDARD_LATTICE_CONSTANT = 3.6150  # Å

# Lattice constant thresholds
LATTICE_MINOR_THRESHOLD = 3.616  # Å
LATTICE_SEVERE_THRESHOLD = 3.618  # Å

# Wavelength
WAVELENGTH_CU_KA1 = 1.54056  # Å

# Poisson's ratio for Cu
POISSON_RATIO_CU = 0.34


# =============================================================================
# Enums
# =============================================================================

class StackingFaultSeverity(Enum):
    """Stacking fault severity classification."""
    NORMAL = "normal"           # α < 0.5%
    MILD = "mild"               # 0.5% ≤ α < 1%
    MODERATE = "moderate"       # 1% ≤ α < 2%
    SEVERE = "severe"           # α ≥ 2%


class LatticeStatus(Enum):
    """Lattice constant status classification."""
    NORMAL = "normal"           # a ≈ 3.615 Å
    MINOR_EXPANSION = "minor"   # 3.616-3.618 Å
    SEVERE_EXPANSION = "severe" # > 3.618 Å


class AnnealingState(Enum):
    """Self-annealing state classification."""
    AS_DEPOSITED = "as-deposited"   # < 1 hour
    PARTIAL = "partial"             # 1-24 hours
    ANNEALED = "annealed"           # 1-7 days
    STABLE = "stable"               # > 7 days
    UNKNOWN = "unknown"             # Not specified


class StressType(Enum):
    """Residual stress type."""
    TENSILE = "tensile"         # d > d₀
    COMPRESSIVE = "compressive" # d < d₀
    NEUTRAL = "neutral"         # d ≈ d₀


# =============================================================================
# Result Dataclasses
# =============================================================================

@dataclass
class StackingFaultResult:
    """Result of stacking fault analysis."""
    peak_separation_deg: float
    deviation_deg: float
    alpha_probability: float  # Stacking fault probability (0-1)
    alpha_percent: float      # As percentage
    severity: StackingFaultSeverity
    sps_warning: bool
    message: str
    
    def __repr__(self) -> str:
        return (
            f"StackingFault: Δ2θ={self.peak_separation_deg:.3f}°, "
            f"α={self.alpha_percent:.2f}% [{self.severity.value}]"
        )


@dataclass
class LatticeConstantResult:
    """Result of lattice constant analysis."""
    lattice_constant: float   # Å
    deviation: float          # From standard
    d_spacing: float          # Å
    hkl_used: Tuple[int, int, int]
    two_theta_used: float
    status: LatticeStatus
    message: str
    
    def __repr__(self) -> str:
        return (
            f"Lattice: a={self.lattice_constant:.4f} Å "
            f"(Δ={self.deviation:+.4f} Å) [{self.status.value}]"
        )


@dataclass
class ResidualStressResult:
    """Result of residual stress estimation."""
    stress_mpa: float
    stress_type: StressType
    d_measured: float
    d_standard: float
    hkl: Tuple[int, int, int]
    youngs_modulus_gpa: float
    message: str


@dataclass
class DefectAnalysisResult:
    """Complete defect and stress analysis result."""
    # Stacking fault
    stacking_fault: Optional[StackingFaultResult] = None
    
    # Lattice constant
    lattice_constant: Optional[LatticeConstantResult] = None
    
    # Residual stress
    residual_stress: Optional[ResidualStressResult] = None
    
    # Self-annealing
    annealing_state: AnnealingState = AnnealingState.UNKNOWN
    sample_age_hours: Optional[float] = None
    annealing_note: str = ""
    
    # Overall
    overall_status: str = ""
    recommendations: List[str] = field(default_factory=list)


# =============================================================================
# Module Q: Stacking Fault Analyzer
# =============================================================================

class StackingFaultAnalyzer:
    """
    Warren-based stacking fault analyzer.
    
    Reference: 文件 07 §2-3
    
    For FCC metals with intrinsic stacking faults:
    - (111) peak shifts to higher angles
    - (200) peak shifts to lower angles
    - Peak separation decreases
    
    Warren formula:
        α = (Δ2θ_exp - Δ2θ_std) / G
        
    where G ≈ -20.0 for FCC (111)-(200) when using degrees
    """
    
    def __init__(self, g_coefficient: float = WARREN_G_COEFFICIENT):
        """
        Initialize stacking fault analyzer.
        
        Args:
            g_coefficient: Warren geometric coefficient (default: -20.0)
        """
        self.g_coefficient = g_coefficient
    
    def analyze(
        self,
        two_theta_111: float,
        two_theta_200: float
    ) -> StackingFaultResult:
        """
        Analyze stacking fault probability.
        
        Args:
            two_theta_111: Measured 2θ position of (111) peak
            two_theta_200: Measured 2θ position of (200) peak
            
        Returns:
            StackingFaultResult with α probability
        """
        # Q.1.3: Calculate peak separation
        peak_separation = two_theta_200 - two_theta_111
        
        # Q.1.4: Calculate deviation from standard
        deviation = peak_separation - STANDARD_PEAK_SEPARATION
        
        # Q.2.1-2.2: Warren formula for α (文件 07 §3.1)
        # α = deviation / G
        # G is negative, deviation is negative when SF present
        # So α comes out positive
        alpha = deviation / self.g_coefficient if self.g_coefficient != 0 else 0
        alpha = max(0, alpha)  # α cannot be negative
        
        alpha_percent = alpha * 100
        
        # Determine severity (文件 07 §4.1)
        if alpha_percent < 0.5:
            severity = StackingFaultSeverity.NORMAL
        elif alpha_percent < 1.0:
            severity = StackingFaultSeverity.MILD
        elif alpha_percent < 2.0:
            severity = StackingFaultSeverity.MODERATE
        else:
            severity = StackingFaultSeverity.SEVERE
        
        # Q.2.3: SPS warning (文件 07 §4.2)
        sps_warning = peak_separation < 7.0
        
        # Generate message
        if sps_warning:
            message = (
                f"警告：峰間距 {peak_separation:.3f}° < 7.0°，"
                "可能 SPS 加速劑濃度過高"
            )
        elif severity == StackingFaultSeverity.NORMAL:
            message = "堆垛層錯在正常範圍內"
        else:
            message = f"檢測到堆垛層錯 α ≈ {alpha_percent:.2f}%"
        
        return StackingFaultResult(
            peak_separation_deg=peak_separation,
            deviation_deg=deviation,
            alpha_probability=alpha,
            alpha_percent=alpha_percent,
            severity=severity,
            sps_warning=sps_warning,
            message=message
        )


# =============================================================================
# Module R: Lattice Constant Calculator
# =============================================================================

def calculate_d_spacing(
    two_theta: float,
    wavelength: float = WAVELENGTH_CU_KA1
) -> float:
    """
    Calculate d-spacing from Bragg's law.
    
    d = λ / (2 sin θ)
    """
    theta_rad = two_theta / 2 * np.pi / 180
    return wavelength / (2 * np.sin(theta_rad))


def calculate_lattice_constant(
    two_theta: float,
    hkl: Tuple[int, int, int],
    wavelength: float = WAVELENGTH_CU_KA1
) -> float:
    """
    Calculate lattice constant from peak position.
    
    a = d × √(h² + k² + l²)
    
    Reference: 文件 07 §5.3
    """
    d = calculate_d_spacing(two_theta, wavelength)
    h, k, l = hkl
    return d * np.sqrt(h**2 + k**2 + l**2)


class LatticeMonitor:
    """
    Lattice constant and residual stress monitor.
    
    Reference: 文件 07 §5
    
    Prefers high-angle peaks (311, 220) for better accuracy.
    """
    
    def __init__(
        self,
        standard_a: float = STANDARD_LATTICE_CONSTANT,
        wavelength: float = WAVELENGTH_CU_KA1
    ):
        self.standard_a = standard_a
        self.wavelength = wavelength
        
        # Standard d-spacings for Cu
        self.standard_d = {
            (1, 1, 1): 2.0871,
            (2, 0, 0): 1.8075,
            (2, 2, 0): 1.2780,
            (3, 1, 1): 1.0902,
        }
    
    def analyze_lattice(
        self,
        two_theta: float,
        hkl: Tuple[int, int, int]
    ) -> LatticeConstantResult:
        """
        Analyze lattice constant from peak position.
        
        Args:
            two_theta: Measured 2θ position
            hkl: Miller indices of the peak
            
        Returns:
            LatticeConstantResult
        """
        # R.1.2-3: Calculate d and a
        d = calculate_d_spacing(two_theta, self.wavelength)
        a = calculate_lattice_constant(two_theta, hkl, self.wavelength)
        
        # Calculate deviation
        deviation = a - self.standard_a
        
        # R.1.4: Determine status
        if a <= LATTICE_MINOR_THRESHOLD:
            status = LatticeStatus.NORMAL
            message = "晶格常數在正常範圍"
        elif a <= LATTICE_SEVERE_THRESHOLD:
            status = LatticeStatus.MINOR_EXPANSION
            message = "輕微晶格擴張，可能存在雜質或應力"
        else:
            status = LatticeStatus.SEVERE_EXPANSION
            message = "嚴重晶格擴張，需檢查添加劑純度或雜質固溶"
        
        return LatticeConstantResult(
            lattice_constant=a,
            deviation=deviation,
            d_spacing=d,
            hkl_used=hkl,
            two_theta_used=two_theta,
            status=status,
            message=message
        )
    
    def estimate_residual_stress(
        self,
        two_theta: float,
        hkl: Tuple[int, int, int],
        youngs_modulus_gpa: float
    ) -> ResidualStressResult:
        """
        Estimate residual stress from peak shift.
        
        σ = E / (1+ν) × (d - d₀) / d₀
        
        Reference: 文件 07 §7.2
        """
        d_measured = calculate_d_spacing(two_theta, self.wavelength)
        d_standard = self.standard_d.get(hkl, d_measured)
        
        delta_d_ratio = (d_measured - d_standard) / d_standard
        
        # Convert GPa to MPa
        E_mpa = youngs_modulus_gpa * 1000
        
        stress_mpa = E_mpa / (1 + POISSON_RATIO_CU) * delta_d_ratio
        
        # Determine stress type
        if stress_mpa > 10:
            stress_type = StressType.TENSILE
            message = f"拉伸應力 {stress_mpa:.0f} MPa"
        elif stress_mpa < -10:
            stress_type = StressType.COMPRESSIVE
            message = f"壓縮應力 {abs(stress_mpa):.0f} MPa"
        else:
            stress_type = StressType.NEUTRAL
            message = "應力接近中性"
        
        return ResidualStressResult(
            stress_mpa=stress_mpa,
            stress_type=stress_type,
            d_measured=d_measured,
            d_standard=d_standard,
            hkl=hkl,
            youngs_modulus_gpa=youngs_modulus_gpa,
            message=message
        )


# =============================================================================
# Self-Annealing State Machine
# =============================================================================

def determine_annealing_state(
    sample_age_hours: Optional[float] = None,
    fwhm_narrowing_detected: bool = False
) -> Tuple[AnnealingState, str]:
    """
    Determine self-annealing state from sample age.
    
    Reference: 文件 07 §6.4
    
    Args:
        sample_age_hours: Time since deposition (hours), None for unknown
        fwhm_narrowing_detected: True if FWHM is narrower than expected
        
    Returns:
        Tuple of (AnnealingState, recommendation_note)
    """
    if sample_age_hours is None:
        return (
            AnnealingState.UNKNOWN,
            "樣品存放時間未知，無法判定自退火狀態"
        )
    
    if sample_age_hours < 1:
        return (
            AnnealingState.AS_DEPOSITED,
            "鍍態樣品：預期細晶粒、高缺陷密度。建議 7 天後重測以獲得穩定結構"
        )
    elif sample_age_hours < 24:
        state = AnnealingState.PARTIAL
        if fwhm_narrowing_detected:
            note = "自退火進行中：已觀察到 FWHM 窄化"
        else:
            note = "自退火初期：晶粒可能正在長大"
        return (state, note)
    elif sample_age_hours < 168:  # 7 days
        return (
            AnnealingState.ANNEALED,
            "自退火中：結構尚未完全穩定"
        )
    else:
        return (
            AnnealingState.STABLE,
            "穩定態：結構已達穩定，分析結果可靠"
        )


# =============================================================================
# Convenience Functions
# =============================================================================

def analyze_stacking_faults(
    two_theta_111: float,
    two_theta_200: float
) -> StackingFaultResult:
    """Convenience function for stacking fault analysis."""
    analyzer = StackingFaultAnalyzer()
    return analyzer.analyze(two_theta_111, two_theta_200)


def analyze_lattice(
    two_theta: float,
    hkl: Tuple[int, int, int]
) -> LatticeConstantResult:
    """Convenience function for lattice constant analysis."""
    monitor = LatticeMonitor()
    return monitor.analyze_lattice(two_theta, hkl)
