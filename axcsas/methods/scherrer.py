"""
Scherrer Crystallite Size Analysis
==================================

Implements the Scherrer equation with dynamic K values, validity flags,
and complete unit conversion safety.

Reference:
    Langford, J. I., & Wilson, A. J. C. (1978).
    Scherrer after sixty years. J. Appl. Cryst., 11, 102-113.
"""

import numpy as np
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple

from axcsas.core.constants import (
    CU_KA1,
    SCHERRER_K,
    MIN_RELIABLE_SIZE,
    MAX_RELIABLE_SIZE,
    MIN_BROADENING_RATIO,
)
from axcsas.core.copper_crystal import get_k_for_hkl
from axcsas.fitting.hkl_assignment import assign_hkl


# =============================================================================
# Constants
# =============================================================================

WAVELENGTH_CU_KA1 = CU_KA1  # Å
FWHM_RATIO_THRESHOLD = MIN_BROADENING_RATIO


# =============================================================================
# Validity Flag System
# =============================================================================

class ValidityFlag(Enum):
    """
    Scherrer calculation validity flags.
    Scherrer 計算有效性旗標。
    """
    VALID = "VALID"           # Normal calculation 正常計算
    UNRELIABLE = "UNRELIABLE" # FWHM ratio below threshold 寬化比值過低
    WARNING = "WARNING"       # Size exceeds limits 尺寸超出限制
    ERROR = "ERROR"           # Calculation failed 計算失敗


class GrainShape(Enum):
    """
    Grain shape for Scherrer constant selection.
    晶粒形狀，用於選擇 Scherrer 常數。
    """
    SPHERICAL = "spherical"   # 球形
    CUBIC = "cubic"           # 立方
    CUSTOM = "custom"         # 自訂


# =============================================================================
# Result Dataclass
# =============================================================================

@dataclass
class ScherrerResult:
    """
    Result from Scherrer crystallite size calculation.
    Scherrer 晶粒尺寸計算結果。

    Includes validity flags and complete metadata.
    包含有效性旗標與完整中繼資料。

    Attributes:
        size_nm: Crystallite size in nanometers. 晶粒尺寸（奈米）
        size_angstrom: Crystallite size in Ångströms. 晶粒尺寸（埃）
        two_theta: Peak position (2θ) in degrees. 峰位（度）
        hkl: Miller indices if assigned. Miller 指數
        k_factor: Scherrer constant used. 使用的 Scherrer 常數
        fwhm_observed: Observed FWHM in degrees. 觀測 FWHM（度）
        fwhm_instrumental: Instrumental FWHM in degrees. 儀器 FWHM（度）
        fwhm_sample: Sample broadening in degrees. 樣品寬化（度）
        fwhm_sample_rad: Sample broadening in radians. 樣品寬化（弧度）
        validity_flag: Calculation validity status. 計算有效性狀態
        warning_message: Warning or error message. 警告或錯誤訊息
        is_reliable: True if result is reliable. 結果是否可靠
    """
    size_nm: float
    size_angstrom: float
    two_theta: float
    hkl: Optional[Tuple[int, int, int]] = None
    k_factor: float = 0.89
    fwhm_observed: float = 0.0
    fwhm_instrumental: float = 0.0
    fwhm_sample: float = 0.0
    fwhm_sample_rad: float = 0.0
    validity_flag: ValidityFlag = ValidityFlag.VALID
    warning_message: str = ""
    is_reliable: bool = True

    def __repr__(self) -> str:
        hkl_str = f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]})" if self.hkl else "N/A"
        return (
            f"ScherrerResult(hkl={hkl_str}, size={self.size_nm:.1f} nm, "
            f"K={self.k_factor:.3f}, flag={self.validity_flag.value})"
        )


# =============================================================================
# Scherrer Calculator
# =============================================================================

class ScherrerCalculator:
    """
    Scherrer equation for crystallite size calculation.
    使用 Scherrer 方程式計算晶粒尺寸。

    D = K × λ / (β × cos θ)

    Features:
        - Dynamic K values based on hkl (cubic habit)
        - Validity flag system (VALID/UNRELIABLE/WARNING)
        - Unit conversion safety (degrees ↔ radians)
        - Caglioti instrumental broadening integration

    Args:
        wavelength: X-ray wavelength in Å (default: Cu Kα1 = 1.54056)
        use_cubic_habit: Use ED-Cu specific K values (default: True)
        caglioti_params: (U, V, W) tuple for instrumental broadening

    Example:
        >>> calc = ScherrerCalculator()
        >>> result = calc.calculate(43.32, 0.25, fwhm_instrumental=0.08)
        >>> print(f"D = {result.size_nm:.1f} nm")
        D = 49.0 nm
    """

    def __init__(
        self,
        wavelength: float = WAVELENGTH_CU_KA1,
        use_cubic_habit: bool = True,
        caglioti_params: Optional[Tuple[float, float, float]] = None
    ) -> None:
        self.wavelength = wavelength
        self.use_cubic_habit = use_cubic_habit
        self.caglioti_params = caglioti_params

    def calculate(
        self,
        two_theta: float,
        fwhm_observed: float,
        fwhm_instrumental: Optional[float] = None,
        hkl: Optional[Tuple[int, int, int]] = None
    ) -> ScherrerResult:
        """
        Calculate crystallite size using Scherrer equation.
        使用 Scherrer 方程式計算晶粒尺寸。

        Args:
            two_theta: Peak position in degrees (2θ). 峰位角度。
            fwhm_observed: Observed FWHM in degrees. 觀測半高寬。
            fwhm_instrumental: Instrumental FWHM in degrees (optional).
                儀器寬化，若 None 則使用 Caglioti。
            hkl: Miller indices for K selection (auto-assigned if None).
                Miller 指數，若 None 則自動指派。

        Returns:
            ScherrerResult with calculated size and metadata.
            包含計算尺寸與中繼資料的 ScherrerResult。

        Raises:
            ValueError: If fwhm_observed <= 0. 當 fwhm_observed <= 0。
        """
        warnings = []

        # Auto-assign hkl if not provided
        if hkl is None:
            hkl = assign_hkl(two_theta)

        # Get K factor
        if self.use_cubic_habit and hkl:
            k_factor = get_k_for_hkl(hkl[0], hkl[1], hkl[2])
        else:
            k_factor = SCHERRER_K.default

        # Calculate instrumental FWHM if Caglioti params available
        if fwhm_instrumental is None and self.caglioti_params:
            fwhm_instrumental = self._calculate_caglioti(two_theta)
        elif fwhm_instrumental is None:
            fwhm_instrumental = 0.0

        # Check validity threshold
        validity_flag = ValidityFlag.VALID
        is_reliable = True

        if fwhm_instrumental > 0:
            ratio = fwhm_observed / fwhm_instrumental
            if ratio < FWHM_RATIO_THRESHOLD:
                validity_flag = ValidityFlag.UNRELIABLE
                warnings.append(f"FWHM ratio {ratio:.2f} < {FWHM_RATIO_THRESHOLD}")
                is_reliable = False

        # Correct for instrumental broadening using quadratic subtraction:
        # β_sample² = β_observed² - β_instrumental²
        if fwhm_instrumental > 0:
            fwhm_sq_diff = fwhm_observed ** 2 - fwhm_instrumental ** 2
            
            if fwhm_sq_diff <= 0:
                # FWHM_obs <= FWHM_inst: physically impossible to extract sample broadening
                # This means crystallite size is beyond instrument detection limit
                validity_flag = ValidityFlag.UNRELIABLE
                warnings.append(
                    f"FWHM_obs ({fwhm_observed:.4f}°) ≤ FWHM_inst ({fwhm_instrumental:.4f}°): "
                    "crystallite size exceeds instrument detection limit"
                )
                is_reliable = False
                fwhm_sample = np.nan  # Mark as invalid
            else:
                fwhm_sample = np.sqrt(fwhm_sq_diff)
        else:
            # No instrumental correction available
            fwhm_sample = fwhm_observed

        # CRITICAL: Convert to radians
        theta_rad = np.radians(two_theta / 2)
        fwhm_sample_rad = np.radians(fwhm_sample)

        # Scherrer equation: D = K × λ / (β × cos θ)
        cos_theta = np.cos(theta_rad)
        size_angstrom = (k_factor * self.wavelength) / (fwhm_sample_rad * cos_theta)
        size_nm = size_angstrom / 10

        # Check size limits (only if valid calculation)
        if not np.isnan(size_nm):
            if size_nm > MAX_RELIABLE_SIZE:
                if validity_flag == ValidityFlag.VALID:
                    validity_flag = ValidityFlag.WARNING
                warnings.append(f"Size {size_nm:.1f} nm exceeds detection limit")
            elif size_nm < MIN_RELIABLE_SIZE:
                if validity_flag == ValidityFlag.VALID:
                    validity_flag = ValidityFlag.WARNING
                warnings.append(f"Size {size_nm:.1f} nm below precision limit")

        return ScherrerResult(
            size_nm=size_nm,
            size_angstrom=size_angstrom,
            two_theta=two_theta,
            hkl=hkl,
            k_factor=k_factor,
            fwhm_observed=fwhm_observed,
            fwhm_instrumental=fwhm_instrumental,
            fwhm_sample=fwhm_sample,
            fwhm_sample_rad=fwhm_sample_rad,
            validity_flag=validity_flag,
            warning_message="; ".join(warnings),
            is_reliable=is_reliable,
        )

    def _calculate_caglioti(self, two_theta: float) -> float:
        """
        Calculate instrumental FWHM using Caglioti equation.
        使用 Caglioti 方程式計算儀器 FWHM。

        FWHM²_inst = U·tan²θ + V·tanθ + W
        """
        if not self.caglioti_params:
            return 0.0

        U, V, W = self.caglioti_params
        theta_rad = np.radians(two_theta / 2)
        tan_theta = np.tan(theta_rad)

        fwhm_sq = U * tan_theta**2 + V * tan_theta + W
        return np.sqrt(max(fwhm_sq, 0.0))

    def batch_calculate(
        self,
        peaks: List[Tuple[float, float]],
        fwhm_instrumental: Optional[float] = None
    ) -> List[ScherrerResult]:
        """
        Calculate crystallite sizes for multiple peaks.
        批次計算多個峰的晶粒尺寸。

        Args:
            peaks: List of (two_theta, fwhm_observed) tuples.
            fwhm_instrumental: Common instrumental FWHM (optional).

        Returns:
            List of ScherrerResult objects.
        """
        return [
            self.calculate(two_theta, fwhm, fwhm_instrumental)
            for two_theta, fwhm in peaks
        ]

    def average_size(
        self,
        results: List[ScherrerResult],
        include_unreliable: bool = False
    ) -> Tuple[float, float]:
        """
        Calculate average crystallite size from multiple peaks.
        從多個峰計算平均晶粒尺寸。

        Args:
            results: List of ScherrerResult objects.
            include_unreliable: Include UNRELIABLE results in average.

        Returns:
            Tuple of (average_size_nm, std_dev_nm).
        """
        sizes = [
            r.size_nm for r in results
            if (r.is_reliable or include_unreliable) and r.validity_flag != ValidityFlag.ERROR
        ]

        if not sizes:
            return 0.0, 0.0

        return float(np.mean(sizes)), float(np.std(sizes))


# =============================================================================
# Convenience Functions
# =============================================================================

def calculate_crystallite_size(
    two_theta: float,
    fwhm: float,
    wavelength: float = WAVELENGTH_CU_KA1,
    k_factor: float = 0.89,
    fwhm_instrumental: float = 0.0
) -> float:
    """
    Quick crystallite size calculation.
    快速計算晶粒尺寸。

    Args:
        two_theta: Peak position (degrees). 峰位（度）
        fwhm: FWHM in degrees. 半高寬（度）
        wavelength: X-ray wavelength (Å), default Cu Kα1. X 射線波長
        k_factor: Scherrer constant, default 0.89. Scherrer 常數
        fwhm_instrumental: Instrumental FWHM (optional). 儀器寬化

    Returns:
        Crystallite size in nanometers. 晶粒尺寸（奈米）
    """
    calc = ScherrerCalculator(wavelength=wavelength, use_cubic_habit=False)
    result = calc.calculate(two_theta, fwhm, fwhm_instrumental)
    return result.size_nm


def calculate_scherrer(
    two_theta: float,
    fwhm_observed: float,
    fwhm_instrumental: float = 0.0,
    use_cubic_habit: bool = True
) -> ScherrerResult:
    """
    Convenience function for Scherrer calculation with full metadata.
    完整中繼資料的 Scherrer 計算便利函式。

    Example:
        >>> result = calculate_scherrer(43.32, 0.25, 0.08)
        >>> print(f"D = {result.size_nm:.1f} nm")
        D = 49.0 nm
    """
    calc = ScherrerCalculator(use_cubic_habit=use_cubic_habit)
    return calc.calculate(two_theta, fwhm_observed, fwhm_instrumental)


def generate_scherrer_report(
    results: List[ScherrerResult],
    sample_name: str = "Unknown"
) -> str:
    """
    Generate formatted Scherrer analysis report.
    產生格式化的 Scherrer 分析報告。
    """
    lines = [
        "=" * 70,
        "Scherrer Crystallite Size Analysis",
        f"Sample: {sample_name}",
        "=" * 70,
        "",
        f"{'Peak':^8} {'2θ (°)':>8} {'FWHM (°)':>10} {'K':>6} {'D (nm)':>10} {'Flag':>12}",
        "-" * 70,
    ]

    for r in results:
        hkl_str = f"({r.hkl[0]}{r.hkl[1]}{r.hkl[2]})" if r.hkl else "N/A"
        lines.append(
            f"{hkl_str:^8} {r.two_theta:>8.2f} {r.fwhm_sample:>10.4f} "
            f"{r.k_factor:>6.3f} {r.size_nm:>10.1f} {r.validity_flag.value:>12}"
        )

    lines.append("-" * 70)

    valid_sizes = [r.size_nm for r in results if r.is_reliable]
    if valid_sizes:
        avg = np.mean(valid_sizes)
        std = np.std(valid_sizes)
        lines.append(f"Average (reliable): {avg:.1f} ± {std:.1f} nm")

    lines.append("=" * 70)

    return "\n".join(lines)


# =============================================================================
# Backward Compatibility Aliases
# =============================================================================

# For backward compatibility with enhanced module
ScherrerResultEnhanced = ScherrerResult
ScherrerCalculatorEnhanced = ScherrerCalculator
calculate_scherrer_enhanced = calculate_scherrer

