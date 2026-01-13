"""
Williamson-Hall Analysis Module
===============================

Separates crystallite size and microstrain contributions to peak broadening.
分離晶粒尺寸與微應變對峰展寬的貢獻。

Reference:
    Williamson, G. K., & Hall, W. H. (1953).
    X-ray line broadening from filed aluminium and wolfram.
    Acta Metallurgica, 1(1), 22-31.
"""

import numpy as np
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple
from scipy.stats import linregress

from axcsas.core.constants import CU_KA1, SCHERRER_K


# =============================================================================
# Constants
# =============================================================================

WAVELENGTH_CU_KA1 = CU_KA1  # Å
WH_K_FACTOR = 0.9  # Average K value for W-H

# R² quality thresholds
R2_EXCELLENT = 0.95
R2_ACCEPTABLE = 0.85
MIN_PEAKS = 3

# Copper elastic anisotropy
MODULUS_MAP: Dict[Tuple[int, int, int], float] = {
    (1, 1, 1): 191.0,  # GPa, hardest direction
    (2, 0, 0): 67.0,   # GPa, softest direction
    (2, 2, 0): 130.0,  # GPa, intermediate
    (3, 1, 1): 128.0,  # GPa
}

ZENER_ANISOTROPY = 3.21


# =============================================================================
# Quality Level Enum
# =============================================================================

class WHQualityLevel(Enum):
    """
    Williamson-Hall analysis quality classification.
    Williamson-Hall 分析品質分類。
    """
    EXCELLENT = "excellent"     # R² > 0.95
    ACCEPTABLE = "acceptable"   # 0.85 ≤ R² ≤ 0.95
    POOR = "poor"               # R² < 0.85


# =============================================================================
# Result Dataclass
# =============================================================================

@dataclass
class WHResult:
    """
    Result from Williamson-Hall analysis.
    Williamson-Hall 分析結果。

    Attributes:
        crystallite_size_nm: Crystallite size in nm. 晶粒尺寸（奈米）
        microstrain: Dimensionless microstrain ε. 微應變
        r_squared: Goodness of fit R². 擬合優度
        intercept: y-intercept (Kλ/D). Y 軸截距
        slope: Slope (4ε). 斜率
        intercept_stderr: Standard error of intercept. 截距標準誤差
        slope_stderr: Standard error of slope. 斜率標準誤差
        size_error_nm: Error in size estimate. 尺寸誤差
        strain_error: Error in strain estimate. 應變誤差
        quality_level: Quality classification. 品質分類
        is_reliable: Whether fit is reliable. 是否可靠
        n_peaks: Number of peaks used. 使用的峰數
        peak_hkls: List of hkl indices. hkl 索引列表
        warning_message: Warning message. 警告訊息
        anisotropy_note: Anisotropy diagnostic. 異向性診斷
    """
    crystallite_size_nm: float
    microstrain: float
    r_squared: float
    intercept: float
    slope: float
    intercept_stderr: float = 0.0
    slope_stderr: float = 0.0
    size_error_nm: float = 0.0
    strain_error: float = 0.0
    quality_level: WHQualityLevel = WHQualityLevel.ACCEPTABLE
    is_reliable: bool = True
    n_peaks: int = 0
    peak_hkls: List[Tuple[int, int, int]] = field(default_factory=list)
    warning_message: str = ""
    anisotropy_note: str = ""

    def __repr__(self) -> str:
        return (
            f"WHResult(D={self.crystallite_size_nm:.1f} nm, "
            f"ε={self.microstrain:.2e}, R²={self.r_squared:.3f} "
            f"[{self.quality_level.value}])"
        )


# =============================================================================
# Williamson-Hall Analyzer
# =============================================================================

class WilliamsonHallAnalyzer:
    """
    Williamson-Hall analysis for separating size and strain broadening.
    用於分離尺寸與應變展寬的 Williamson-Hall 分析。

    W-H Equation:
        β cos θ = (K λ / D) + 4 ε sin θ

    where:
        β: Total peak broadening (FWHM in radians)
        θ: Bragg angle
        K: Scherrer constant
        λ: X-ray wavelength
        D: Crystallite size
        ε: Microstrain

    Features:
        - R² quality assessment with thresholds
        - Anisotropy diagnostics using MODULUS_MAP
        - Error propagation for D and ε

    Args:
        wavelength: X-ray wavelength in Å (default: Cu Kα1)
        k_factor: Scherrer constant for W-H (default: 0.9)

    Example:
        >>> analyzer = WilliamsonHallAnalyzer()
        >>> result = analyzer.analyze(two_theta, fwhm_sample)
        >>> print(f"D = {result.crystallite_size_nm:.1f} nm")
    """

    def __init__(
        self,
        wavelength: float = WAVELENGTH_CU_KA1,
        k_factor: float = WH_K_FACTOR
    ) -> None:
        self.wavelength = wavelength
        self.k_factor = k_factor

    def analyze(
        self,
        two_theta: np.ndarray,
        fwhm_sample: np.ndarray,
        hkl_list: Optional[List[Tuple[int, int, int]]] = None,
        fwhm_in_radians: bool = False
    ) -> WHResult:
        """
        Perform Williamson-Hall analysis.
        執行 Williamson-Hall 分析。

        Args:
            two_theta: Array of 2θ peak positions (degrees). 峰位陣列（度）
            fwhm_sample: Array of sample FWHM values. 樣品 FWHM 陣列
            hkl_list: Optional list of (h,k,l) for each peak. hkl 列表
            fwhm_in_radians: If True, FWHM is in radians. FWHM 是否為弧度

        Returns:
            WHResult with size, strain, and quality metrics.
        """
        two_theta = np.asarray(two_theta)
        fwhm_sample = np.asarray(fwhm_sample)
        n_peaks = len(two_theta)

        # Validate input
        if n_peaks < MIN_PEAKS:
            return self._create_failed_result(
                f"Need at least {MIN_PEAKS} peaks for W-H analysis", n_peaks
            )

        if len(fwhm_sample) != n_peaks:
            return self._create_failed_result(
                "2θ and FWHM arrays length mismatch", n_peaks
            )

        # Convert units
        theta_rad = np.radians(two_theta / 2.0)
        beta_rad = fwhm_sample if fwhm_in_radians else np.radians(fwhm_sample)

        # W-H coordinates: X = sin(θ), Y = β × cos(θ)
        x_data = np.sin(theta_rad)
        y_data = beta_rad * np.cos(theta_rad)

        # Linear regression
        reg_result = linregress(x_data, y_data)
        slope = reg_result.slope
        intercept = reg_result.intercept
        r_squared = reg_result.rvalue ** 2
        slope_stderr = reg_result.stderr
        intercept_stderr = getattr(reg_result, 'intercept_stderr', 0.0)

        # Calculate physical quantities: D = Kλ / intercept, ε = slope / 4
        if intercept > 0:
            size_angstrom = self.k_factor * self.wavelength / intercept
            size_nm = size_angstrom / 10.0
            size_error_nm = (size_angstrom / 10.0) * (intercept_stderr / intercept) if intercept_stderr > 0 else 0.0
        else:
            size_nm = float('inf')
            size_error_nm = 0.0

        microstrain = slope / 4.0
        strain_error = slope_stderr / 4.0 if slope_stderr else 0.0

        # Quality assessment
        quality_level, warning = self._assess_quality(r_squared)

        # Anisotropy diagnostics
        anisotropy_note = self._generate_anisotropy_note(r_squared, hkl_list)

        return WHResult(
            crystallite_size_nm=size_nm,
            microstrain=microstrain,
            r_squared=r_squared,
            intercept=intercept,
            slope=slope,
            intercept_stderr=intercept_stderr,
            slope_stderr=slope_stderr,
            size_error_nm=size_error_nm,
            strain_error=strain_error,
            quality_level=quality_level,
            is_reliable=(quality_level != WHQualityLevel.POOR),
            n_peaks=n_peaks,
            peak_hkls=hkl_list or [],
            warning_message=warning,
            anisotropy_note=anisotropy_note
        )

    def analyze_with_correction(
        self,
        two_theta: np.ndarray,
        fwhm_observed: np.ndarray,
        fwhm_instrumental: np.ndarray
    ) -> WHResult:
        """
        Perform W-H analysis with instrumental broadening correction.
        執行含儀器寬化校正的 W-H 分析。

        Uses geometric approximation: β_sample = β_obs - β²_inst / β_obs
        """
        fwhm_obs = np.asarray(fwhm_observed)
        fwhm_inst = np.asarray(fwhm_instrumental)
        fwhm_corrected = fwhm_obs - (fwhm_inst**2 / fwhm_obs)
        fwhm_corrected = np.maximum(fwhm_corrected, 0.001)
        return self.analyze(two_theta, fwhm_corrected)

    def get_plot_data(
        self,
        two_theta: np.ndarray,
        fwhm_sample: np.ndarray,
        fwhm_in_radians: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get data for W-H plot.
        取得 W-H 圖資料。

        Returns:
            Tuple of (x_data, y_data, x_fit_line, y_fit_line).
        """
        theta_rad = np.radians(np.asarray(two_theta) / 2.0)
        beta_rad = fwhm_sample if fwhm_in_radians else np.radians(np.asarray(fwhm_sample))

        x_data = np.sin(theta_rad)
        y_data = beta_rad * np.cos(theta_rad)

        reg = linregress(x_data, y_data)
        x_fit = np.linspace(x_data.min() * 0.9, x_data.max() * 1.1, 100)
        y_fit = reg.intercept + reg.slope * x_fit

        return x_data, y_data, x_fit, y_fit

    def _assess_quality(self, r_squared: float) -> Tuple[WHQualityLevel, str]:
        """Assess analysis quality based on R²."""
        if r_squared > R2_EXCELLENT:
            return WHQualityLevel.EXCELLENT, ""
        elif r_squared > R2_ACCEPTABLE:
            return WHQualityLevel.ACCEPTABLE, f"R² = {r_squared:.3f} (acceptable)"
        else:
            return WHQualityLevel.POOR, (
                f"R² = {r_squared:.3f} < {R2_ACCEPTABLE} - "
                "significant anisotropy detected"
            )

    def _generate_anisotropy_note(
        self,
        r_squared: float,
        hkl_list: Optional[List[Tuple[int, int, int]]]
    ) -> str:
        """Generate anisotropy diagnostic note."""
        if r_squared > R2_ACCEPTABLE:
            return ""

        lines = [
            "【Elastic Anisotropy Diagnostic 彈性異向性診斷】",
            f"Cu Zener ratio A = {ZENER_ANISOTROPY} (extreme anisotropy)",
            "",
            "Direction-dependent Young's modulus 各方向楊氏模數:"
        ]

        for hkl, E in sorted(MODULUS_MAP.items(), key=lambda x: -x[1]):
            hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
            lines.append(f"  E{hkl_str} = {E:.0f} GPa")

        lines.extend([
            "",
            f"E(111)/E(200) = {MODULUS_MAP[(1,1,1)]/MODULUS_MAP[(2,0,0)]:.1f}",
            "  → (200) strain broadening may be overestimated ~3x"
        ])

        return "\n".join(lines)

    def _create_failed_result(self, message: str, n_peaks: int) -> WHResult:
        """Create a failed result object."""
        return WHResult(
            crystallite_size_nm=0.0,
            microstrain=0.0,
            r_squared=0.0,
            intercept=0.0,
            slope=0.0,
            quality_level=WHQualityLevel.POOR,
            is_reliable=False,
            n_peaks=n_peaks,
            warning_message=message
        )


# =============================================================================
# Convenience Functions
# =============================================================================

def analyze_williamson_hall(
    two_theta: np.ndarray,
    fwhm_sample: np.ndarray,
    hkl_list: Optional[List[Tuple[int, int, int]]] = None
) -> WHResult:
    """
    Convenience function for Williamson-Hall analysis.
    Williamson-Hall 分析便利函式。

    Example:
        >>> two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        >>> fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        >>> result = analyze_williamson_hall(two_theta, fwhm)
        >>> print(f"D = {result.crystallite_size_nm:.1f} nm")
    """
    analyzer = WilliamsonHallAnalyzer()
    return analyzer.analyze(two_theta, fwhm_sample, hkl_list)


def generate_wh_report(result: WHResult, sample_name: str = "Unknown") -> str:
    """
    Generate formatted W-H analysis report.
    產生格式化的 W-H 分析報告。
    """
    lines = [
        "=" * 60,
        "Williamson-Hall Analysis",
        f"Sample: {sample_name}",
        "=" * 60,
        "",
        "Linear Regression:",
        f"  Y-intercept: {result.intercept:.5f}",
        f"  Slope: {result.slope:.5f}",
        f"  R²: {result.r_squared:.3f} ({result.quality_level.value.upper()})",
        "",
        "Results:",
        f"  Crystallite Size (D): {result.crystallite_size_nm:.1f} ± {result.size_error_nm:.1f} nm",
        f"  Microstrain (ε): {result.microstrain:.2e} ± {result.strain_error:.2e}",
        f"  Number of peaks: {result.n_peaks}",
    ]

    if result.warning_message:
        lines.extend(["", f"Note: {result.warning_message}"])

    if result.anisotropy_note:
        lines.extend(["", result.anisotropy_note])

    lines.append("=" * 60)
    return "\n".join(lines)


def get_modulus_for_hkl(hkl: Tuple[int, int, int]) -> float:
    """
    Get Young's modulus for a given hkl direction.
    取得指定 hkl 方向的楊氏模數。

    Returns:
        Young's modulus in GPa, or 120 GPa (average) if not found.
    """
    return MODULUS_MAP.get(hkl, 120.0)


# =============================================================================
# Backward Compatibility Aliases
# =============================================================================

WHResultEnhanced = WHResult
WilliamsonHallEnhanced = WilliamsonHallAnalyzer
analyze_williamson_hall_enhanced = analyze_williamson_hall
williamson_hall_analysis = analyze_williamson_hall
