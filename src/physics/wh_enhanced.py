"""
Enhanced Williamson-Hall Analysis Module
========================================

Advanced W-H analysis with R² quality assessment, anisotropy diagnostics,
and integration with Phase 04 Scherrer results.

Reference: 計劃書/05_晶粒尺寸計算_下篇_WH分析.md
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict
from enum import Enum
from scipy.stats import linregress
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from core.copper_crystal import CU_ELASTIC


# =============================================================================
# Constants
# =============================================================================

# W-H uses average K value (文件 05 §7.4)
WH_K_FACTOR = 0.9
WAVELENGTH_CU_KA1 = 1.54056  # Å

# R² quality thresholds (文件 05 §8.1)
R2_EXCELLENT = 0.95
R2_ACCEPTABLE = 0.85

# Copper elastic anisotropy (文件 05 §4.2)
MODULUS_MAP: Dict[Tuple[int, int, int], float] = {
    (1, 1, 1): 191.0,  # GPa, hardest direction
    (2, 0, 0): 67.0,   # GPa, softest direction
    (2, 2, 0): 130.0,  # GPa, intermediate
    (3, 1, 1): 128.0,  # GPa
}

ZENER_ANISOTROPY = 3.21  # 文件 05 §3.3


# =============================================================================
# Quality Level Enum
# =============================================================================

class WHQualityLevel(Enum):
    """
    Williamson-Hall analysis quality classification.
    
    Reference: 文件 05 §8.1
    """
    EXCELLENT = "excellent"       # R² > 0.95
    ACCEPTABLE = "acceptable"     # 0.85 ≤ R² ≤ 0.95
    POOR = "poor"                 # R² < 0.85


# =============================================================================
# Enhanced W-H Result
# =============================================================================

@dataclass
class WHResultEnhanced:
    """
    Enhanced result from Williamson-Hall analysis.
    
    Includes quality assessment and anisotropy diagnostics.
    """
    # Physical results
    crystallite_size_nm: float
    microstrain: float
    
    # Regression parameters
    intercept: float
    slope: float
    r_squared: float
    
    # Error estimates
    intercept_stderr: float = 0.0
    slope_stderr: float = 0.0
    size_error_nm: float = 0.0
    strain_error: float = 0.0
    
    # Quality assessment
    quality_level: WHQualityLevel = WHQualityLevel.ACCEPTABLE
    is_reliable: bool = True
    
    # Data points info
    n_peaks: int = 0
    peak_hkls: List[Tuple[int, int, int]] = field(default_factory=list)
    
    # Warnings and diagnostics
    warning_message: str = ""
    anisotropy_note: str = ""
    
    def __repr__(self) -> str:
        return (
            f"WHResult: D = {self.crystallite_size_nm:.1f} nm, "
            f"ε = {self.microstrain:.2e}, "
            f"R² = {self.r_squared:.3f} [{self.quality_level.value}]"
        )


@dataclass
class WHDataPoint:
    """Single data point for W-H plot."""
    hkl: Optional[Tuple[int, int, int]]
    two_theta: float
    sin_theta: float
    beta_cos_theta: float
    fwhm_sample: float


# =============================================================================
# Enhanced W-H Analyzer
# =============================================================================

class WilliamsonHallEnhanced:
    """
    Enhanced Williamson-Hall analyzer for electrodeposited copper.
    
    Features:
    - R² quality assessment with clear thresholds
    - Anisotropy diagnostics using MODULUS_MAP
    - Integration with ScherrerResultEnhanced
    - Error propagation for D and ε
    
    Reference: 文件 05 (全章)
    
    W-H Equation:
        β cos θ = (K λ / D) + 4 ε sin θ
        
        Y = intercept + slope × X
        where:
            X = sin θ
            Y = β cos θ
    """
    
    def __init__(
        self,
        wavelength: float = WAVELENGTH_CU_KA1,
        k_factor: float = WH_K_FACTOR
    ):
        """
        Initialize enhanced W-H analyzer.
        
        Args:
            wavelength: X-ray wavelength in Å (default: Cu Kα1)
            k_factor: Scherrer constant for W-H (default: 0.9)
        """
        self.wavelength = wavelength
        self.k_factor = k_factor
    
    def analyze(
        self,
        two_theta: np.ndarray,
        fwhm_sample: np.ndarray,
        hkl_list: Optional[List[Tuple[int, int, int]]] = None,
        fwhm_in_radians: bool = False
    ) -> WHResultEnhanced:
        """
        Perform enhanced W-H analysis.
        
        Args:
            two_theta: Array of 2θ peak positions (degrees)
            fwhm_sample: Array of sample FWHM values (corrected for instrumental)
            hkl_list: Optional list of (h,k,l) for each peak
            fwhm_in_radians: If True, FWHM is already in radians
            
        Returns:
            WHResultEnhanced with size, strain, and quality metrics
        """
        # Step 1: Validate input
        n_peaks = len(two_theta)
        if n_peaks < 3:
            return self._create_failed_result(
                "需要至少 3 個峰進行 W-H 分析",
                n_peaks
            )
        
        if len(fwhm_sample) != n_peaks:
            return self._create_failed_result(
                "2θ 與 FWHM 數組長度不匹配",
                n_peaks
            )
        
        # Step 2: Convert units
        # K.1.2: θ_rad = 2θ / 2 * π / 180
        theta_rad = two_theta / 2.0 * np.pi / 180.0
        
        # K.1.4: β_rad (if not already in radians)
        if fwhm_in_radians:
            beta_rad = fwhm_sample
        else:
            beta_rad = fwhm_sample * np.pi / 180.0
        
        # Step 3: Calculate W-H coordinates
        # K.1.3: X = sin(θ)
        x_data = np.sin(theta_rad)
        
        # K.1.5: Y = β × cos(θ)
        y_data = beta_rad * np.cos(theta_rad)
        
        # Step 4: Linear regression
        reg_result = linregress(x_data, y_data)
        slope = reg_result.slope
        intercept = reg_result.intercept
        r_squared = reg_result.rvalue ** 2
        slope_stderr = reg_result.stderr
        intercept_stderr = reg_result.intercept_stderr if hasattr(reg_result, 'intercept_stderr') else 0.0
        
        # Step 5: Calculate physical quantities (文件 05 §7.4)
        # L.2.2: D = Kλ / intercept
        if intercept > 0:
            size_angstrom = self.k_factor * self.wavelength / intercept
            size_nm = size_angstrom / 10.0
            
            # Error propagation
            if intercept_stderr > 0:
                size_error_angstrom = size_angstrom * (intercept_stderr / intercept)
                size_error_nm = size_error_angstrom / 10.0
            else:
                size_error_nm = 0.0
        else:
            size_nm = float('inf')
            size_error_nm = 0.0
        
        # L.2.1: ε = slope / 4
        microstrain = slope / 4.0
        strain_error = slope_stderr / 4.0 if slope_stderr else 0.0
        
        # Step 6: Quality assessment (文件 05 §8.1)
        quality_level, warning = self._assess_quality(r_squared)
        
        # Step 7: Anisotropy diagnostics
        anisotropy_note = self._generate_anisotropy_note(r_squared, hkl_list)
        
        return WHResultEnhanced(
            crystallite_size_nm=size_nm,
            microstrain=microstrain,
            intercept=intercept,
            slope=slope,
            r_squared=r_squared,
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
    
    def _assess_quality(self, r_squared: float) -> Tuple[WHQualityLevel, str]:
        """
        Assess W-H analysis quality based on R².
        
        Reference: 文件 05 §8.1
        """
        if r_squared > R2_EXCELLENT:
            return WHQualityLevel.EXCELLENT, ""
        elif r_squared > R2_ACCEPTABLE:
            return WHQualityLevel.ACCEPTABLE, (
                f"R² = {r_squared:.3f} (輕微彈性異方性，結果可接受)"
            )
        else:
            return WHQualityLevel.POOR, (
                f"R² = {r_squared:.3f} < {R2_ACCEPTABLE} - "
                "檢測到顯著彈性異方性，建議使用 Modified W-H (MWH) 分析"
            )
    
    def _generate_anisotropy_note(
        self,
        r_squared: float,
        hkl_list: Optional[List[Tuple[int, int, int]]]
    ) -> str:
        """
        Generate anisotropy diagnostic note.
        
        Reference: 文件 05 §4
        """
        if r_squared > R2_ACCEPTABLE:
            return ""
        
        lines = [
            "【彈性異方性診斷】",
            f"銅的 Zener 異方性比 A = {ZENER_ANISOTROPY}，屬極端異方性材料",
            "",
            "各方向楊氏模數:"
        ]
        
        for hkl, E in sorted(MODULUS_MAP.items(), key=lambda x: -x[1]):
            hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
            lines.append(f"  E{hkl_str} = {E:.0f} GPa")
        
        lines.extend([
            "",
            f"影響：E(111)/E(200) = {MODULUS_MAP[(1,1,1)]/MODULUS_MAP[(2,0,0)]:.1f}",
            "  → (200) 峰的應變展寬可能被高估 ~3 倍"
        ])
        
        return "\n".join(lines)
    
    def _create_failed_result(self, message: str, n_peaks: int) -> WHResultEnhanced:
        """Create a failed result object."""
        return WHResultEnhanced(
            crystallite_size_nm=0.0,
            microstrain=0.0,
            intercept=0.0,
            slope=0.0,
            r_squared=0.0,
            quality_level=WHQualityLevel.POOR,
            is_reliable=False,
            n_peaks=n_peaks,
            warning_message=message
        )
    
    def get_plot_data(
        self,
        two_theta: np.ndarray,
        fwhm_sample: np.ndarray,
        fwhm_in_radians: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get data for W-H plot.
        
        Returns:
            Tuple of (x_data, y_data, x_fit_line, y_fit_line)
        """
        theta_rad = two_theta / 2.0 * np.pi / 180.0
        
        if fwhm_in_radians:
            beta_rad = fwhm_sample
        else:
            beta_rad = fwhm_sample * np.pi / 180.0
        
        x_data = np.sin(theta_rad)
        y_data = beta_rad * np.cos(theta_rad)
        
        # Fit line
        reg = linregress(x_data, y_data)
        x_fit = np.linspace(x_data.min() * 0.9, x_data.max() * 1.1, 100)
        y_fit = reg.intercept + reg.slope * x_fit
        
        return x_data, y_data, x_fit, y_fit


# =============================================================================
# Convenience Functions
# =============================================================================

def analyze_williamson_hall_enhanced(
    two_theta: np.ndarray,
    fwhm_sample: np.ndarray,
    hkl_list: Optional[List[Tuple[int, int, int]]] = None
) -> WHResultEnhanced:
    """
    Convenience function for enhanced W-H analysis.
    
    Example (文件 05 §7):
        >>> two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        >>> fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        >>> result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        >>> print(f"D = {result.crystallite_size_nm:.1f} nm, R² = {result.r_squared:.3f}")
    """
    analyzer = WilliamsonHallEnhanced()
    return analyzer.analyze(two_theta, fwhm_sample, hkl_list)


def generate_wh_report(result: WHResultEnhanced, sample_name: str = "Unknown") -> str:
    """
    Generate formatted W-H analysis report.
    
    Reference: 文件 05 §9
    """
    lines = [
        "=== Williamson-Hall Analysis ===",
        f"Sample: {sample_name}",
        "",
        "Linear Regression:",
        f"  Y-intercept: {result.intercept:.5f}",
        f"  Slope: {result.slope:.5f}",
        f"  R²: {result.r_squared:.3f} ({result.quality_level.value.upper()})",
        "",
        "Results:",
        f"  Crystallite Size (D): {result.crystallite_size_nm:.1f} ± {result.size_error_nm:.1f} nm",
        f"  Microstrain (ε): {result.microstrain:.2e} ± {result.strain_error:.2e}",
        f"  Number of peaks used: {result.n_peaks}",
    ]
    
    if result.warning_message:
        lines.extend(["", f"Note: {result.warning_message}"])
    
    if result.anisotropy_note:
        lines.extend(["", result.anisotropy_note])
    
    return "\n".join(lines)


def get_modulus_for_hkl(hkl: Tuple[int, int, int]) -> float:
    """
    Get Young's modulus for a given hkl direction.
    
    Reference: 文件 05 §4.2
    
    Returns:
        Young's modulus in GPa, or average (120 GPa) if not found
    """
    return MODULUS_MAP.get(hkl, 120.0)
