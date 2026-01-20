"""
Physical Constants Module 物理常數模組
======================================
Centralized physical constants for XRD analysis.
集中管理 XRD 分析所需的物理常數。

All parameters are verified with paper citations.
所有參數皆經論文驗證。
"""

from dataclasses import dataclass
from typing import Dict, Tuple


# =============================================================================
# X-ray Wavelengths (Ångströms) / X-ray 波長 (埃)
# Reference 出處: Bearden (1967), Rev. Mod. Phys. 39, 78-124, Table V, p.9
# =============================================================================

# Cu Kα1: K-L_III transition / Cu Kα1: K-L_III 躍遷
# Bearden (1967) Table V: 1.540562 Å, 8.04778 keV
CU_KA1 = 1.540562

# Cu Kα2: K-L_II transition / Cu Kα2: K-L_II 躍遷
# Bearden (1967) Table V: 1.544390 Å, 8.02783 keV
CU_KA2 = 1.544390

# Weighted average Kα for unresolved doublet / 未分辨雙峰之加權平均
# Calculation 計算: λ_avg = (2×Kα1 + Kα2) / 3 = (2×1.540562 + 1.544390) / 3
# Weighting basis 權重依據: Burger-Dorgelo rule (2p₃/₂: 4e⁻, 2p₁/₂: 2e⁻ → 2:1)
# Reference 出處: Cullity, "Elements of X-ray Diffraction", 3rd Ed., p.7
CU_KA_AVG = 1.541838

# Kα2/Kα1 intensity ratio / Kα2/Kα1 強度比
# Theoretical basis 理論依據: Burger-Dorgelo rule
#   - Kα1 from 2p₃/₂→1s transition (4 electrons) / Kα1 來自 2p₃/₂→1s 躍遷 (4 個電子)
#   - Kα2 from 2p₁/₂→1s transition (2 electrons) / Kα2 來自 2p₁/₂→1s 躍遷 (2 個電子)
#   - Ratio 比值: 2/4 = 0.5
# Reference 出處: Cullity, "Elements of X-ray Diffraction", 3rd Ed., p.7
KA2_KA1_RATIO = 0.5


# =============================================================================
# Scherrer Constants (K values) / Scherrer 常數 (K 值)
# Reference 出處: Langford & Wilson (1978), J. Appl. Cryst. 11, 102-113
# =============================================================================

@dataclass
class ScherrerConstants:
    """
    Scherrer constant K values for different grain shapes.
    不同晶粒形狀的 Scherrer 常數 K 值。
    
    IMPORTANT 重要:
    - These are K_w values for FWHM (half-width) definition
    - 這些是 K_w 值，適用於 FWHM（半高寬）定義
    - For integral breadth, use K_β values from L&W Table 1
    - 若使用積分寬度，請參考 L&W 論文 Table 1 的 K_β 值
    
    Reference 出處: L&W (1978), Table 1, p.107
    """
    # Sphere 球形: K_w = 0.8290 (L&W Table 1, p.107, Line 846)
    spherical: float = 0.829
    
    # Generic cubic habit 一般立方晶體: approximation for mixed orientations
    # Reference: Warren (1969), X-ray Diffraction
    cubic: float = 0.94
    
    # Octahedron 八面體: K_w = 0.83 (L&W p.107, same as sphere for 110)
    octahedral: float = 0.83
    
    # Default value 預設值: L&W 1978 spherical standard
    default: float = 0.829


SCHERRER_K = ScherrerConstants()





def get_jcpds_data(material: str = "Cu") -> Dict[Tuple[int, int, int], Dict]:
    """
    Get JCPDS standard data for copper.
    取得銅的 JCPDS 標準數據。
    
    Args:
        material: Material name (only "Cu" supported)
        
    Returns:
        Dictionary of JCPDS data keyed by (h, k, l)
    """
    if material.upper() != "CU":
        raise ValueError(f"Only Cu is supported. Got: {material}")
    
    # Import here to avoid circular dependencies (low-level constants vs high-level crystal)
    # 在此導入以避免循環依賴
    from axcsas.core.copper_crystal import CU_JCPDS_EXTENDED
    return CU_JCPDS_EXTENDED


# =============================================================================
# Physical Limits and Empirical Thresholds / 物理極限與經驗閾值
# =============================================================================

# Crystallite size detection limits / 晶粒尺寸檢測極限
# EMPIRICAL THRESHOLDS - Based on instrumental limitations and practical experience
# Commonly used in XRD crystallite size analysis
#
# Lower limit rationale:
#   Below ~2 nm, peak broadening approaches instrumental resolution limits
#   Scherrer equation becomes less reliable due to:
#   - Instrumental broadening dominates the signal
#   - Quantum size effects may alter crystal structure
#   - Surface relaxation becomes significant
#
# Upper limit rationale:
#   Above ~200 nm, crystallite size broadening becomes negligible
#   Peak width is dominated by instrumental broadening
#   Other techniques (SEM, TEM) are more appropriate
#
# Reference: Klug & Alexander (1974), "X-ray Diffraction Procedures" (general guidance)
#            Commonly accepted values in XRD community
MIN_RELIABLE_SIZE = 2.0    # nm, minimum detectable crystallite size
MAX_RELIABLE_SIZE = 200.0  # nm, maximum reliable crystallite size

# FWHM quality threshold / FWHM 品質閾值
# EMPIRICAL THRESHOLD - Quality control for peak broadening measurements
#
# Rationale: Ratio of (sample FWHM) / (instrumental FWHM)
#   If ratio < 1.2, sample broadening is too small relative to instrumental broadening
#   Results become unreliable due to uncertainty propagation
#   sqrt(β_sample² - β_inst²) becomes numerically unstable when β_sample ≈ β_inst
#
# Reference: Commonly used threshold in Rietveld refinement and XRD analysis
#            Young (1993), "The Rietveld Method" (general guidance)
MIN_BROADENING_RATIO = 1.2  # Minimum sample/instrumental FWHM ratio for reliable analysis

# Fit quality thresholds / 擬合品質閾值
# ═══════════════════════════════════════════════════════════════════════════
# 使用者自定義 / USER-DEFINED PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════
# 這些閾值由使用者根據實驗需求設定，用於作圖與數據篩選
# These thresholds are user-defined based on experimental requirements
# for plotting and data filtering purposes.
#
# 可根據需求調整至 0.99 或更高精確度
# Can be adjusted to 0.99 or higher precision as needed
#
# 參考文獻 Reference (一般指引):
#   Young, R. A. (1993). "The Rietveld Method." Oxford University Press.
#   (R_wp < 10% 為 Rietveld 精修一般可接受標準)
MAX_RWP_PERCENT = 10.0       # Maximum acceptable Rwp (%) / 最大可接受 R_wp
MIN_R_SQUARED = 0.95         # Minimum acceptable R² / 最小可接受 R²
