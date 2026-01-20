"""
AXCSAS Methods Module
=====================

Analysis methods for XRD data.
XRD 資料分析方法。
"""

from axcsas.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
    GrainShape,
    calculate_scherrer,
    calculate_crystallite_size,
    generate_scherrer_report,
)

from axcsas.methods.williamson_hall import (
    WilliamsonHallAnalyzer,
    WHResult,
    WHQualityLevel,
    analyze_williamson_hall,
    generate_wh_report,
    get_modulus_for_hkl,
)

from axcsas.methods.texture import (
    TextureAnalyzer,
    TextureAnalysisResult,
    TCResult,
    OrientationType,
    analyze_texture,
    generate_texture_report,
    get_standard_intensity,
    calculate_texture_coefficient,
)

from axcsas.methods.caglioti import (
    CagliotiCorrection,
    CagliotiParams,
    calculate_instrumental_broadening,
)

__all__ = [
    # Scherrer
    "ScherrerCalculator",
    "ScherrerResult",
    "ValidityFlag",
    "GrainShape",
    "calculate_scherrer",
    "calculate_crystallite_size",
    "generate_scherrer_report",
    # W-H
    "WilliamsonHallAnalyzer",
    "WHResult",
    "WHQualityLevel",
    "analyze_williamson_hall",
    "generate_wh_report",
    "get_modulus_for_hkl",
    # Texture
    "TextureAnalyzer",
    "TextureAnalysisResult",
    "TCResult",
    "OrientationType",
    "analyze_texture",
    "generate_texture_report",
    "get_standard_intensity",
    "calculate_texture_coefficient",
    # Caglioti
    "CagliotiCorrection",
    "CagliotiParams",
    "calculate_instrumental_broadening",
]

