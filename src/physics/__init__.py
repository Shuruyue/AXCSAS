# Physics Module
# Physics Calculation Core Module
"""
Physics Calculation Core Module
Contains physical equations for XRD crystallography analysis.
"""

from .caglioti import (
    CagliotiCorrection,
    CagliotiParams,
    calculate_instrumental_broadening
)
from .scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    GrainShape,
    calculate_crystallite_size
)
from .scherrer_enhanced import (
    ScherrerCalculatorEnhanced,
    ScherrerResultEnhanced,
    ValidityFlag,
    calculate_scherrer_enhanced,
    generate_scherrer_report,
    WAVELENGTH_CU_KA1,
    FWHM_RATIO_THRESHOLD,
)
from .williamson_hall import (
    WilliamsonHallAnalyzer,
    WHResult,
    williamson_hall_analysis
)
from .wh_enhanced import (
    WilliamsonHallEnhanced,
    WHResultEnhanced,
    WHQualityLevel,
    analyze_williamson_hall_enhanced,
    generate_wh_report,
    MODULUS_MAP,
    WH_K_FACTOR,
    R2_EXCELLENT,
    R2_ACCEPTABLE,
)
from .texture import (
    TextureAnalyzer,
    TextureResult,
    TextureAnalysisResult,
    calculate_texture_coefficient,
    CU_JCPDS_STANDARD
)
from .texture_enhanced import (
    TextureAnalyzerEnhanced,
    TextureResultEnhanced,
    TCResult,
    OrientationType,
    analyze_texture_enhanced,
    generate_texture_report,
    JCPDS_STANDARD_INTENSITY,
)
from .defect_analysis import (
    StackingFaultAnalyzer,
    StackingFaultResult,
    StackingFaultSeverity,
    LatticeMonitor,
    LatticeConstantResult,
    LatticeStatus,
    ResidualStressResult,
    StressType,
    AnnealingState,
    DefectAnalysisResult,
    analyze_stacking_faults,
    analyze_lattice,
    determine_annealing_state,
    calculate_lattice_constant,
    STANDARD_PEAK_SEPARATION,
    WARREN_G_COEFFICIENT,
    STANDARD_LATTICE_CONSTANT,
)
from .report_generator import (
    ComprehensiveResult,
    generate_comprehensive_report,
    generate_csv_summary,
    generate_process_recommendations,
)

__all__ = [
    # Caglioti
    "CagliotiCorrection",
    "CagliotiParams",
    "calculate_instrumental_broadening",
    # Scherrer (original)
    "ScherrerCalculator",
    "ScherrerResult",
    "GrainShape",
    "calculate_crystallite_size",
    # Scherrer (enhanced)
    "ScherrerCalculatorEnhanced",
    "ScherrerResultEnhanced",
    "ValidityFlag",
    "calculate_scherrer_enhanced",
    "generate_scherrer_report",
    "WAVELENGTH_CU_KA1",
    "FWHM_RATIO_THRESHOLD",
    # Williamson-Hall (original)
    "WilliamsonHallAnalyzer",
    "WHResult",
    "williamson_hall_analysis",
    # Williamson-Hall (enhanced)
    "WilliamsonHallEnhanced",
    "WHResultEnhanced",
    "WHQualityLevel",
    "analyze_williamson_hall_enhanced",
    "generate_wh_report",
    "MODULUS_MAP",
    "WH_K_FACTOR",
    "R2_EXCELLENT",
    "R2_ACCEPTABLE",
    # Texture (original)
    "TextureAnalyzer",
    "TextureResult",
    "TextureAnalysisResult",
    "calculate_texture_coefficient",
    "CU_JCPDS_STANDARD",
    # Texture (enhanced)
    "TextureAnalyzerEnhanced",
    "TextureResultEnhanced",
    "TCResult",
    "OrientationType",
    "analyze_texture_enhanced",
    "generate_texture_report",
    "JCPDS_STANDARD_INTENSITY",
    # Defect Analysis
    "StackingFaultAnalyzer",
    "StackingFaultResult",
    "StackingFaultSeverity",
    "LatticeMonitor",
    "LatticeConstantResult",
    "LatticeStatus",
    "ResidualStressResult",
    "StressType",
    "AnnealingState",
    "DefectAnalysisResult",
    "analyze_stacking_faults",
    "analyze_lattice",
    "determine_annealing_state",
    "calculate_lattice_constant",
    "STANDARD_PEAK_SEPARATION",
    "WARREN_G_COEFFICIENT",
    "STANDARD_LATTICE_CONSTANT",
    # Report Generator
    "ComprehensiveResult",
    "generate_comprehensive_report",
    "generate_csv_summary",
    "generate_process_recommendations",
]
