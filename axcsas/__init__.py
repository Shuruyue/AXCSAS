"""
AXCSAS: Advanced XRD Crystallite Size Analysis System
=====================================================

A comprehensive toolkit for XRD data analysis including:
- Scherrer crystallite size calculation
- Williamson-Hall size/strain separation
- Harris Texture Coefficient analysis
- Pseudo-Voigt peak fitting

Example:
    >>> import axcsas
    >>> print(axcsas.__version__)
    0.1.0
"""

from axcsas.__version__ import __version__, __version_info__

from axcsas.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
    calculate_scherrer,
    calculate_crystallite_size,
)

from axcsas.methods.williamson_hall import (
    WilliamsonHallAnalyzer,
    WHResult,
    WHQualityLevel,
    analyze_williamson_hall,
)

from axcsas.methods.texture import (
    TextureAnalyzer,
    TextureAnalysisResult,
    analyze_texture,
)

__all__ = [
    # Version
    "__version__",
    "__version_info__",
    # Scherrer
    "ScherrerCalculator",
    "ScherrerResult",
    "ValidityFlag",
    "calculate_scherrer",
    "calculate_crystallite_size",
    # Williamson-Hall
    "WilliamsonHallAnalyzer",
    "WHResult",
    "WHQualityLevel",
    "analyze_williamson_hall",
    # Texture
    "TextureAnalyzer",
    "TextureAnalysisResult",
    "analyze_texture",
]
