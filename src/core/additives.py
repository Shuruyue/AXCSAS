"""
Electroplating Additives Classification Module
================================================

Classification and properties of acidic copper sulfate plating additives.

Reference: 計劃書/01_樣品製備與背景知識.md § 4

Three-component additive system:
1. Suppressor (抑制劑): PEG + Cl⁻
2. Accelerator (加速劑): SPS
3. Leveler (平整劑): JGB, IMEP, SH110
"""

from dataclasses import dataclass
from enum import Enum
from typing import Tuple, Optional


class AdditiveType(Enum):
    """Classification of copper electroplating additives."""
    SUPPRESSOR = "suppressor"    # 抑制劑
    ACCELERATOR = "accelerator"  # 加速劑
    LEVELER = "leveler"          # 平整劑


@dataclass
class Additive:
    """
    Electroplating additive with expected microstructural effects.
    
    Attributes:
        name: Common name or abbreviation
        type: Additive classification
        common_chemicals: Typical chemical compositions
        texture_tendency: Expected preferred orientation in XRD
        xrd_effect: Description of XRD peak characteristics
        grain_size_effect: Effect on crystallite size
    """
    name: str
    type: AdditiveType
    common_chemicals: Tuple[str, ...]
    texture_tendency: str
    xrd_effect: str
    grain_size_effect: str
    
    def __repr__(self) -> str:
        return f"Additive({self.name}, {self.type.value})"


# =============================================================================
# Standard Additive Definitions
# =============================================================================

SUPPRESSOR = Additive(
    name="Suppressor",
    type=AdditiveType.SUPPRESSOR,
    common_chemicals=("PEG", "PPG", "Cl⁻"),
    texture_tendency="(111)",
    xrd_effect="Significant peak broadening due to grain refinement",
    grain_size_effect="Fine grains, reduces crystallite size"
)

ACCELERATOR = Additive(
    name="Accelerator",
    type=AdditiveType.ACCELERATOR,
    common_chemicals=("SPS", "MPS", "MPSA"),
    texture_tendency="(111)",
    xrd_effect="Peak spacing may decrease (stacking faults), promotes nanotwins",
    grain_size_effect="Moderate grain size, promotes twin formation"
)

# Leveler variants
LEVELER_JGB = Additive(
    name="JGB (Janus Green B)",
    type=AdditiveType.LEVELER,
    common_chemicals=("Janus Green B",),
    texture_tendency="random",
    xrd_effect="Very broad peaks, high hardness",
    grain_size_effect="Extremely fine grains"
)

LEVELER_IMEP = Additive(
    name="IMEP",
    type=AdditiveType.LEVELER,
    common_chemicals=("IMEP",),
    texture_tendency="random",
    xrd_effect="Broad peaks, may suppress (111) texture",
    grain_size_effect="Fine grains"
)

LEVELER_SH110 = Additive(
    name="SH110",
    type=AdditiveType.LEVELER,
    common_chemicals=("SH110",),
    texture_tendency="random",
    xrd_effect="Broad peaks, promotes surface planarization",
    grain_size_effect="Fine grains with improved uniformity"
)


# =============================================================================
# Additive Effect Lookup
# =============================================================================

ALL_ADDITIVES = {
    "suppressor": SUPPRESSOR,
    "accelerator": ACCELERATOR,
    "jgb": LEVELER_JGB,
    "imep": LEVELER_IMEP,
    "sh110": LEVELER_SH110,
}


def get_additive(name: str) -> Optional[Additive]:
    """
    Get additive by name (case-insensitive).
    
    Args:
        name: Additive name or abbreviation
        
    Returns:
        Additive object or None if not found
    """
    return ALL_ADDITIVES.get(name.lower())


def predict_xrd_effects(additive_type: AdditiveType) -> str:
    """
    Predict XRD characteristics based on additive type.
    
    Args:
        additive_type: Type of additive used
        
    Returns:
        Description of expected XRD effects
    """
    effects = {
        AdditiveType.SUPPRESSOR: (
            "Expected XRD effects:\n"
            "- Peak broadening (smaller crystallites)\n"
            "- (111) texture enhancement\n"
            "- Possible lattice expansion (Cl incorporation)"
        ),
        AdditiveType.ACCELERATOR: (
            "Expected XRD effects:\n"
            "- Moderate peak width\n"
            "- (111) texture\n"
            "- Possible peak splitting or asymmetry (twins/stacking faults)\n"
            "- Lattice expansion from S incorporation"
        ),
        AdditiveType.LEVELER: (
            "Expected XRD effects:\n"
            "- Very broad peaks (nano-crystalline)\n"
            "- Random or suppressed texture\n"
            "- High background (amorphous content)\n"
            "- Significant lattice strain"
        ),
    }
    return effects.get(additive_type, "Unknown additive type")
