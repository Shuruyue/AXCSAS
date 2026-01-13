"""
Physical Constants Module
Centralized physical constants for XRD analysis.
"""

from dataclasses import dataclass
from typing import Dict, Tuple


# =============================================================================
# X-ray Wavelengths (Ångströms)
# =============================================================================

# Copper X-ray tube
CU_KA1 = 1.54056      # Cu Kα1 (IUPAC Standard)
CU_KA2 = 1.54439      # Cu Kα2
CU_KA_AVG = 1.54184   # Weighted average Kα
CU_KB = 1.39222       # Cu Kβ

# Cobalt X-ray tube
CO_KA1 = 1.78897      # Co Kα1
CO_KA2 = 1.79285      # Co Kα2

# Molybdenum X-ray tube
MO_KA1 = 0.70930      # Mo Kα1
MO_KA2 = 0.71359      # Mo Kα2

# Chromium X-ray tube
CR_KA1 = 2.28970      # Cr Kα1
CR_KA2 = 2.29361      # Cr Kα2

# Kα2/Kα1 intensity ratio (for Rachinger correction)
KA2_KA1_RATIO = 0.5


# =============================================================================
# Scherrer Constants (Langford & Wilson, 1978)
# =============================================================================

@dataclass
class ScherrerConstants:
    """Scherrer constant K values for different grain shapes."""
    spherical: float = 0.89
    cubic: float = 0.94
    tetrahedral: float = 0.94
    octahedral: float = 0.83
    disk_100: float = 0.89   # Disk with 1:100 aspect ratio
    disk_10: float = 0.84    # Disk with 1:10 aspect ratio
    columnar: float = 0.91   # Columnar grains
    default: float = 0.89    # Default (spherical assumption)


SCHERRER_K = ScherrerConstants()


# =============================================================================
# JCPDS Standard Data
# =============================================================================

# Copper (PDF 04-0836)
CU_JCPDS = {
    (1, 1, 1): {"two_theta": 43.298, "intensity": 100.0, "d_spacing": 2.0872},
    (2, 0, 0): {"two_theta": 50.434, "intensity": 46.0, "d_spacing": 1.8079},
    (2, 2, 0): {"two_theta": 74.130, "intensity": 20.0, "d_spacing": 1.2780},
    (3, 1, 1): {"two_theta": 89.931, "intensity": 17.0, "d_spacing": 1.0899},
    (2, 2, 2): {"two_theta": 95.139, "intensity": 5.0, "d_spacing": 1.0436},
}

# Silver (PDF 04-0783)
AG_JCPDS = {
    (1, 1, 1): {"two_theta": 38.116, "intensity": 100.0, "d_spacing": 2.3588},
    (2, 0, 0): {"two_theta": 44.277, "intensity": 40.0, "d_spacing": 2.0442},
    (2, 2, 0): {"two_theta": 64.426, "intensity": 25.0, "d_spacing": 1.4450},
    (3, 1, 1): {"two_theta": 77.472, "intensity": 26.0, "d_spacing": 1.2313},
    (2, 2, 2): {"two_theta": 81.536, "intensity": 12.0, "d_spacing": 1.1794},
}

# Gold (PDF 04-0784)
AU_JCPDS = {
    (1, 1, 1): {"two_theta": 38.184, "intensity": 100.0, "d_spacing": 2.3544},
    (2, 0, 0): {"two_theta": 44.392, "intensity": 52.0, "d_spacing": 2.0390},
    (2, 2, 0): {"two_theta": 64.576, "intensity": 32.0, "d_spacing": 1.4420},
    (3, 1, 1): {"two_theta": 77.547, "intensity": 36.0, "d_spacing": 1.2297},
    (2, 2, 2): {"two_theta": 81.721, "intensity": 12.0, "d_spacing": 1.1772},
}


def get_jcpds_data(material: str) -> Dict[Tuple[int, int, int], Dict]:
    """
    Get JCPDS standard data for a material.
    
    Args:
        material: Material name (Cu, Ag, Au)
        
    Returns:
        Dictionary of JCPDS data
    """
    material = material.upper()
    data_map = {
        "CU": CU_JCPDS,
        "AG": AG_JCPDS,
        "AU": AU_JCPDS,
    }
    
    if material not in data_map:
        raise ValueError(f"Unknown material: {material}. Available: {list(data_map.keys())}")
    
    return data_map[material]


# =============================================================================
# Physical Limits
# =============================================================================

# Crystallite size detection limits (nm)
MIN_RELIABLE_SIZE = 2.0      # Below: precision issues
MAX_RELIABLE_SIZE = 200.0    # Above: exceeds XRD detection limit

# Instrumental broadening correction threshold
MIN_BROADENING_RATIO = 1.2   # β_obs / β_inst

# Fit quality thresholds
MAX_RWP_PERCENT = 10.0       # Maximum acceptable Rwp (%)
MIN_R_SQUARED = 0.95         # Minimum acceptable R²
