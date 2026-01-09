"""
Core Module Initialization
Central module for AXCSAS physical constants and crystal parameters.
"""

from .copper_crystal import (
    CopperCrystal,
    CU_JCPDS_EXTENDED,
    ScherrerCubicK,
    CopperElasticModuli,
    get_k_for_hkl,
    get_youngs_modulus,
    is_fcc_allowed,
    validate_lattice_constant,
    explain_lattice_deviation,
)

from .additives import (
    AdditiveType,
    Additive,
    SUPPRESSOR,
    ACCELERATOR,
    LEVELER_JGB,
    LEVELER_IMEP,
    LEVELER_SH110,
)

__all__ = [
    # Copper crystal
    "CopperCrystal",
    "CU_JCPDS_EXTENDED",
    "ScherrerCubicK",
    "CopperElasticModuli",
    "get_k_for_hkl",
    "get_youngs_modulus",
    "is_fcc_allowed",
    "validate_lattice_constant",
    "explain_lattice_deviation",
    # Additives
    "AdditiveType",
    "Additive",
    "SUPPRESSOR",
    "ACCELERATOR",
    "LEVELER_JGB",
    "LEVELER_IMEP",
    "LEVELER_SH110",
]
