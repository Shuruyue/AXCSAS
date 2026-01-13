"""
AXCSAS Core Module
==================

Core domain logic including constants, configuration, and crystal data.
核心領域邏輯，包括常數、配置和晶體資料。
"""

from axcsas.core.constants import (
    CU_KA1,
    CU_KA2,
    SCHERRER_K,
    MIN_RELIABLE_SIZE,
    MAX_RELIABLE_SIZE,
    MIN_BROADENING_RATIO,
    CU_JCPDS,
)

from axcsas.core.copper_crystal import (
    get_k_for_hkl,
    SCHERRER_CUBIC_K,
    CU_ELASTIC,
)

from axcsas.core.config_loader import load_config

from axcsas.core.units import (
    deg_to_rad,
    rad_to_deg,
    angstrom_to_nm,
    nm_to_angstrom,
)

__all__ = [
    # Constants
    "CU_KA1",
    "CU_KA2",
    "SCHERRER_K",
    "MIN_RELIABLE_SIZE",
    "MAX_RELIABLE_SIZE",
    "MIN_BROADENING_RATIO",
    "CU_JCPDS",
    # Copper crystal
    "get_k_for_hkl",
    "SCHERRER_CUBIC_K",
    "CU_ELASTIC",
    # Config
    "load_config",
    # Units
    "deg_to_rad",
    "rad_to_deg",
    "angstrom_to_nm",
    "nm_to_angstrom",
]
