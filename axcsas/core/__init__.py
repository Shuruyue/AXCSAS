"""
AXCSAS Core Module
==================

Core domain logic including constants, configuration, and crystal data.
核心領域邏輯，包括常數、配置和晶體資料。
"""

from axcsas.core.constants import (
    CU_KA1,
    CU_KA2,
    CU_KA_AVG,
    KA2_KA1_RATIO,
    SCHERRER_K,
    MIN_RELIABLE_SIZE,
    MAX_RELIABLE_SIZE,
    MIN_BROADENING_RATIO,
)

from axcsas.core.copper_crystal import (
    get_k_for_hkl,
    get_standard_peaks,
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
    # Constants - X-ray wavelengths
    "CU_KA1",
    "CU_KA2",
    "CU_KA_AVG",
    "KA2_KA1_RATIO",
    # Constants - Scherrer
    "SCHERRER_K",
    # Constants - Limits
    "MIN_RELIABLE_SIZE",
    "MAX_RELIABLE_SIZE",
    "MIN_BROADENING_RATIO",
    # Constants - JCPDS
    # CU_JCPDS removed (Use copper_crystal.CU_JCPDS_EXTENDED)
    # Copper crystal
    "get_k_for_hkl",
    "get_standard_peaks",
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
