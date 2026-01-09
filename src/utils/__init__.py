"""
Utilities Module
Common utilities for XRD analysis.
"""

from .constants import (
    # Wavelengths
    CU_KA1, CU_KA2, CU_KA_AVG, CU_KB,
    CO_KA1, MO_KA1, CR_KA1,
    KA2_KA1_RATIO,
    # Scherrer
    SCHERRER_K, ScherrerConstants,
    # JCPDS
    CU_JCPDS, AG_JCPDS, AU_JCPDS,
    get_jcpds_data,
    # Limits
    MIN_RELIABLE_SIZE, MAX_RELIABLE_SIZE,
    MIN_BROADENING_RATIO, MAX_RWP_PERCENT, MIN_R_SQUARED
)

from .unit_converter import (
    # Angle conversions
    deg_to_rad, rad_to_deg,
    two_theta_to_theta, theta_to_two_theta,
    # Length conversions
    angstrom_to_nm, nm_to_angstrom,
    nm_to_meter, meter_to_nm,
    # FWHM conversions
    fwhm_deg_to_rad, fwhm_rad_to_deg,
    # Crystallographic
    d_spacing_to_two_theta, two_theta_to_d_spacing
)

__all__ = [
    # Wavelengths
    "CU_KA1", "CU_KA2", "CU_KA_AVG", "CU_KB",
    "CO_KA1", "MO_KA1", "CR_KA1",
    "KA2_KA1_RATIO",
    # Scherrer
    "SCHERRER_K", "ScherrerConstants",
    # JCPDS
    "CU_JCPDS", "AG_JCPDS", "AU_JCPDS",
    "get_jcpds_data",
    # Limits
    "MIN_RELIABLE_SIZE", "MAX_RELIABLE_SIZE",
    "MIN_BROADENING_RATIO", "MAX_RWP_PERCENT", "MIN_R_SQUARED",
    # Unit conversions
    "deg_to_rad", "rad_to_deg",
    "two_theta_to_theta", "theta_to_two_theta",
    "angstrom_to_nm", "nm_to_angstrom",
    "nm_to_meter", "meter_to_nm",
    "fwhm_deg_to_rad", "fwhm_rad_to_deg",
    "d_spacing_to_two_theta", "two_theta_to_d_spacing",
]
