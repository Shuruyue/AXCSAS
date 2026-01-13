"""
Copper Crystal Physical Constants Module
=========================================

Comprehensive physical parameters for FCC copper crystallography,
specifically tailored for electrodeposited (ED) copper thin films.

References:
- JCPDS 04-0836 (Copper Standard)
- Langford & Wilson, J. Appl. Cryst. 11, 102-113 (1978)
- Electrodeposition of Copper, M. Schlesinger (2000)

Physical Constraints (from 計劃書/01):
1. Electrodeposited Cu exhibits columnar grain growth
2. Scherrer K must use cubic habit correction, NOT spherical assumption
3. Elastic anisotropy: E varies ~3x between <111> and <100>
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List
from math import gcd, sqrt


# =============================================================================
# FCC Copper Crystal Structure Constants (298K Standard)
# =============================================================================

@dataclass(frozen=True)
class CopperCrystal:
    """
    Copper FCC crystal structure constants at 298K.
    
    Reference: JCPDS 04-0836
    
    Attributes:
        space_group: Hermann-Mauguin symbol (Fm-3m)
        space_group_number: International Tables number (225)
        lattice_constant: Standard lattice parameter in Ångströms
        density: Density in g/cm³
        packing_factor: Atomic packing factor (FCC theoretical limit)
    """
    space_group: str = "Fm-3m"
    space_group_number: int = 225
    lattice_constant: float = 3.6150  # Å (298K standard)
    density: float = 8.92  # g/cm³
    packing_factor: float = 0.74  # APF = π√2/6

    def __repr__(self) -> str:
        return (
            f"CopperCrystal(a₀={self.lattice_constant} Å, "
            f"space_group={self.space_group}, ρ={self.density} g/cm³)"
        )


# Default instance
CU_CRYSTAL = CopperCrystal()


# =============================================================================
# JCPDS 04-0836 Extended Standard Diffraction Data
# =============================================================================

CU_JCPDS_EXTENDED: Dict[Tuple[int, int, int], Dict] = {
    (1, 1, 1): {
        "two_theta": 43.297,      # 2θ in degrees (Cu Kα1, λ=1.54056 Å)
        "d_spacing": 2.088,       # Å, d = a₀/√3
        "intensity": 100,         # Relative intensity (I/I₀)
        "multiplicity": 8,        # Number of equivalent planes {111}
        "description": "Strongest peak, most commonly observed"
    },
    (2, 0, 0): {
        "two_theta": 50.433,      # 2θ in degrees
        "d_spacing": 1.808,       # Å, d = a₀/2
        "intensity": 46,
        "multiplicity": 6,        # {200}
        "description": "Second strongest peak"
    },
    (2, 2, 0): {
        "two_theta": 74.130,      # 2θ in degrees
        "d_spacing": 1.278,       # Å, d = a₀/√8
        "intensity": 20,
        "multiplicity": 12,       # {220}
        "description": "Third major peak"
    },
    (3, 1, 1): {
        "two_theta": 89.931,      # 2θ in degrees
        "d_spacing": 1.090,       # Å
        "intensity": 17,
        "multiplicity": 24,       # {311}
        "description": "Fourth peak (often weak in ED-Cu)"
    },
    (2, 2, 2): {
        "two_theta": 95.139,      # 2θ in degrees
        "d_spacing": 1.044,       # Å, d = a₀/√12
        "intensity": 5,
        "multiplicity": 8,        # {222}
        "description": "Fifth peak (may be absent in textured films)"
    },
}


def is_fcc_allowed(h: int, k: int, l: int) -> bool:
    """
    Check if (hkl) reflection is allowed by FCC extinction rules.
    
    FCC Selection Rule:
    - Diffraction occurs only when h, k, l are ALL ODD or ALL EVEN
    - Mixed indices are forbidden (systematically absent)
    
    Args:
        h, k, l: Miller indices
        
    Returns:
        True if reflection is allowed, False if forbidden
        
    Examples:
        >>> is_fcc_allowed(1, 1, 1)  # All odd
        True
        >>> is_fcc_allowed(2, 0, 0)  # All even
        True
        >>> is_fcc_allowed(1, 0, 0)  # Mixed - forbidden
        False
    """
    parities = [x % 2 for x in (h, k, l)]
    return len(set(parities)) == 1  # All same parity


# =============================================================================
# Scherrer K Constants for Cubic Habit Grains
# =============================================================================

@dataclass(frozen=True)
class ScherrerCubicK:
    """
    Scherrer constant K values for cubic habit crystallites.
    
    CRITICAL: Electrodeposited copper forms COLUMNAR grains with cubic habit,
    NOT spherical grains. The K value varies with diffraction direction.
    
    Reference: 
    - Langford & Wilson, J. Appl. Cryst. 11, 102-113 (1978)
    - 計劃書 04 §2.2 電鍍銅動態 K 值表
    
    Physical Meaning:
    - K relates the measured FWHM to crystallite dimension
    - For cubic grains, the projected shape varies with viewing direction
    - (111): View along cube body diagonal → K = 1.155
    - (200): View along cube edge → K = 1.000
    - (220): View along face diagonal → K = 1.061
    - (311): Complex direction → K = 1.116
    
    Attributes:
        K_111: K value for (111) reflection = 1.155
        K_200: K value for (200) reflection = 1.000
        K_220: K value for (220) reflection = 1.061 (updated from 0.707)
        K_311: K value for (311) reflection = 1.116
        K_SPHERICAL: Traditional spherical assumption (reference only)
    """
    K_111: float = 1.155     # Body diagonal, 文件 04 §2.2
    K_200: float = 1.000     # Cube edge direction
    K_220: float = 1.061     # Face diagonal, 文件 04 §2.2 (更新自 0.707)
    K_311: float = 1.116     # Complex direction, 文件 04 §2.2
    K_222: float = 1.155     # Same as (111) parallel planes
    K_SPHERICAL: float = 0.89  # Traditional assumption (NOT for ED-Cu)
    K_CUBIC_GENERAL: float = 0.94  # Average for cubic grains


# Default instance
SCHERRER_CUBIC_K = ScherrerCubicK()


def get_k_for_hkl(
    h: int, k: int, l: int, 
    use_cubic_habit: bool = True,
    fallback_value: float = 0.89
) -> float:
    """
    Get appropriate Scherrer K value for given (hkl) direction.
    
    For electrodeposited copper with cubic habit grains, the Scherrer
    constant varies with crystallographic direction due to the 
    non-spherical grain shape.
    
    Args:
        h, k, l: Miller indices of the diffraction peak
        use_cubic_habit: If True, use direction-dependent K for cubic grains
                        If False, return traditional spherical K=0.89
        fallback_value: K value for unmapped directions
        
    Returns:
        Appropriate Scherrer K value (dimensionless)
        
    Physical Rationale (Langford & Wilson 1978):
        For a cube-shaped crystallite:
        - Viewing along <111>: See hexagonal projection → K = 2/√3 = 1.155
        - Viewing along <100>: See square projection → K = 1.000
        - Viewing along <110>: K = 1.061 (updated per 文件 04 §2.2)
        
    Examples:
        >>> get_k_for_hkl(1, 1, 1)  # Cubic habit
        1.155
        >>> get_k_for_hkl(2, 0, 0)  # Cubic habit
        1.0
        >>> get_k_for_hkl(1, 1, 1, use_cubic_habit=False)  # Spherical
        0.89
    """
    if not use_cubic_habit:
        return 0.89  # Traditional spherical assumption
    
    # Normalize to simplest form for matching
    g = gcd(gcd(abs(h), abs(k)), abs(l)) if l != 0 else gcd(abs(h), abs(k))
    if g == 0:
        return fallback_value
    h_n, k_n, l_n = abs(h) // g, abs(k) // g, abs(l) // g
    
    # Sort for canonical representation
    hkl_sorted = tuple(sorted([h_n, k_n, l_n]))
    
    # K value mapping for Scherrer size calculation
    # Reference: 文件 04 §2.2 電鍍銅動態 K 值表
    K_MAP = {
        (1, 1, 1): 1.155,   # Body diagonal
        (0, 0, 1): 1.000,   # (100), (200), (400), etc. - cube edge
        (0, 1, 1): 1.061,   # (110), (220), etc. - face diagonal (updated)
        (1, 1, 3): 1.116,   # (311) - updated from 0.89
        (1, 2, 2): 1.155,   # (222) parallel to (111)
    }
    
    return K_MAP.get(hkl_sorted, fallback_value)


# =============================================================================
# Elastic Anisotropy Parameters (for Williamson-Hall Analysis)
# =============================================================================

@dataclass(frozen=True)
class CopperElasticModuli:
    """
    Direction-dependent Young's modulus for copper single crystal.
    
    CRITICAL for Williamson-Hall Analysis:
    Copper is elastically anisotropic. The elastic modulus varies by
    nearly 3x between the softest <100> and hardest <111> directions.
    
    Using isotropic approximation (E~130 GPa) introduces significant
    errors in microstrain calculations for textured electrodeposited films.
    
    Reference: Copper single crystal elastic constants
              (Simmons & Wang, 1971)
              
    Zener Anisotropy Ratio: A = 2C₄₄/(C₁₁-C₁₂) ≈ 3.2 for Cu
    """
    E_111: float = 191.0    # GPa, hardest direction (close-packed)
    E_100: float = 66.0     # GPa, softest direction
    E_110: float = 130.0    # GPa, intermediate
    E_isotropic: float = 130.0  # GPa, polycrystalline average (Voigt-Reuss-Hill)


# Default instance
CU_ELASTIC = CopperElasticModuli()


def get_youngs_modulus(h: int, k: int, l: int) -> float:
    """
    Get direction-dependent Young's modulus for copper.
    
    Required for proper Williamson-Hall analysis of textured films.
    
    Args:
        h, k, l: Miller indices of the diffraction peak
        
    Returns:
        Young's modulus in GPa
        
    Note:
        For directions not explicitly mapped, returns isotropic average.
        More precise calculations require the full elastic tensor.
    """
    # Normalize
    g = gcd(gcd(abs(h), abs(k)), abs(l)) if l != 0 else gcd(abs(h), abs(k))
    if g == 0:
        return CU_ELASTIC.E_isotropic
    h_n, k_n, l_n = abs(h) // g, abs(k) // g, abs(l) // g
    hkl_sorted = tuple(sorted([h_n, k_n, l_n]))
    
    E_MAP = {
        (1, 1, 1): 191.0,   # <111> - hardest
        (0, 0, 1): 66.0,    # <100> - softest
        (0, 1, 1): 130.0,   # <110> - intermediate
        (1, 1, 3): 130.0,   # <311> - approximation
        (1, 2, 2): 191.0,   # <222> parallel to <111>
    }
    
    return E_MAP.get(hkl_sorted, CU_ELASTIC.E_isotropic)


# =============================================================================
# Electrodeposited Copper Lattice Anomaly Detection
# =============================================================================

# Standard lattice constant range for electrodeposited copper
ELECTROPLATED_A_STANDARD = 3.6150  # Å (reference)
ELECTROPLATED_A_MIN = 3.6150       # Å (pure, stress-free)
ELECTROPLATED_A_MAX = 3.6200       # Å (with impurity expansion)


@dataclass
class LatticeValidationResult:
    """Result of lattice constant validation."""
    measured_value: float
    is_normal: bool
    deviation_percent: float
    warning_level: str  # "none", "minor", "significant", "critical"
    explanation: str


def validate_lattice_constant(measured_a: float) -> LatticeValidationResult:
    """
    Validate measured lattice constant against expected range.
    
    Electrodeposited copper often exhibits lattice expansion due to:
    1. Impurity incorporation (S, Cl, C from additives)
    2. Vacancy accumulation
    3. Residual stress
    
    Args:
        measured_a: Measured lattice constant in Ångströms
        
    Returns:
        LatticeValidationResult with validation status and explanation
    """
    deviation = measured_a - ELECTROPLATED_A_STANDARD
    deviation_percent = (deviation / ELECTROPLATED_A_STANDARD) * 100
    
    if abs(deviation_percent) < 0.05:
        return LatticeValidationResult(
            measured_value=measured_a,
            is_normal=True,
            deviation_percent=deviation_percent,
            warning_level="none",
            explanation="Lattice constant within normal range for pure copper."
        )
    elif deviation_percent > 0 and deviation_percent < 0.12:
        return LatticeValidationResult(
            measured_value=measured_a,
            is_normal=True,
            deviation_percent=deviation_percent,
            warning_level="minor",
            explanation=(
                f"Slight lattice expansion (+{deviation_percent:.3f}%). "
                "Common in freshly deposited copper due to impurity incorporation "
                "(S, Cl from additives) or vacancy accumulation."
            )
        )
    elif deviation_percent >= 0.12 and deviation_percent < 0.3:
        return LatticeValidationResult(
            measured_value=measured_a,
            is_normal=False,
            deviation_percent=deviation_percent,
            warning_level="significant",
            explanation=(
                f"Significant lattice expansion (+{deviation_percent:.3f}%). "
                "Indicates high impurity content or strong residual stress. "
                "Sample may be in 'as-deposited' state before self-annealing. "
                "Consider recording storage time."
            )
        )
    elif deviation_percent < 0 and deviation_percent >= -0.3:
        return LatticeValidationResult(
            measured_value=measured_a,
            is_normal=False,
            deviation_percent=deviation_percent,
            warning_level="significant",
            explanation=(
                f"Unusual lattice contraction ({deviation_percent:.3f}%). "
                "May indicate tensile residual stress, measurement error, "
                "or sample displacement in XRD geometry. Verify instrument alignment."
            )
        )
    elif deviation_percent < -0.3:
        return LatticeValidationResult(
            measured_value=measured_a,
            is_normal=False,
            deviation_percent=deviation_percent,
            warning_level="critical",
            explanation=(
                f"Extreme lattice contraction ({deviation_percent:.3f}%). "
                "Possible causes: severe tensile stress, instrument miscalibration, "
                "sample misalignment, or measurement artifact. Re-examine setup."
            )
        )
    else:
        return LatticeValidationResult(
            measured_value=measured_a,
            is_normal=False,
            deviation_percent=deviation_percent,
            warning_level="critical",
            explanation=(
                f"Extreme lattice deviation (+{deviation_percent:.3f}%). "
                "Possible causes: severe contamination, phase impurity, "
                "or measurement artifact. Re-examine sample purity and "
                "instrument calibration."
            )
        )


def explain_lattice_deviation(
    measured_a: float,
    sample_age_hours: Optional[float] = None
) -> str:
    """
    Generate detailed explanation for observed lattice constant deviation.
    
    This function provides physical interpretation when the measured
    lattice constant differs from the standard value (3.6150 Å).
    
    Args:
        measured_a: Measured lattice constant in Ångströms
        sample_age_hours: Hours since electrodeposition (for self-annealing context)
        
    Returns:
        Human-readable explanation string
    """
    result = validate_lattice_constant(measured_a)
    
    explanation_parts = [
        f"Measured Lattice Constant: {measured_a:.4f} Å",
        f"Reference Value: {ELECTROPLATED_A_STANDARD:.4f} Å",
        f"Deviation: {result.deviation_percent:+.4f}%",
        "",
        "Physical Interpretation:",
        "-" * 40,
        result.explanation,
    ]
    
    # Add self-annealing context if sample age is provided
    if sample_age_hours is not None:
        explanation_parts.extend([
            "",
            "Self-Annealing Context:",
            "-" * 40,
        ])
        if sample_age_hours < 1:
            explanation_parts.append(
                "Sample is in 'as-deposited' state. Expect fine grains (50-100 nm) "
                "and high internal stress. Lattice expansion is typical."
            )
        elif sample_age_hours < 24:
            explanation_parts.append(
                "Early self-annealing phase. Grain growth is beginning. "
                "Lattice may still show expansion but typically decreasing."
            )
        else:
            explanation_parts.append(
                f"Sample aged {sample_age_hours:.1f} hours. Self-annealing likely "
                "significant. Grains may have grown to micron scale. "
                "Lattice should approach standard value."
            )
    
    # Add possible causes for expansion
    if result.deviation_percent > 0.05:
        explanation_parts.extend([
            "",
            "Possible Causes of Lattice Expansion:",
            "1. Sulfur incorporation from SPS accelerator",
            "2. Chloride incorporation from suppressor system",
            "3. Carbon/organic residue from leveler molecules",
            "4. High vacancy concentration from rapid deposition",
            "5. Residual compressive stress (apparent expansion)",
        ])
    
    return "\n".join(explanation_parts)


# =============================================================================
# Convenience Functions
# =============================================================================

def get_jcpds_peak(hkl: Tuple[int, int, int]) -> Optional[Dict]:
    """Get JCPDS data for a specific (hkl) reflection."""
    return CU_JCPDS_EXTENDED.get(hkl)


def get_all_peaks() -> List[Tuple[int, int, int]]:
    """Get list of all standard Cu peak indices."""
    return list(CU_JCPDS_EXTENDED.keys())


def calculate_d_spacing(h: int, k: int, l: int, a: float = 3.6150) -> float:
    """
    Calculate d-spacing for cubic crystal.
    
    d = a / √(h² + k² + l²)
    """
    return a / sqrt(h**2 + k**2 + l**2)
