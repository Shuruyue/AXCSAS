"""
Enhanced Scherrer Calculator Module
====================================

Crystallite size calculation with dynamic K values, validity flags,
and complete unit conversion safety.

Reference: 計劃書/04_晶粒尺寸計算_上篇_Scherrer.md
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, List
from enum import Enum
import sys
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.copper_crystal import get_k_for_hkl, SCHERRER_CUBIC_K
from fitting.hkl_assignment import assign_hkl


# =============================================================================
# Constants
# =============================================================================

# Cu Kα1 wavelength in Ångströms
WAVELENGTH_CU_KA1 = 1.54056

# Validity thresholds (文件 04 §3.2)
FWHM_RATIO_THRESHOLD = 1.2  # β_obs must be > 1.2 × β_inst

# Size limits (文件 04 §7.1)
SIZE_MIN_NM = 2.0      # Minimum detectable size
SIZE_MAX_NM = 200.0    # Beyond this, results unreliable


# =============================================================================
# Validity Flag System (文件 04 §3.2-3.3)
# =============================================================================

class ValidityFlag(Enum):
    """
    Scherrer calculation validity flags.
    
    Reference: 文件 04 §3.2, §3.3
    """
    VALID = "VALID"              # Normal calculation, trustworthy
    UNRELIABLE = "UNRELIABLE"    # FWHM_obs ≤ 1.2 × FWHM_inst
    WARNING = "WARNING"          # Size > 200 nm (beyond method limit)
    ERROR = "ERROR"              # Calculation failed


# =============================================================================
# Enhanced Scherrer Result
# =============================================================================

@dataclass
class ScherrerResultEnhanced:
    """
    Enhanced result from Scherrer crystallite size calculation.
    
    Includes validity flags and complete metadata per 文件 04 §6.1.
    """
    # Size results
    size_nm: float
    size_angstrom: float
    
    # Peak information
    two_theta: float
    hkl: Optional[Tuple[int, int, int]]
    k_factor: float
    
    # Broadening values
    fwhm_observed: float     # Observed FWHM (degrees)
    fwhm_instrumental: float  # Instrumental FWHM (degrees)
    fwhm_sample: float       # Sample broadening (degrees)
    fwhm_sample_rad: float   # Sample broadening (radians) - CRITICAL
    
    # Validity
    validity_flag: ValidityFlag
    warning_message: str
    is_reliable: bool
    
    def __repr__(self) -> str:
        hkl_str = f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]})" if self.hkl else "(?)"
        return (
            f"ScherrerResult: {hkl_str} @ {self.two_theta:.2f}° → "
            f"D = {self.size_nm:.1f} nm [{self.validity_flag.value}]"
        )


# =============================================================================
# Enhanced Scherrer Calculator
# =============================================================================

class ScherrerCalculatorEnhanced:
    """
    Enhanced Scherrer calculator for electrodeposited copper.
    
    Features:
    - Dynamic K values based on hkl (rejects spherical K=0.89)
    - Validity flag system (VALID/UNRELIABLE/WARNING)
    - Explicit unit conversion steps with safety checks
    - Integration with Caglioti instrumental broadening
    
    Reference: 文件 04 (全章)
    
    Calculation Pipeline (文件 04 §3.1):
    1. Get FWHM_inst from Caglioti
    2. Check validity: FWHM_obs > 1.2 × FWHM_inst
    3. Calculate β_sample (geometric deconvolution)
    4. [CRITICAL] Convert degrees → radians
    5. Calculate cos(θ)
    6. Look up K value by hkl
    7. D = Kλ / (β × cos θ)
    8. Convert Å → nm
    9. Check size limits
    """
    
    def __init__(
        self,
        wavelength: float = WAVELENGTH_CU_KA1,
        use_cubic_habit: bool = True,
        caglioti_params: Optional[Tuple[float, float, float]] = None
    ):
        """
        Initialize enhanced Scherrer calculator.
        
        Args:
            wavelength: X-ray wavelength in Å (default: Cu Kα1)
            use_cubic_habit: Use ED-Cu specific K values (default: True)
            caglioti_params: (U, V, W) tuple for instrumental broadening
        """
        self.wavelength = wavelength
        self.use_cubic_habit = use_cubic_habit
        self.caglioti_uvw = caglioti_params
    
    def calculate(
        self,
        two_theta: float,
        fwhm_observed: float,
        fwhm_instrumental: Optional[float] = None,
        hkl: Optional[Tuple[int, int, int]] = None
    ) -> ScherrerResultEnhanced:
        """
        Calculate crystallite size using Scherrer equation.
        
        Args:
            two_theta: Peak position in degrees (2θ)
            fwhm_observed: Observed FWHM in degrees
            fwhm_instrumental: Instrumental FWHM in degrees (optional)
            hkl: Miller indices (auto-detected if not provided)
            
        Returns:
            ScherrerResultEnhanced with size and metadata
        """
        warnings = []
        
        # Step 1: Auto-detect hkl if not provided
        if hkl is None:
            hkl = assign_hkl(two_theta)
        
        # Step 2: Get K value for this hkl
        if hkl is not None:
            k_factor = get_k_for_hkl(
                hkl[0], hkl[1], hkl[2],
                use_cubic_habit=self.use_cubic_habit
            )
        else:
            k_factor = SCHERRER_CUBIC_K.K_SPHERICAL if not self.use_cubic_habit else 0.94
            warnings.append("Unknown hkl, using fallback K value")
        
        # Step 3: Handle instrumental broadening
        if fwhm_instrumental is None:
            if self.caglioti_uvw is not None:
                fwhm_instrumental = self._calculate_caglioti(two_theta)
            else:
                # No instrumental correction - use observed FWHM directly
                fwhm_instrumental = 0.0
                warnings.append("No instrumental correction applied")
        
        # Step 4: Check validity threshold (文件 04 §3.2)
        validity_flag = ValidityFlag.VALID
        is_reliable = True
        
        if fwhm_instrumental > 0:
            ratio = fwhm_observed / fwhm_instrumental
            if ratio <= FWHM_RATIO_THRESHOLD:
                validity_flag = ValidityFlag.UNRELIABLE
                is_reliable = False
                warnings.append(
                    f"FWHM ratio {ratio:.2f} ≤ {FWHM_RATIO_THRESHOLD} - "
                    "result may significantly underestimate size"
                )
        
        # Step 5: Calculate sample broadening (geometric deconvolution)
        # β_sample = FWHM_obs - FWHM_inst² / FWHM_obs
        if fwhm_instrumental > 0 and fwhm_observed > fwhm_instrumental:
            fwhm_sample = fwhm_observed - (fwhm_instrumental**2 / fwhm_observed)
        else:
            fwhm_sample = fwhm_observed
        
        # Ensure positive sample broadening
        fwhm_sample = max(fwhm_sample, 0.001)
        
        # ==================================================================
        # STEP 6: CRITICAL - Convert degrees to radians
        # WARNING: Using degrees directly will give results 57x too large!
        # Reference: 文件 04 §1.3
        # ==================================================================
        fwhm_sample_rad = fwhm_sample * np.pi / 180.0
        theta_rad = (two_theta / 2.0) * np.pi / 180.0
        
        # Step 7: Calculate cos(θ)
        cos_theta = np.cos(theta_rad)
        
        # Step 8: Scherrer equation: D = Kλ / (β × cos θ)
        size_angstrom = (k_factor * self.wavelength) / (fwhm_sample_rad * cos_theta)
        
        # Step 9: Convert Å to nm
        size_nm = size_angstrom / 10.0
        
        # Step 10: Check size limits (文件 04 §7.1)
        if size_nm > SIZE_MAX_NM:
            if validity_flag == ValidityFlag.VALID:
                validity_flag = ValidityFlag.WARNING
            warnings.append(
                f"Size {size_nm:.1f} nm > {SIZE_MAX_NM} nm - "
                "beyond Scherrer method applicability"
            )
        elif size_nm < SIZE_MIN_NM:
            warnings.append(
                f"Size {size_nm:.1f} nm < {SIZE_MIN_NM} nm - "
                "approaching XRD detection limit"
            )
        
        return ScherrerResultEnhanced(
            size_nm=size_nm,
            size_angstrom=size_angstrom,
            two_theta=two_theta,
            hkl=hkl,
            k_factor=k_factor,
            fwhm_observed=fwhm_observed,
            fwhm_instrumental=fwhm_instrumental,
            fwhm_sample=fwhm_sample,
            fwhm_sample_rad=fwhm_sample_rad,
            validity_flag=validity_flag,
            warning_message="; ".join(warnings) if warnings else "",
            is_reliable=is_reliable
        )
    
    def _calculate_caglioti(self, two_theta: float) -> float:
        """
        Calculate instrumental FWHM using Caglioti equation.
        
        FWHM²_inst = U·tan²θ + V·tanθ + W
        """
        if self.caglioti_uvw is None:
            return 0.0
        
        U, V, W = self.caglioti_uvw
        theta_rad = (two_theta / 2.0) * np.pi / 180.0
        tan_theta = np.tan(theta_rad)
        
        fwhm_sq = U * tan_theta**2 + V * tan_theta + W
        
        if fwhm_sq < 0:
            return 0.0
        
        return np.sqrt(fwhm_sq)
    
    def batch_calculate(
        self,
        peaks: List[Tuple[float, float]],
        fwhm_instrumental: Optional[float] = None
    ) -> List[ScherrerResultEnhanced]:
        """
        Calculate crystallite sizes for multiple peaks.
        
        Args:
            peaks: List of (two_theta, fwhm_observed) tuples
            fwhm_instrumental: Common instrumental FWHM (or None for per-peak calculation)
            
        Returns:
            List of ScherrerResultEnhanced objects
        """
        results = []
        for two_theta, fwhm_obs in peaks:
            result = self.calculate(
                two_theta=two_theta,
                fwhm_observed=fwhm_obs,
                fwhm_instrumental=fwhm_instrumental
            )
            results.append(result)
        return results
    
    def average_size(
        self,
        results: List[ScherrerResultEnhanced],
        include_unreliable: bool = False
    ) -> Tuple[float, float]:
        """
        Calculate average crystallite size from multiple peaks.
        
        Args:
            results: List of ScherrerResultEnhanced objects
            include_unreliable: Include UNRELIABLE results in average
            
        Returns:
            Tuple of (average_size_nm, std_dev_nm)
        """
        sizes = []
        for r in results:
            if r.is_reliable or include_unreliable:
                if r.validity_flag != ValidityFlag.ERROR:
                    sizes.append(r.size_nm)
        
        if not sizes:
            return 0.0, 0.0
        
        return float(np.mean(sizes)), float(np.std(sizes))


# =============================================================================
# Convenience Functions
# =============================================================================

def calculate_scherrer_enhanced(
    two_theta: float,
    fwhm_observed: float,
    fwhm_instrumental: float = 0.0,
    use_cubic_habit: bool = True
) -> ScherrerResultEnhanced:
    """
    Convenience function for enhanced Scherrer calculation.
    
    Example (文件 04 §4):
        >>> result = calculate_scherrer_enhanced(43.32, 0.25, 0.08)
        >>> print(f"D = {result.size_nm:.1f} nm")
        D = 49.0 nm
    """
    calc = ScherrerCalculatorEnhanced(use_cubic_habit=use_cubic_habit)
    return calc.calculate(two_theta, fwhm_observed, fwhm_instrumental)


def generate_scherrer_report(
    results: List[ScherrerResultEnhanced],
    sample_name: str = "Unknown"
) -> str:
    """
    Generate formatted Scherrer analysis report.
    
    Reference: 文件 04 §6.1
    """
    lines = [
        "=== Scherrer Crystallite Size Analysis ===",
        f"Sample: {sample_name}",
        "Analysis: Using Electroplated Cu-specific K values",
        "",
        "| Peak   | 2θ (°)  | FWHM_obs | FWHM_sample | K     | D (nm) | Flag        |",
        "|--------|---------|----------|-------------|-------|--------|-------------|",
    ]
    
    valid_sizes = []
    for r in results:
        hkl_str = f"({r.hkl[0]}{r.hkl[1]}{r.hkl[2]})" if r.hkl else "(?)"
        lines.append(
            f"| {hkl_str:6s} | {r.two_theta:7.3f} | {r.fwhm_observed:8.3f} | "
            f"{r.fwhm_sample:11.3f} | {r.k_factor:.3f} | {r.size_nm:6.1f} | "
            f"{r.validity_flag.value:11s} |"
        )
        if r.is_reliable:
            valid_sizes.append(r.size_nm)
    
    if valid_sizes:
        avg = np.mean(valid_sizes)
        std = np.std(valid_sizes)
        lines.append("")
        lines.append(f"Average D: {avg:.1f} ± {std:.1f} nm")
        
        if std / avg > 0.1:
            lines.append("Note: Variation indicates anisotropic growth")
    
    return "\n".join(lines)
