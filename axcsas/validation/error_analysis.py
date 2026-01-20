"""
Error Analysis Module 誤差分析模組
===================================
Implements instrumental limits and reliability validation.
實現儀器限制和可靠性驗證。
"""

import numpy as np
from typing import Optional, List, NamedTuple
from dataclasses import dataclass
from enum import Enum


class WarningLevel(Enum):
    """Warning severity levels."""
    INFO = "info"
    WARNING = "warning"
    CRITICAL = "critical"


@dataclass
class ValidationWarning:
    """A validation warning message."""
    level: WarningLevel
    message: str
    parameter: str
    value: float
    threshold: float
    
    def __repr__(self):
        return f"[{self.level.value.upper()}] {self.message}"


@dataclass
class ValidationResult:
    """Result of validation check."""
    is_valid: bool
    warnings: List[ValidationWarning]
    
    def __repr__(self):
        status = "✓" if self.is_valid else "✗"
        return f"Validation({status}, {len(self.warnings)} warnings)"


class ErrorAnalyzer:
    """
    Analyzer for XRD analysis error and reliability.
    
    Checks include:
    - Crystallite size range validation (2-200 nm)
    - Instrumental broadening ratio (β_obs > 1.2 × β_inst)
    - Goodness of fit thresholds
    - Physical consistency
    """
    
    # Default thresholds
    MIN_SIZE_NM = 2.0          # Below: precision issues
    MAX_SIZE_NM = 200.0        # Above: exceeds detection limit
    MIN_BROADENING_RATIO = 1.2  # β_obs / β_inst threshold
    MAX_RWP = 10.0             # Maximum acceptable Rwp (%)
    MIN_R_SQUARED = 0.95       # Minimum acceptable R²
    
    def __init__(
        self,
        min_size: float = MIN_SIZE_NM,
        max_size: float = MAX_SIZE_NM,
        min_broadening_ratio: float = MIN_BROADENING_RATIO,
        max_rwp: float = MAX_RWP,
        min_r_squared: float = MIN_R_SQUARED
    ):
        """
        Initialize error analyzer with custom thresholds.
        
        Args:
            min_size: Minimum reliable size (nm)
            max_size: Maximum reliable size (nm)
            min_broadening_ratio: β_obs / β_inst threshold
            max_rwp: Maximum acceptable Rwp (%)
            min_r_squared: Minimum acceptable R²
        """
        self.min_size = min_size
        self.max_size = max_size
        self.min_broadening_ratio = min_broadening_ratio
        self.max_rwp = max_rwp
        self.min_r_squared = min_r_squared
    
    def validate_size(self, size_nm: float) -> ValidationResult:
        """
        Validate crystallite size is within reliable range.
        
        Args:
            size_nm: Crystallite size in nanometers
            
        Returns:
            ValidationResult
        """
        warnings = []
        is_valid = True
        
        if size_nm < self.min_size:
            warnings.append(ValidationWarning(
                level=WarningLevel.WARNING,
                message=f"Size {size_nm:.1f} nm below precision limit",
                parameter="crystallite_size",
                value=size_nm,
                threshold=self.min_size
            ))
            is_valid = False
        
        if size_nm > self.max_size:
            warnings.append(ValidationWarning(
                level=WarningLevel.CRITICAL,
                message=f"Size {size_nm:.1f} nm exceeds detection limit",
                parameter="crystallite_size",
                value=size_nm,
                threshold=self.max_size
            ))
            is_valid = False
        
        if np.isinf(size_nm) or np.isnan(size_nm):
            warnings.append(ValidationWarning(
                level=WarningLevel.CRITICAL,
                message="Invalid size value (inf or nan)",
                parameter="crystallite_size",
                value=size_nm,
                threshold=0
            ))
            is_valid = False
        
        return ValidationResult(is_valid=is_valid, warnings=warnings)
    
    def validate_broadening(
        self,
        fwhm_observed: float,
        fwhm_instrumental: float
    ) -> ValidationResult:
        """
        Validate broadening ratio for reliable correction.
        
        Args:
            fwhm_observed: Observed FWHM (any units)
            fwhm_instrumental: Instrumental FWHM (same units)
            
        Returns:
            ValidationResult
        """
        warnings = []
        is_valid = True
        
        if fwhm_instrumental <= 0:
            warnings.append(ValidationWarning(
                level=WarningLevel.CRITICAL,
                message="Instrumental FWHM must be positive",
                parameter="fwhm_instrumental",
                value=fwhm_instrumental,
                threshold=0
            ))
            return ValidationResult(is_valid=False, warnings=warnings)
        
        ratio = fwhm_observed / fwhm_instrumental
        
        if ratio < self.min_broadening_ratio:
            warnings.append(ValidationWarning(
                level=WarningLevel.WARNING,
                message=f"Broadening ratio {ratio:.2f} below threshold, "
                       f"instrumental broadening may dominate",
                parameter="broadening_ratio",
                value=ratio,
                threshold=self.min_broadening_ratio
            ))
            is_valid = False
        
        if fwhm_observed <= fwhm_instrumental:
            warnings.append(ValidationWarning(
                level=WarningLevel.CRITICAL,
                message="Observed FWHM smaller than instrumental, "
                       "correction impossible",
                parameter="fwhm_observed",
                value=fwhm_observed,
                threshold=fwhm_instrumental
            ))
            is_valid = False
        
        return ValidationResult(is_valid=is_valid, warnings=warnings)
    
    def validate_fit_quality(
        self,
        rwp: Optional[float] = None,
        r_squared: Optional[float] = None
    ) -> ValidationResult:
        """
        Validate fit quality metrics.
        
        Args:
            rwp: Weighted profile R-factor (%)
            r_squared: Coefficient of determination
            
        Returns:
            ValidationResult
        """
        warnings = []
        is_valid = True
        
        if rwp is not None and rwp > self.max_rwp:
            warnings.append(ValidationWarning(
                level=WarningLevel.WARNING,
                message=f"Rwp {rwp:.1f}% exceeds threshold",
                parameter="rwp",
                value=rwp,
                threshold=self.max_rwp
            ))
            is_valid = False
        
        if r_squared is not None and r_squared < self.min_r_squared:
            warnings.append(ValidationWarning(
                level=WarningLevel.WARNING,
                message=f"R² {r_squared:.4f} below threshold",
                parameter="r_squared",
                value=r_squared,
                threshold=self.min_r_squared
            ))
            is_valid = False
        
        return ValidationResult(is_valid=is_valid, warnings=warnings)
    
    def validate_all(
        self,
        size_nm: Optional[float] = None,
        fwhm_observed: Optional[float] = None,
        fwhm_instrumental: Optional[float] = None,
        rwp: Optional[float] = None,
        r_squared: Optional[float] = None
    ) -> ValidationResult:
        """
        Perform all applicable validations.
        
        Args:
            size_nm: Crystallite size
            fwhm_observed: Observed FWHM
            fwhm_instrumental: Instrumental FWHM
            rwp: Weighted R-factor
            r_squared: R² value
            
        Returns:
            Combined ValidationResult
        """
        all_warnings = []
        is_valid = True
        
        if size_nm is not None:
            result = self.validate_size(size_nm)
            all_warnings.extend(result.warnings)
            is_valid = is_valid and result.is_valid
        
        if fwhm_observed is not None and fwhm_instrumental is not None:
            result = self.validate_broadening(fwhm_observed, fwhm_instrumental)
            all_warnings.extend(result.warnings)
            is_valid = is_valid and result.is_valid
        
        if rwp is not None or r_squared is not None:
            result = self.validate_fit_quality(rwp, r_squared)
            all_warnings.extend(result.warnings)
            is_valid = is_valid and result.is_valid
        
        return ValidationResult(is_valid=is_valid, warnings=all_warnings)


# Convenience functions

def validate_size_range(
    size_nm: float,
    min_size: float = 2.0,
    max_size: float = 200.0
) -> bool:
    """
    Quick check if size is within reliable range.
    
    Args:
        size_nm: Crystallite size (nm)
        min_size: Minimum reliable size (nm)
        max_size: Maximum reliable size (nm)
        
    Returns:
        True if size is within range
    """
    return min_size <= size_nm <= max_size


def check_broadening_ratio(
    fwhm_observed: float,
    fwhm_instrumental: float,
    min_ratio: float = 1.2
) -> bool:
    """
    Quick check if broadening ratio is sufficient.
    
    Args:
        fwhm_observed: Observed FWHM
        fwhm_instrumental: Instrumental FWHM
        min_ratio: Minimum ratio threshold
        
    Returns:
        True if ratio is sufficient
    """
    if fwhm_instrumental <= 0:
        return False
    return (fwhm_observed / fwhm_instrumental) >= min_ratio
