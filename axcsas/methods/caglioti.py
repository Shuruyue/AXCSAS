"""
Caglioti Equation Module Caglioti方程模組
==========================================
Implements instrumental broadening correction using Caglioti equation.
使用 Caglioti 方程實現儀器展寬校正。

Reference 出處:
    Caglioti, G., Paoletti, A., & Ricci, F. P. (1958).
    Nucl. Instr., 3, 223-228.
"""

import numpy as np
from typing import Optional, Tuple
from dataclasses import dataclass


@dataclass
class CagliotiParams:
    """Caglioti equation parameters."""
    U: float  # deg²
    V: float  # deg²
    W: float  # deg²
    
    def __post_init__(self):
        if any(p is None for p in [self.U, self.V, self.W]):
            raise ValueError("All Caglioti parameters must be calibrated")


class CagliotiCorrection:
    """
    Instrumental broadening correction using Caglioti equation.
    
    FWHM²_inst = U·tan²θ + V·tanθ + W
    
    The Caglioti equation describes how instrumental broadening
    varies with diffraction angle 2θ.
    
    Reference: Caglioti, Paoletti, and Ricci (1958)
    """
    
    def __init__(
        self,
        U: Optional[float] = None,
        V: Optional[float] = None,
        W: Optional[float] = None
    ):
        """
        Initialize Caglioti correction.
        
        Args:
            U, V, W: Caglioti parameters (calibrated from standard)
        """
        self.U = U
        self.V = V
        self.W = W
        self._calibrated = all(p is not None for p in [U, V, W])
    
    @property
    def is_calibrated(self) -> bool:
        """Check if instrument has been calibrated."""
        return self._calibrated
    
    def calculate_fwhm_inst(self, two_theta: float) -> float:
        """
        Calculate instrumental FWHM at given 2θ.
        
        Args:
            two_theta: Diffraction angle in degrees
            
        Returns:
            Instrumental FWHM in degrees
        """
        if not self._calibrated:
            raise ValueError(
                "Instrument not calibrated. "
                "Run calibrate() with standard data first."
            )
        
        theta_rad = np.radians(two_theta / 2)
        tan_theta = np.tan(theta_rad)
        
        fwhm_sq = self.U * tan_theta**2 + self.V * tan_theta + self.W
        
        if fwhm_sq < 0:
            raise ValueError(f"Negative FWHM² at 2θ={two_theta}°. Check calibration.")
        
        return np.sqrt(fwhm_sq)
    
    def correct_broadening(
        self,
        fwhm_observed: float,
        two_theta: float,
        method: str = "geometric"
    ) -> Tuple[float, bool, str]:
        """
        Correct observed FWHM for instrumental broadening.
        
        For Pseudo-Voigt profiles, simple subtraction is incorrect.
        We use the geometric approximation:
        
        β_sample = β_obs - β²_inst / β_obs
        
        This is valid when β_obs > 1.2 × β_inst (error < 1%)
        
        Note 注意:
            "忽視儀器展寬會嚴重低估晶粒尺寸（特別是 D > 50 nm）"
        
        Args:
            fwhm_observed: Observed FWHM (degrees)
            two_theta: Diffraction angle (degrees)
            method: Correction method ("geometric" or "quadratic")
            
        Returns:
            Tuple of (corrected_fwhm, is_reliable, warning_message)
        """
        fwhm_inst = self.calculate_fwhm_inst(two_theta)
        
        # Check reliability threshold
        ratio = fwhm_observed / fwhm_inst
        is_reliable = ratio >= 1.2
        
        # Generate appropriate warning message
        if ratio < 1.0:
            warning = (
                f"CRITICAL: β_obs ({fwhm_observed:.4f}°) < β_inst ({fwhm_inst:.4f}°). "
                "Observed width smaller than instrument - check data or calibration."
            )
        elif ratio < 1.2:
            warning = (
                f"WARNING: β_obs/β_inst = {ratio:.2f} (< 1.2). "
                "Result may significantly underestimate crystallite size (D > 50 nm range). "
                "Consider this result unreliable."
            )
        elif ratio < 1.5:
            warning = (
                f"CAUTION: β_obs/β_inst = {ratio:.2f} (borderline). "
                "Result acceptable but approaching detection limit."
            )
        else:
            warning = ""
        
        if method == "geometric":
            # Geometric approximation for Pseudo-Voigt
            fwhm_sample = fwhm_observed - (fwhm_inst**2 / fwhm_observed)
        elif method == "quadratic":
            # Quadratic subtraction (valid for pure Gaussian)
            fwhm_sq_diff = fwhm_observed**2 - fwhm_inst**2
            if fwhm_sq_diff < 0:
                fwhm_sample = 0.001
            else:
                fwhm_sample = np.sqrt(fwhm_sq_diff)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Ensure positive result
        fwhm_sample = max(fwhm_sample, 0.001)
        
        return fwhm_sample, is_reliable, warning
    
    def calibrate(
        self,
        two_theta_std: np.ndarray,
        fwhm_std: np.ndarray
    ) -> CagliotiParams:
        """
        Calibrate Caglioti parameters from standard sample data.
        
        Args:
            two_theta_std: 2θ positions of standard peaks (degrees)
            fwhm_std: FWHM of standard peaks (degrees)
            
        Returns:
            CagliotiParams with fitted U, V, W values
        """
        # Convert to radians and calculate tan(θ)
        theta_rad = np.radians(two_theta_std / 2)
        tan_theta = np.tan(theta_rad)
        
        # Fit: FWHM² = U·tan²θ + V·tanθ + W
        # Design matrix: [tan²θ, tanθ, 1]
        A = np.column_stack([
            tan_theta**2,
            tan_theta,
            np.ones_like(tan_theta)
        ])
        
        fwhm_sq = fwhm_std**2
        
        # Least squares fit
        coeffs, residuals, rank, s = np.linalg.lstsq(A, fwhm_sq, rcond=None)
        
        self.U, self.V, self.W = coeffs
        self._calibrated = True
        
        return CagliotiParams(U=self.U, V=self.V, W=self.W)
    
    def get_params(self) -> CagliotiParams:
        """Return current Caglioti parameters."""
        return CagliotiParams(U=self.U, V=self.V, W=self.W)


def calculate_instrumental_broadening(
    two_theta: float,
    U: float,
    V: float,
    W: float
) -> float:
    """
    Convenience function to calculate instrumental FWHM.
    
    Args:
        two_theta: Diffraction angle (degrees)
        U, V, W: Caglioti parameters
        
    Returns:
        Instrumental FWHM (degrees)
    """
    correction = CagliotiCorrection(U, V, W)
    return correction.calculate_fwhm_inst(two_theta)
