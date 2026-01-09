"""
XRD Preprocessing Pipeline
==========================

Orchestrates the complete preprocessing workflow for XRD data:
1. Data loading and validation
2. Savitzky-Golay smoothing
3. Background subtraction
4. Kα2 stripping (conditional)

Reference: 計劃書/02_XRD數據採集與預處理.md §6
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple
from pathlib import Path
import time

from .data_loader import XRDDataLoader, load_xrd_data
from .validation import (
    XRDDataset, 
    DataValidationResult, 
    ValidationWarning,
    WarningLevel,
    validate_xrd_data,
    check_negative_values
)
from .smoothing import SavitzkyGolayFilter
from .background import BackgroundSubtractor
from .kalpha_strip import KalphaStripper


@dataclass
class PreprocessingStep:
    """Record of a single preprocessing step."""
    name: str
    duration_ms: float
    applied: bool
    notes: str = ""


@dataclass
class PreprocessingResult:
    """
    Result of the complete preprocessing pipeline.
    
    Attributes:
        two_theta: 2θ angle array
        intensity: Preprocessed intensity array
        raw_intensity: Original intensity (before processing)
        background: Estimated background (if computed)
        validation: Data validation result
        steps: List of applied preprocessing steps
        warnings: Accumulated warnings from all steps
    """
    two_theta: np.ndarray
    intensity: np.ndarray
    raw_intensity: np.ndarray
    background: Optional[np.ndarray] = None
    validation: Optional[DataValidationResult] = None
    steps: List[PreprocessingStep] = field(default_factory=list)
    warnings: List[ValidationWarning] = field(default_factory=list)
    
    def summary(self) -> str:
        """Generate human-readable processing summary."""
        lines = [
            "=== Preprocessing Summary ===",
            f"Data points: {len(self.two_theta)}",
            f"2θ range: {self.two_theta.min():.2f}° - {self.two_theta.max():.2f}°",
            "",
            "Steps performed:"
        ]
        for step in self.steps:
            status = "✓" if step.applied else "○"
            lines.append(f"  {status} {step.name} ({step.duration_ms:.1f} ms)")
            if step.notes:
                lines.append(f"      {step.notes}")
        
        if self.warnings:
            lines.append("")
            lines.append(f"Warnings ({len(self.warnings)}):")
            for w in self.warnings:
                lines.append(f"  [{w.level.value}] {w.message}")
        
        return "\n".join(lines)


class PreprocessingPipeline:
    """
    XRD data preprocessing pipeline.
    
    Orchestrates the complete preprocessing workflow with configurable
    parameters and optional steps.
    
    Processing Order (per 文件 02 §6):
    1. Data validation
    2. Savitzky-Golay smoothing
    3. Background subtraction
    4. Kα2 stripping (if 2θ_max > 40°)
    
    Example:
        >>> pipeline = PreprocessingPipeline()
        >>> result = pipeline.run(two_theta, intensity)
        >>> print(result.summary())
    """
    
    # Default parameters from 文件 02 §7
    DEFAULT_WINDOW_SIZE = 11
    DEFAULT_POLY_ORDER = 3
    DEFAULT_BG_METHOD = "chebyshev"
    DEFAULT_BG_DEGREE = 5
    KA2_THRESHOLD_ANGLE = 40.0  # 2θ below which Kα2 stripping is skipped
    
    def __init__(
        self,
        window_size: int = DEFAULT_WINDOW_SIZE,
        poly_order: int = DEFAULT_POLY_ORDER,
        background_method: str = DEFAULT_BG_METHOD,
        background_degree: int = DEFAULT_BG_DEGREE,
        enable_smoothing: bool = True,
        enable_background: bool = True,
        enable_kalpha_strip: bool = True,
        auto_correct_negative: bool = True
    ):
        """
        Initialize preprocessing pipeline.
        
        Args:
            window_size: Savitzky-Golay window size (must be odd)
            poly_order: Savitzky-Golay polynomial order
            background_method: "chebyshev" or "sonneveld_visser"
            background_degree: Polynomial degree for Chebyshev method
            enable_smoothing: Whether to apply smoothing
            enable_background: Whether to subtract background
            enable_kalpha_strip: Whether to strip Kα2
            auto_correct_negative: Automatically set negative values to zero
        """
        self.window_size = window_size
        self.poly_order = poly_order
        self.background_method = background_method
        self.background_degree = background_degree
        self.enable_smoothing = enable_smoothing
        self.enable_background = enable_background
        self.enable_kalpha_strip = enable_kalpha_strip
        self.auto_correct_negative = auto_correct_negative
        
        # Initialize processors
        self._init_processors()
    
    def _init_processors(self):
        """Initialize processing components."""
        self.smoother = SavitzkyGolayFilter(
            window_size=self.window_size,
            poly_order=self.poly_order
        )
        self.background_subtractor = BackgroundSubtractor(
            method=self.background_method,
            poly_degree=self.background_degree
        )
        self.kalpha_stripper = KalphaStripper()
    
    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> "PreprocessingPipeline":
        """
        Create pipeline from configuration dictionary.
        
        Expected config keys (matching config.yaml):
            preprocessing.smoothing.window_size
            preprocessing.smoothing.poly_order
            preprocessing.background.method
            preprocessing.background.poly_degree
            preprocessing.kalpha_strip.enable
        """
        smoothing_cfg = config.get('preprocessing', {}).get('smoothing', {})
        bg_cfg = config.get('preprocessing', {}).get('background', {})
        ka_cfg = config.get('preprocessing', {}).get('kalpha_strip', {})
        
        return cls(
            window_size=smoothing_cfg.get('window_size', cls.DEFAULT_WINDOW_SIZE),
            poly_order=smoothing_cfg.get('poly_order', cls.DEFAULT_POLY_ORDER),
            background_method=bg_cfg.get('method', cls.DEFAULT_BG_METHOD),
            background_degree=bg_cfg.get('poly_degree', cls.DEFAULT_BG_DEGREE),
            enable_smoothing=smoothing_cfg.get('enable', True),
            enable_background=bg_cfg.get('enable', True),
            enable_kalpha_strip=ka_cfg.get('enable', True)
        )
    
    def run(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray
    ) -> PreprocessingResult:
        """
        Run the complete preprocessing pipeline.
        
        Args:
            two_theta: 2θ angle array in degrees
            intensity: Raw intensity array
            
        Returns:
            PreprocessingResult with processed data and metadata
        """
        steps = []
        warnings = []
        raw_intensity = intensity.copy()
        current_intensity = intensity.copy()
        background = None
        
        # Step 1: Validation
        t0 = time.perf_counter()
        validation = validate_xrd_data(two_theta, intensity)
        t1 = time.perf_counter()
        
        steps.append(PreprocessingStep(
            name="Data Validation",
            duration_ms=(t1 - t0) * 1000,
            applied=True,
            notes=f"{'PASS' if validation.is_valid else 'FAIL'} - {len(validation.warnings)} warnings"
        ))
        warnings.extend(validation.warnings)
        
        # Step 2: Smoothing
        if self.enable_smoothing:
            t0 = time.perf_counter()
            try:
                current_intensity = self.smoother.apply(current_intensity)
                t1 = time.perf_counter()
                steps.append(PreprocessingStep(
                    name="Savitzky-Golay Smoothing",
                    duration_ms=(t1 - t0) * 1000,
                    applied=True,
                    notes=f"window={self.window_size}, order={self.poly_order}"
                ))
            except ValueError as e:
                t1 = time.perf_counter()
                steps.append(PreprocessingStep(
                    name="Savitzky-Golay Smoothing",
                    duration_ms=(t1 - t0) * 1000,
                    applied=False,
                    notes=f"Skipped: {e}"
                ))
                warnings.append(ValidationWarning(
                    level=WarningLevel.WARNING,
                    code="SMOOTHING_SKIPPED",
                    message=str(e)
                ))
        
        # Step 3: Background subtraction
        if self.enable_background:
            t0 = time.perf_counter()
            current_intensity, background = self.background_subtractor.subtract(
                two_theta, current_intensity
            )
            t1 = time.perf_counter()
            steps.append(PreprocessingStep(
                name="Background Subtraction",
                duration_ms=(t1 - t0) * 1000,
                applied=True,
                notes=f"method={self.background_method}"
            ))
            
            # Check for negative values after background subtraction
            n_negative = np.sum(current_intensity < 0)
            if n_negative > 0:
                warnings.append(ValidationWarning(
                    level=WarningLevel.WARNING,
                    code="POST_BG_NEGATIVE",
                    message=f"{n_negative} negative values after background subtraction"
                ))
                if self.auto_correct_negative:
                    current_intensity, _ = check_negative_values(
                        current_intensity, auto_correct=True
                    )
        
        # Step 4: Kα2 stripping (conditional)
        if self.enable_kalpha_strip:
            should_strip = should_apply_kalpha_stripping(two_theta)
            
            if should_strip:
                t0 = time.perf_counter()
                current_intensity = self.kalpha_stripper.strip(
                    two_theta, current_intensity
                )
                t1 = time.perf_counter()
                steps.append(PreprocessingStep(
                    name="Kα2 Stripping",
                    duration_ms=(t1 - t0) * 1000,
                    applied=True,
                    notes="Rachinger correction applied"
                ))
            else:
                steps.append(PreprocessingStep(
                    name="Kα2 Stripping",
                    duration_ms=0,
                    applied=False,
                    notes=f"Skipped: 2θ_max < {self.KA2_THRESHOLD_ANGLE}°"
                ))
        
        return PreprocessingResult(
            two_theta=two_theta,
            intensity=current_intensity,
            raw_intensity=raw_intensity,
            background=background,
            validation=validation,
            steps=steps,
            warnings=warnings
        )
    
    def run_from_file(self, filepath: str) -> PreprocessingResult:
        """
        Load and preprocess XRD data from file.
        
        Args:
            filepath: Path to XRD data file (.xy, .csv, .txt)
            
        Returns:
            PreprocessingResult
        """
        two_theta, intensity = load_xrd_data(filepath)
        return self.run(two_theta, intensity)


def should_apply_kalpha_stripping(
    two_theta: np.ndarray,
    threshold: float = 40.0
) -> bool:
    """
    Determine if Kα2 stripping should be applied.
    
    Physical Rationale (引用文件 02 §4.1):
        At low angles (2θ < 40°), the Kα1/Kα2 doublet is not resolved
        and stripping has minimal effect. At high angles (2θ > 60°),
        the doublet becomes clearly separated and must be corrected.
    
    Args:
        two_theta: 2θ angle array
        threshold: Minimum 2θ to consider stripping (default: 40°)
        
    Returns:
        True if stripping should be applied
    """
    return float(two_theta.max()) > threshold


def get_kalpha_shift_table() -> Dict[int, float]:
    """
    Return reference table of Kα1-Kα2 angular shifts.
    
    Reference: 文件 02 §4.3
    
    Returns:
        Dictionary mapping 2θ (degrees) to Δ2θ (degrees)
    """
    return {
        40: 0.05,
        60: 0.10,
        80: 0.20,
        100: 0.40,
    }
