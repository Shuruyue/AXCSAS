"""
AXCSAS Complete Analysis Pipeline 完整分析管道
==========================================

Unified pipeline integrating all analysis phases.
統一的分析流程，整合所有分析階段。

- Phase 02: Preprocessing 預處理
- Phase 03: Peak Fitting 峰値擬合
- Phase 04: Scherrer Size Calculation 晶粒尺寸計算
- Phase 05: Williamson-Hall Strain Analysis 應變分析
- Phase 06: Texture Analysis 織構分析
- Phase 07: Defect and Stress Diagnosis 缺陷與應力診斷
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict, Any
from pathlib import Path
import re

from axcsas.core.constants import CU_KA1
from axcsas.core.copper_crystal import CU_JCPDS_EXTENDED

from axcsas.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
)
from axcsas.methods.williamson_hall import (
    WilliamsonHallAnalyzer,
    WHResult,
)
from axcsas.methods.texture import (
    TextureAnalyzer,
    TextureAnalysisResult,
)
from axcsas.methods.defect_analysis import (
    StackingFaultAnalyzer,
    StackingFaultResult,
    LatticeMonitor,
    LatticeConstantResult,
    AnnealingState,
    determine_annealing_state,
)
from axcsas.analysis.report_generator import (
    ComprehensiveResult,
    generate_comprehensive_report,
    generate_csv_summary,
)
from axcsas.fitting.hkl_assignment import assign_hkl
from axcsas.fitting.lm_optimizer import LMOptimizer
from axcsas.fitting.pseudo_voigt import PseudoVoigt, PseudoVoigtParams


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class AnalysisConfig:
    """
    Configuration for AXCSAS analysis pipeline.
    AXCSAS 分析管道配置。
    """
    
    # X-ray parameters / X-ray 參數
    # Default: Cu Kα1 from constants (Bearden 1967)
    wavelength: float = CU_KA1
    
    # Scherrer parameters / Scherrer 參數
    use_cubic_habit: bool = True
    
    # Peak detection / 峰值偵測
    peak_window: float = 2.0  # degrees around expected position
    min_intensity: float = 100  # minimum counts
    
    # Instrumental broadening (Caglioti U, V, W) / 儀器展寬
    # Typical empirical values for lab diffractometers
    # FWHM_inst² = U·tan²θ + V·tanθ + W
    caglioti_u: float = 0.0
    caglioti_v: float = 0.0
    caglioti_w: float = 0.003  # → FWHM_inst = √0.003 ≈ 0.055°
    
    # Cu peak positions from JCPDS 04-0836 / 銅峰位從 JCPDS 標準
    # Dynamically generated from CU_JCPDS constants
    EXPECTED_PEAKS: Dict[Tuple[int, int, int], float] = field(
        default_factory=lambda: {
            hkl: data["two_theta"] for hkl, data in CU_JCPDS_EXTENDED.items()
        }
    )


# =============================================================================
# Result Container
# =============================================================================

@dataclass
class PeakData:
    """Single peak data."""
    hkl: Tuple[int, int, int]
    two_theta: float
    intensity: float
    fwhm: float
    area: float = 0.0


@dataclass
class PipelineResult:
    """Complete pipeline analysis result."""
    
    # Input
    filepath: str
    sample_name: str
    
    # Sample metadata
    leveler_concentration: Optional[float] = None
    plating_time_hours: Optional[float] = None
    sample_age_hours: Optional[float] = None
    
    # Peak data
    peaks: List[PeakData] = field(default_factory=list)
    
    # Phase 04: Scherrer
    scherrer_results: List[ScherrerResult] = field(default_factory=list)
    average_size_nm: Optional[float] = None
    
    # Phase 05: W-H
    wh_result: Optional[WHResult] = None
    
    # Phase 06: Texture
    texture_result: Optional[TextureAnalysisResult] = None
    
    # Phase 07: Defects
    stacking_fault: Optional[StackingFaultResult] = None
    lattice_result: Optional[LatticeConstantResult] = None
    annealing_state: AnnealingState = AnnealingState.UNKNOWN
    
    # Comprehensive
    comprehensive: Optional[ComprehensiveResult] = None
    report: str = ""


# =============================================================================
# Data Loader
# =============================================================================

def load_bruker_txt(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load Bruker TXT format XRD data.
    
    Returns:
        Tuple of (two_theta, intensity) arrays
    """
    two_theta = []
    intensity = []
    in_data_section = False
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            
            if '[Data]' in line:
                in_data_section = True
                continue
            
            if in_data_section and line:
                # Skip header row
                if 'Angle' in line or 'PSD' in line:
                    continue
                
                # Parse data
                parts = line.replace(',', ' ').split()
                if len(parts) >= 2:
                    try:
                        theta = float(parts[0])
                        counts = float(parts[1])
                        two_theta.append(theta)
                        intensity.append(counts)
                    except ValueError:
                        continue
    
    return np.array(two_theta), np.array(intensity)


def parse_filename(filepath: str) -> Dict[str, Any]:
    """
    Parse sample info from filename.
    
    Format: YYYYMMDD_Xml_Xh.txt or YYYYMMDD_Xml_Xh_Xmin.txt
    """
    name = Path(filepath).stem
    
    result = {
        'name': name,
        'concentration_ml': None,
        'time_hours': None,
    }
    
    # Extract concentration (e.g., "0ml", "4.5ml", "9ml", "18ml")
    conc_match = re.search(r'(\d+\.?\d*)ml', name)
    if conc_match:
        result['concentration_ml'] = float(conc_match.group(1))
    
    # Extract time (e.g., "2h", "0h_15min", "24h")
    time_hours = 0
    hour_match = re.search(r'(\d+)h', name)
    if hour_match:
        time_hours = int(hour_match.group(1))
    
    min_match = re.search(r'(\d+)min', name)
    if min_match:
        time_hours += int(min_match.group(1)) / 60
    
    result['time_hours'] = time_hours
    
    return result


# =============================================================================
# Peak Finding with Pseudo-Voigt Fitting
# =============================================================================

def _estimate_fwhm_simple(
    theta_range: np.ndarray,
    int_range: np.ndarray,
    idx_max: int,
    peak_int: float
) -> float:
    """Estimate FWHM using half-maximum method."""
    half_max = peak_int / 2
    
    left_idx = idx_max
    while left_idx > 0 and int_range[left_idx] > half_max:
        left_idx -= 1
    
    right_idx = idx_max
    while right_idx < len(int_range) - 1 and int_range[right_idx] > half_max:
        right_idx += 1
    
    return max(theta_range[right_idx] - theta_range[left_idx], 0.1)


def _fit_peak_pseudo_voigt(
    theta_range: np.ndarray,
    int_range: np.ndarray,
    peak_theta: float,
    peak_int: float,
    initial_fwhm: float
) -> Tuple[bool, float, float, float, float, float]:
    """Try Pseudo-Voigt fitting. Returns (success, theta, intensity, fwhm, eta, area)."""
    try:
        optimizer = LMOptimizer(max_iterations=500, tolerance=1e-6)
        
        initial_guess = PseudoVoigtParams(
            center=peak_theta,
            amplitude=peak_int,
            fwhm=initial_fwhm,
            eta=0.5
        )
        
        fit_result = optimizer.fit_single_peak(theta_range, int_range, initial_guess=initial_guess)
        
        if fit_result.success and fit_result.r_squared > 0.8:
            fitted_curve = PseudoVoigt.profile(
                theta_range, fit_result.params.center,
                fit_result.params.amplitude, fit_result.params.fwhm, fit_result.params.eta
            )
            try:
                area = np.trapezoid(fitted_curve, theta_range)
            except AttributeError:
                area = np.trapz(fitted_curve, theta_range)
            
            return (
                True,
                fit_result.params.center,
                fit_result.params.amplitude,
                fit_result.params.fwhm,
                fit_result.params.eta,
                area
            )
    except Exception:
        pass
    
    return False, peak_theta, peak_int, initial_fwhm, 0.5, 0.0


def _calculate_peak_area_simple(
    theta_range: np.ndarray,
    int_range: np.ndarray
) -> float:
    """Calculate peak area using trapezoidal integration."""
    try:
        return np.trapezoid(int_range, theta_range)
    except AttributeError:
        return np.trapz(int_range, theta_range)



def find_peak_in_range(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    center: float,
    window: float = 2.0,
    use_pv_fitting: bool = True
) -> Optional[PeakData]:
    """
    Find peak near expected position using Pseudo-Voigt fitting.
    
    Args:
        two_theta: 2θ array
        intensity: Intensity array
        center: Expected peak center position
        window: Search window (degrees)
        use_pv_fitting: If True, use Pseudo-Voigt fitting
        
    Returns:
        PeakData with fitted parameters, or None if no peak found
    """
    # Select range
    mask = (two_theta >= center - window) & (two_theta <= center + window)
    if not np.any(mask):
        return None
    
    theta_range = two_theta[mask]
    int_range = intensity[mask]
    
    # Find maximum
    idx_max = np.argmax(int_range)
    peak_theta = theta_range[idx_max]
    peak_int = int_range[idx_max]
    
    if peak_int < 50:
        return None
    
    # Estimate initial FWHM
    initial_fwhm = _estimate_fwhm_simple(theta_range, int_range, idx_max, peak_int)
    
    # Try Pseudo-Voigt fitting
    if use_pv_fitting:
        success, peak_theta, peak_int, fwhm, eta, area = _fit_peak_pseudo_voigt(
            theta_range, int_range, peak_theta, peak_int, initial_fwhm
        )
        if success:
            return PeakData(hkl=(0, 0, 0), two_theta=peak_theta, intensity=peak_int, fwhm=fwhm, area=area)
    
    # Fallback to simple method
    fwhm = max(initial_fwhm, 0.05)
    area = _calculate_peak_area_simple(theta_range, int_range)
    
    return PeakData(hkl=(0, 0, 0), two_theta=peak_theta, intensity=peak_int, fwhm=fwhm, area=area)


# =============================================================================
# Main Pipeline
# =============================================================================

class AXCSASPipeline:
    """
    Complete AXCSAS analysis pipeline.
    
    Integrates all Phase 04-07 analysis modules into a unified workflow.
    """
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        """Initialize pipeline with configuration."""
        self.config = config or AnalysisConfig()
        
        # Initialize analyzers
        self.scherrer = ScherrerCalculator()
        
        # Initialize Williamson-Hall analyzer
        # 初始化 Williamson-Hall 分析器
        self.wh = WilliamsonHallAnalyzer()
        
        # Initialize Texture analyzer
        # 初始化紋理分析器
        self.texture = TextureAnalyzer()
        self.sf_analyzer = StackingFaultAnalyzer()
        self.lattice = LatticeMonitor()
    
    def analyze(
        self,
        filepath: str,
        sample_age_hours: Optional[float] = None
    ) -> PipelineResult:
        """
        Run complete analysis on XRD data file.
        
        Args:
            filepath: Path to XRD data file
            sample_age_hours: Time since deposition (hours)
            
        Returns:
            PipelineResult with all analysis data
        """
        # Parse filename and initialize result
        file_info = parse_filename(filepath)
        effective_sample_age = sample_age_hours if sample_age_hours is not None else file_info['time_hours']
        
        result = PipelineResult(
            filepath=filepath,
            sample_name=file_info['name'],
            leveler_concentration=file_info['concentration_ml'],
            plating_time_hours=None,
            sample_age_hours=effective_sample_age,
        )
        
        # Step 1: Load data
        two_theta, intensity, error = self._load_and_validate_data(filepath)
        if error:
            result.report = error
            return result
        
        # Step 2: Find peaks
        result.peaks = self._find_peaks_from_data(two_theta, intensity)
        if len(result.peaks) < 2:
            result.report = f"Only {len(result.peaks)} peaks found"
            return result
        
        # Step 3: Scherrer analysis
        result.scherrer_results, result.average_size_nm = self._run_scherrer_analysis(result.peaks)
        
        # Step 4: W-H analysis (if enough peaks)
        if len(result.peaks) >= 3:
            two_theta_arr = np.array([p.two_theta for p in result.peaks])
            fwhm_arr = np.array([p.fwhm for p in result.peaks])
            hkl_list = [p.hkl for p in result.peaks]
            result.wh_result = self.wh.analyze(two_theta_arr, fwhm_arr, hkl_list)
        
        # Step 5: Texture analysis
        intensities = {p.hkl: p.intensity for p in result.peaks}
        result.texture_result = self.texture.analyze(intensities)
        
        # Step 6: Defect analysis
        result.stacking_fault, result.lattice_result = self._run_defect_analysis(result.peaks)
        
        # Step 7: Annealing state and comprehensive report
        result.annealing_state, _ = determine_annealing_state(sample_age_hours)
        result.comprehensive = self._build_comprehensive(result)
        result.report = generate_comprehensive_report(result.comprehensive)
        
        return result
    
    def _build_comprehensive(self, result: PipelineResult) -> ComprehensiveResult:
        """Build ComprehensiveResult from pipeline result."""
        comp = ComprehensiveResult(
            sample_name=result.sample_name,
            sample_age_hours=result.sample_age_hours,
        )
        
        # Scherrer
        if result.average_size_nm:
            comp.scherrer_size_nm = result.average_size_nm
            comp.scherrer_validity = "VALID"
        
        # W-H
        if result.wh_result:
            comp.wh_size_nm = result.wh_result.crystallite_size_nm
            comp.wh_strain = result.wh_result.microstrain
            comp.wh_r_squared = result.wh_result.r_squared
            comp.wh_quality = result.wh_result.quality_level.value
        
        # Texture (DATA ONLY, no diagnosis)
        if result.texture_result:
            comp.dominant_orientation = result.texture_result.dominant_hkl
            comp.dominant_tc = result.texture_result.dominant_tc
            comp.is_random_texture = result.texture_result.is_random
        
        # Defects
        if result.stacking_fault:
            comp.peak_separation_deg = result.stacking_fault.peak_separation_deg
            comp.stacking_fault_alpha = result.stacking_fault.alpha_percent
            comp.stacking_fault_severity = result.stacking_fault.severity.value
        
        if result.lattice_result:
            comp.lattice_constant = result.lattice_result.lattice_constant
            comp.lattice_status = result.lattice_result.status.value
        
        comp.annealing_state = result.annealing_state.value
        
        return comp


# =============================================================================
# Convenience Functions
# =============================================================================

def run_full_analysis(
    filepath: str,
    sample_age_hours: Optional[float] = None,
    config: Optional[AnalysisConfig] = None
) -> PipelineResult:
    """
    Run complete AXCSAS analysis on a single file.
    
    Example:
        >>> result = run_full_analysis("data/raw/sample.txt")
        >>> print(result.report)
    """
    pipeline = AXCSASPipeline(config)
    return pipeline.analyze(filepath, sample_age_hours)


def batch_analyze(
    filepaths: List[str],
    sample_age_hours: Optional[float] = None,
    config: Optional[AnalysisConfig] = None
) -> List[PipelineResult]:
    """
    Run analysis on multiple files.
    """
    pipeline = AXCSASPipeline(config)
    results = []
    
    for fp in filepaths:
        result = pipeline.analyze(fp, sample_age_hours)
        results.append(result)
    
    return results
