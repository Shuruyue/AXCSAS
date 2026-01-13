"""
AXCSAS Complete Analysis Pipeline
=================================

Unified pipeline integrating all analysis phases:
- Phase 02: Preprocessing
- Phase 03: Peak Fitting
- Phase 04: Scherrer Size Calculation
- Phase 05: Williamson-Hall Strain Analysis
- Phase 06: Texture Analysis
- Phase 07: Defect and Stress Diagnosis
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict, Any
from pathlib import Path
import re
import sys

# Add parent path
sys.path.insert(0, str(Path(__file__).parent.parent))

from axcsas.methods.scherrer import (
    ScherrerCalculatorEnhanced,
    ScherrerResultEnhanced,
    ValidityFlag,
)
from axcsas.methods.williamson_hall import (
    WilliamsonHallEnhanced,
    WHResultEnhanced,
    WHQualityLevel,
)
from axcsas.methods.texture import (
    TextureAnalyzerEnhanced,
    TextureResultEnhanced,
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
from axcsas.fitting.lm_optimizer import LMOptimizer, FitResult
from axcsas.fitting.pseudo_voigt import PseudoVoigt, PseudoVoigtParams


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class AnalysisConfig:
    """Configuration for AXCSAS analysis pipeline."""
    
    # X-ray parameters
    wavelength: float = 1.54056  # Cu Kα1
    
    # Scherrer parameters
    use_cubic_habit: bool = True
    
    # Peak detection
    peak_window: float = 2.0  # degrees around expected position
    min_intensity: float = 100  # minimum counts
    
    # Instrumental broadening (Caglioti U, V, W)
    caglioti_u: float = 0.0
    caglioti_v: float = 0.0
    caglioti_w: float = 0.01
    
    # JCPDS Cu peak positions
    peak_positions: Dict[Tuple[int, int, int], float] = field(default_factory=lambda: {
        (1, 1, 1): 43.3,
        (2, 0, 0): 50.4,
        (2, 2, 0): 74.1,
        (3, 1, 1): 89.9,
    })


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
    scherrer_results: List[ScherrerResultEnhanced] = field(default_factory=list)
    average_size_nm: Optional[float] = None
    
    # Phase 05: W-H
    wh_result: Optional[WHResultEnhanced] = None
    
    # Phase 06: Texture
    texture_result: Optional[TextureResultEnhanced] = None
    
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

def find_peak_in_range(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    center: float,
    window: float = 2.0,
    use_pv_fitting: bool = True
) -> Optional[PeakData]:
    """
    Find peak near expected position using Pseudo-Voigt fitting.
    
    使用 Pseudo-Voigt 擬合找到峰位附近的峰值。
    
    Args:
        two_theta: 2θ array
        intensity: Intensity array
        center: Expected peak center position
        window: Search window (degrees)
        use_pv_fitting: If True, use Pseudo-Voigt fitting; otherwise use simple method
        
    Returns:
        PeakData with fitted parameters, or None if no peak found
    """
    # Select range
    mask = (two_theta >= center - window) & (two_theta <= center + window)
    if not np.any(mask):
        return None
    
    theta_range = two_theta[mask]
    int_range = intensity[mask]
    
    # Find maximum for initial guess
    idx_max = np.argmax(int_range)
    peak_theta = theta_range[idx_max]
    peak_int = int_range[idx_max]
    
    # Check minimum intensity
    if peak_int < 50:
        return None
    
    # Default values (will be overwritten by fitting)
    fwhm = 0.2
    eta = 0.5
    fit_r_squared = 0.0
    area = 0.0
    
    if use_pv_fitting:
        # =====================================================================
        # Pseudo-Voigt Fitting (Phase 03)
        # =====================================================================
        try:
            optimizer = LMOptimizer(max_iterations=500, tolerance=1e-6)
            
            # Create initial guess based on simple FWHM estimate
            half_max = peak_int / 2
            left_idx = idx_max
            while left_idx > 0 and int_range[left_idx] > half_max:
                left_idx -= 1
            right_idx = idx_max
            while right_idx < len(int_range) - 1 and int_range[right_idx] > half_max:
                right_idx += 1
            initial_fwhm = max(theta_range[right_idx] - theta_range[left_idx], 0.1)
            
            initial_guess = PseudoVoigtParams(
                center=peak_theta,
                amplitude=peak_int,
                fwhm=initial_fwhm,
                eta=0.5  # Start with mixed Gaussian-Lorentzian
            )
            
            fit_result = optimizer.fit_single_peak(
                theta_range, int_range, initial_guess=initial_guess
            )
            
            if fit_result.success and fit_result.r_squared > 0.8:
                # Use fitted parameters
                peak_theta = fit_result.params.center
                peak_int = fit_result.params.amplitude
                fwhm = fit_result.params.fwhm
                eta = fit_result.params.eta
                fit_r_squared = fit_result.r_squared
                
                # Calculate integrated area from the fit
                fitted_curve = PseudoVoigt.profile(
                    theta_range, peak_theta, peak_int, fwhm, eta
                )
                try:
                    area = np.trapezoid(fitted_curve, theta_range)
                except AttributeError:
                    area = np.trapz(fitted_curve, theta_range)
            else:
                # Fallback to simple method if fitting fails
                use_pv_fitting = False
                
        except Exception:
            # Fallback to simple method
            use_pv_fitting = False
    
    if not use_pv_fitting:
        # =====================================================================
        # Simple Half-Maximum Method (Fallback)
        # =====================================================================
        half_max = peak_int / 2
        
        left_idx = idx_max
        while left_idx > 0 and int_range[left_idx] > half_max:
            left_idx -= 1
        right_idx = idx_max
        while right_idx < len(int_range) - 1 and int_range[right_idx] > half_max:
            right_idx += 1
        
        fwhm = theta_range[right_idx] - theta_range[left_idx]
        fwhm = max(fwhm, 0.05)  # Minimum FWHM (instrument limit)
        
        try:
            area = np.trapezoid(int_range, theta_range)
        except AttributeError:
            area = np.trapz(int_range, theta_range)
    
    return PeakData(
        hkl=(0, 0, 0),  # Will be assigned later
        two_theta=peak_theta,
        intensity=peak_int,
        fwhm=fwhm,
        area=area
    )


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
        self.scherrer = ScherrerCalculatorEnhanced(
            wavelength=self.config.wavelength,
            use_cubic_habit=self.config.use_cubic_habit
        )
        self.wh = WilliamsonHallEnhanced(
            wavelength=self.config.wavelength
        )
        self.texture = TextureAnalyzerEnhanced()
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
        # Parse filename
        file_info = parse_filename(filepath)
        
        result = PipelineResult(
            filepath=filepath,
            sample_name=file_info['name'],
            leveler_concentration=file_info['concentration_ml'],
            plating_time_hours=file_info['time_hours'],
            sample_age_hours=sample_age_hours,
        )
        
        # 1. Load data
        try:
            two_theta, intensity = load_bruker_txt(filepath)
        except Exception as e:
            result.report = f"Error loading file: {e}"
            return result
        
        if len(two_theta) == 0:
            result.report = "No data found in file"
            return result
        
        # 2. Find peaks
        for hkl, expected_pos in self.config.peak_positions.items():
            peak = find_peak_in_range(
                two_theta, intensity, expected_pos, self.config.peak_window
            )
            if peak:
                peak.hkl = hkl
                result.peaks.append(peak)
        
        if len(result.peaks) < 2:
            result.report = f"Only {len(result.peaks)} peaks found"
            return result
        
        # 3. Phase 04: Scherrer analysis
        for peak in result.peaks:
            scherrer_result = self.scherrer.calculate(
                two_theta=peak.two_theta,
                fwhm_observed=peak.fwhm,
                fwhm_instrumental=np.sqrt(self.config.caglioti_w),
                hkl=peak.hkl
            )
            result.scherrer_results.append(scherrer_result)
        
        # Calculate average size
        valid_sizes = [
            r.size_nm for r in result.scherrer_results
            if r.validity_flag != ValidityFlag.UNRELIABLE
        ]
        if valid_sizes:
            result.average_size_nm = np.mean(valid_sizes)
        
        # 4. Phase 05: W-H analysis (if enough peaks)
        if len(result.peaks) >= 3:
            two_theta_arr = np.array([p.two_theta for p in result.peaks])
            fwhm_arr = np.array([p.fwhm for p in result.peaks])
            hkl_list = [p.hkl for p in result.peaks]
            
            result.wh_result = self.wh.analyze(
                two_theta_arr, fwhm_arr, hkl_list
            )
        
        # 5. Phase 06: Texture analysis
        intensities = {p.hkl: p.intensity for p in result.peaks}
        result.texture_result = self.texture.analyze(intensities)
        
        # 6. Phase 07: Defect analysis
        peak_111 = next((p for p in result.peaks if p.hkl == (1, 1, 1)), None)
        peak_200 = next((p for p in result.peaks if p.hkl == (2, 0, 0)), None)
        
        if peak_111 and peak_200:
            result.stacking_fault = self.sf_analyzer.analyze(
                peak_111.two_theta, peak_200.two_theta
            )
        
        # Lattice from (311) or (220) preferably
        high_angle_peak = next(
            (p for p in result.peaks if p.hkl in [(3, 1, 1), (2, 2, 0)]),
            result.peaks[-1] if result.peaks else None
        )
        if high_angle_peak:
            result.lattice_result = self.lattice.analyze_lattice(
                high_angle_peak.two_theta, high_angle_peak.hkl
            )
        
        # Self-annealing state
        result.annealing_state, _ = determine_annealing_state(sample_age_hours)
        
        # 7. Generate comprehensive result
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
