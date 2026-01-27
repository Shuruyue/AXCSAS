#!/usr/bin/env python3
"""
AXCSAS Plot Generator
=====================

Generate FWHM evolution and Scherrer Size plots from XRD data.
使用新的 visualization 模組生成 FWHM 演化圖與 Scherrer 尺寸圖。

Refactored to use axcsas.visualization module.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import re

from axcsas.methods.scherrer import ValidityFlag, calculate_scherrer
from axcsas.analysis.pipeline import AXCSASPipeline, AnalysisConfig, PipelineResult, parse_filename, load_bruker_txt

# Import new visualization module
from axcsas.visualization.style import (
    apply_axcsas_style,
    COLORBLIND_SAFE,
    save_figure,
)
from axcsas.visualization.fwhm_plots import (
    plot_fwhm_evolution,
    plot_fwhm_by_peak,
    plot_fwhm_by_concentration,
)
from axcsas.visualization.scherrer_plots import (
    plot_scherrer_evolution_by_peak,
    plot_scherrer_by_concentration,
)

# Apply unified style
apply_axcsas_style()


@dataclass
class SampleData:
    """Parsed sample data."""
    filepath: str
    filename: str
    concentration_ml: float
    time_hours: float
    result: Optional[PipelineResult] = None


# parse_sample_info removed - use parse_filename from pipeline


def analyze_all_samples(data_dir: Path) -> List[SampleData]:
    """Analyze all XRD samples in directory."""
    samples = []
    pipeline = AXCSASPipeline()
    
    txt_files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(txt_files)} XRD data files")
    
    for filepath in txt_files:
        try:
            file_info = parse_filename(str(filepath))
            concentration = file_info.get('concentration_ml', 0.0)
            time_hours = file_info.get('time_hours', 0.0)
            
            # Run AXCSAS analysis
            result = pipeline.analyze(str(filepath))
            
            sample = SampleData(
                filepath=str(filepath),
                filename=filepath.name,
                concentration_ml=concentration,
                time_hours=time_hours,
                result=result
            )
            samples.append(sample)
            print(f"  ✓ {filepath.name}")
            
        except Exception as e:
            print(f"  ✗ {filepath.name}: {e}")
    
    return samples

# Import fitting function from diagnosis for consistency
from axcsas.core.copper_crystal import get_standard_peaks
from axcsas.visualization.generate_fitting_diagnosis import fit_peak_with_diagnosis

# Get standard peaks
PEAK_POSITIONS = get_standard_peaks()

def convert_samples_to_plot_data(samples: List[SampleData]) -> List[Dict]:
    """
    Convert SampleData list to format expected by visualization module.
    轉換 SampleData 列表為視覺化模組所需格式。
    
    Refitted using fit_peak_with_diagnosis to ensure consistency with diagnosis plots.
    """
    plot_data = []
    
    print("  Refitting peaks for consistency with diagnosis plots...", flush=True)
    
    for sample in samples:
        try:
            # Reload data to perform consistent fitting
            two_theta, intensity = load_bruker_txt(sample.filepath)
        except Exception as e:
            print(f"    Warning: Could not reload {sample.filename}: {e}")
            continue
            
        peaks_data = []
        
        # Use valid peaks from the pipeline to know which HKLs are present, 
        # then fit them using the rigorous method from generate_fitting_diagnosis
        # Actually, let's stick to the standard set of peaks defined in PEAK_POSITIONS
        # to match the diagnosis plots exactly.
        
        for hkl_tuple, expected_pos in PEAK_POSITIONS.items():
             res = fit_peak_with_diagnosis(
                 two_theta, intensity, expected_pos, 
                 window=2.5,
                 use_doublet=True  # Same as fitting_diagnosis for consistency
             )
             
             hkl_str = f"({hkl_tuple[0]}{hkl_tuple[1]}{hkl_tuple[2]})"
             
             if res['success']:
                 # Calculate Scherrer size for consistency
                 sr = calculate_scherrer(
                     two_theta=res['center'], 
                     fwhm_observed=res['fwhm'], 
                     fwhm_instrumental=0.05
                 )
                 
                 # Calculate size error based on FWHM error
                 # D = Kλ / (β cosθ)  => |dD/dβ| = D / β
                 # Error(D) = D * Error(β) / β
                 # Use sample broadening (sr.fwhm_sample_rad) for β, but fwhm_error is in degrees
                 # Convert error to radians: err_rad = np.radians(fwhm_error)
                 # Error(D) = D * err_rad / beta_rad = D * fwhm_error_deg / fwhm_sample_deg
                 
                 fwhm_err = res.get('fwhm_error', 0.0)
                 size_err = 0.0
                 if sr.size_nm and sr.fwhm_sample > 0 and fwhm_err > 0:
                     size_err = sr.size_nm * (fwhm_err / sr.fwhm_sample)
                 
                 peaks_data.append({
                     'hkl': hkl_str,
                     'fwhm': res['fwhm'],
                     'fwhm_error': fwhm_err,
                     'intensity': res['amplitude'],
                     'size_nm': sr.size_nm,
                     'size_err': size_err,
                     'validity': sr.validity_flag.value,
                     'fit_quality': 'high' if not res.get('low_quality', False) else 'low',
                 })
             else:
                 # Fallback: use pipeline result if enhanced fitting fails
                 # This ensures stable-state samples with very narrow (220) peaks are included
                 if sample.result and sample.result.peaks:
                     for pipeline_peak in sample.result.peaks:
                         if pipeline_peak.hkl == hkl_tuple:
                             # Calculate fallback size
                             sr_fallback = calculate_scherrer(
                                 two_theta=pipeline_peak.two_theta,
                                 fwhm_observed=pipeline_peak.fwhm,
                                 fwhm_instrumental=0.05
                             )
                             peaks_data.append({
                                 'hkl': hkl_str,
                                 'fwhm': pipeline_peak.fwhm,
                                 'intensity': pipeline_peak.intensity,
                                 'size_nm': sr_fallback.size_nm,
                                 'validity': sr_fallback.validity_flag.value,
                                 'fit_quality': 'fallback',  # Mark as fallback
                             })
                             break

        # Consistently calculate Scherrer size using the DIAGNOSIS FWHM
        # This ensures data consistency with the FWHM plots
        for p_data in peaks_data:
            # We assume use_cubic_habit=True (default) and instrument broadening=0.05
            # We must pass the EXPECTED position for accurate K-factor calculation
            # But calculate_scherrer calculates K based on observed 2theta if HKL is auto-assigned?
            # Better to let it assign HKL based on 2theta.
            # However, we already KNOW the HKL. 
            # Note: calculate_scherrer takes (two_theta, fwhm, fwhm_instrumental)
            # We should probably pass the observed FWHM.
            
            # Find the peak position (approximate from expected)
            # Actually, `fit_peak_with_diagnosis` returns 'center' in res['center'] if we captured it?
            # It seems fit_peak_with_diagnosis only returns:
            # {'success': True, 'amplitude': ..., 'center': ..., 'sigma': ..., 'fwhm': ..., 'r2': ...}
            # Let's verify what fit_peak_with_diagnosis returns. 
            # It returns a dict.
            
            # Wait, `res` is available inside the loop. But here we are iterating over `peaks_data` which is a processed list.
            # We should calculate sizing inside the loop where `res` is available.
            pass

        # Let's rewrite the loop above to include size calculation immediately

        
        plot_data.append({
            'name': sample.filename,
            'concentration': sample.concentration_ml,
            'time': sample.time_hours,
            'peaks': peaks_data,
        })
    
    return plot_data


def generate_fwhm_plots(samples: List[SampleData], output_dir: Path) -> int:
    """
    Generate FWHM evolution plots using new visualization module.
    使用新視覺化模組生成 FWHM 演化圖。
    """
    plot_data = convert_samples_to_plot_data(samples)
    
    if not plot_data:
        print("No valid data for FWHM plots")
        return 0
    
    count = 0
    
    # Plot 1: FWHM Evolution by Time (Annealing Time, grouped by concentration)
    try:
        fig = plot_fwhm_evolution(
            plot_data,
            x_param='time',
            output_path=str(output_dir / 'fwhm_evolution_by_time.png'),
            show=False,
            dpi=1200,
            instrument_limit=0.05  # degrees, user-defined or from Caglioti fit
                                    # 使用者自定義或從 Caglioti 擬合計算
        )
        count += 1
        print(f"  ✓ fwhm_evolution_by_time.png")
    except Exception as e:
        print(f"  ✗ fwhm_evolution_by_time.png: {e}")
    
    # Plot 2: FWHM by Concentration (4 subplots, 3 peak lines each)
    try:
        fig = plot_fwhm_by_concentration(
            plot_data,
            output_path=str(output_dir / 'fwhm_by_concentration.png'),
            show=False,
            dpi=1200,
            instrument_limit=0.05  # degrees, user-defined or from Caglioti fit
        )
        count += 1
        print(f"  ✓ fwhm_by_concentration.png")
    except Exception as e:
        print(f"  ✗ fwhm_by_concentration.png: {e}")
    
    return count


def generate_scherrer_plots(samples: List[SampleData], output_dir: Path) -> int:
    """
    Generate Scherrer crystallite size plots using new visualization module.
    使用新視覺化模組生成 Scherrer 晶粒尺寸圖。
    """
    plot_data = convert_samples_to_plot_data(samples)
    
    if not plot_data:
        print("No valid data for Scherrer plots")
        return 0
    
    count = 0
    
    # Plot 1: Scherrer evolution by peak (direction)
    try:
        fig = plot_scherrer_evolution_by_peak(
            plot_data,
            output_path=str(output_dir / 'scherrer_size_by_direction.png'),
            show=False,
            dpi=1200,
        )
        count += 1
        print(f"  ✓ scherrer_size_by_direction.png")
    except Exception as e:
        print(f"  ✗ scherrer_size_by_direction.png: {e}")
    
    # Plot 2: Scherrer by concentration (4 subplots, square)
    try:
        fig = plot_scherrer_by_concentration(
            plot_data,
            output_path=str(output_dir / 'scherrer_size_by_concentration.png'),
            show=False,
            dpi=1200,
        )
        count += 1
        print(f"  ✓ scherrer_size_by_concentration.png")
    except Exception as e:
        print(f"  ✗ scherrer_size_by_concentration.png: {e}")
    
    return count


def main():
    """Main entry point."""
    print("=" * 60)
    print("AXCSAS Plot Generator (Refactored)")
    print("Using axcsas.visualization module")
    print("=" * 60)
    
    # Setup paths - go up 2 levels from visualization/ to project root
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data" / "raw" / "202511"
    output_dir = project_root / "outputs" / "plots"
    
    # Create output directories
    fwhm_dir = output_dir / "fwhm"
    scherrer_dir = output_dir / "scherrer"
    fwhm_dir.mkdir(parents=True, exist_ok=True)
    scherrer_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze all samples
    print(f"\nAnalyzing samples from: {data_dir}")
    samples = analyze_all_samples(data_dir)
    
    if not samples:
        print("No samples analyzed!")
        return 1
    
    print(f"\nSuccessfully analyzed {len(samples)} samples")
    
    # Generate FWHM plots
    print("\n" + "-" * 40)
    print("Generating FWHM Evolution Plots...")
    print("-" * 40)
    fwhm_count = generate_fwhm_plots(samples, fwhm_dir)
    
    # Generate Scherrer plots
    print("\n" + "-" * 40)
    print("Generating Scherrer Size Plots...")
    print("-" * 40)
    scherrer_count = generate_scherrer_plots(samples, scherrer_dir)
    
    print("\n" + "=" * 60)
    print(f"Complete! Generated {fwhm_count + scherrer_count} plots")
    print(f"Output directory: {output_dir}")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    exit(main())
