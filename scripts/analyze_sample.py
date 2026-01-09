#!/usr/bin/env python3
"""
Single Sample Analysis Script
Performs complete XRD analysis on a single sample.

Usage:
    python analyze_sample.py --input data/raw/sample.xy --output outputs/results/
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.preprocessing import (
    load_xrd_data, 
    apply_smoothing, 
    subtract_background,
    strip_kalpha2
)
from src.fitting import detect_peaks, fit_peaks_lm, PseudoVoigt
from src.physics import (
    ScherrerCalculator,
    WilliamsonHallAnalyzer,
    TextureAnalyzer,
    CagliotiCorrection
)
from src.validation import assess_fit_quality, ErrorAnalyzer
from src.utils import CU_KA1, SCHERRER_K, CU_JCPDS


def load_config(config_path: str = "config.yaml") -> dict:
    """Load configuration file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Warning: Config file not found, using defaults")
        return {}


def analyze_sample(
    input_file: str,
    output_dir: str,
    config: dict = None,
    generate_plots: bool = True
):
    """
    Perform complete XRD analysis on a single sample.
    
    Args:
        input_file: Path to XRD data file
        output_dir: Directory for output files
        config: Configuration dictionary
        generate_plots: Whether to generate plots
    """
    config = config or {}
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    sample_name = Path(input_file).stem
    
    print("=" * 60)
    print(f"AXCSAS Sample Analysis: {sample_name}")
    print("=" * 60)
    
    # =========================================================================
    # Step 1: Load Data
    # =========================================================================
    print("\n[1/6] Loading XRD data...")
    two_theta, intensity = load_xrd_data(input_file)
    print(f"      Data range: {two_theta[0]:.2f}° - {two_theta[-1]:.2f}°")
    print(f"      Data points: {len(two_theta)}")
    
    # =========================================================================
    # Step 2: Preprocessing
    # =========================================================================
    print("\n[2/6] Preprocessing...")
    
    # Smoothing
    preproc_config = config.get('preprocessing', {})
    smooth_config = preproc_config.get('smoothing', {})
    window = smooth_config.get('window_size', 11)
    poly = smooth_config.get('poly_order', 3)
    
    intensity_smooth = apply_smoothing(intensity, window_size=window, poly_order=poly)
    
    # Background subtraction
    intensity_bg = subtract_background(two_theta, intensity_smooth)
    
    # Kα2 stripping (optional)
    if preproc_config.get('kalpha_strip', {}).get('enable', True):
        intensity_clean = strip_kalpha2(two_theta, intensity_bg)
    else:
        intensity_clean = intensity_bg
    
    print("      ✓ Smoothing applied")
    print("      ✓ Background subtracted")
    
    # =========================================================================
    # Step 3: Peak Detection and Fitting
    # =========================================================================
    print("\n[3/6] Peak detection and fitting...")
    
    fit_config = config.get('fitting', {})
    peak_config = fit_config.get('peak_detection', {})
    min_height = peak_config.get('min_height', 100)
    min_distance = peak_config.get('min_distance', 0.5)
    
    peaks = detect_peaks(
        two_theta, intensity_clean,
        min_height=min_height,
        min_distance_deg=min_distance
    )
    
    fit_results = fit_peaks_lm(two_theta, intensity_clean, peaks)
    
    print(f"      Found {len(fit_results)} peaks")
    
    # =========================================================================
    # Step 4: Crystallite Size Calculation
    # =========================================================================
    print("\n[4/6] Calculating crystallite size...")
    
    phys_config = config.get('physical_constants', {})
    wavelength = phys_config.get('wavelength', CU_KA1)
    
    scherrer_config = phys_config.get('scherrer_k', {})
    k_factor = scherrer_config.get('default', SCHERRER_K.default)
    
    scherrer = ScherrerCalculator(wavelength=wavelength, k_factor=k_factor)
    
    # Check for Caglioti calibration
    inst_config = config.get('instrument', {})
    caglioti_config = inst_config.get('caglioti', {})
    
    has_calibration = all(
        caglioti_config.get(p) is not None 
        for p in ['U', 'V', 'W']
    )
    
    size_results = []
    for peak in fit_results:
        if has_calibration:
            caglioti = CagliotiCorrection(
                U=caglioti_config['U'],
                V=caglioti_config['V'],
                W=caglioti_config['W']
            )
            fwhm_inst = caglioti.calculate_fwhm_inst(peak.center)
            result = scherrer.calculate_with_correction(
                peak.center, peak.fwhm, fwhm_inst
            )
        else:
            result = scherrer.calculate(peak.center, peak.fwhm)
        
        size_results.append({
            'two_theta': peak.center,
            'fwhm': peak.fwhm,
            'eta': peak.eta,
            'size_nm': result.size_nm,
            'reliable': result.is_reliable
        })
    
    avg_size = np.mean([r['size_nm'] for r in size_results if r['reliable']])
    
    print(f"      Average crystallite size: {avg_size:.1f} nm")
    
    # =========================================================================
    # Step 5: Williamson-Hall Analysis
    # =========================================================================
    print("\n[5/6] Williamson-Hall analysis...")
    
    wh_analyzer = WilliamsonHallAnalyzer(wavelength=wavelength, k_factor=k_factor)
    
    peak_positions = np.array([r['two_theta'] for r in size_results])
    peak_fwhms = np.array([r['fwhm'] for r in size_results])
    
    wh_result = wh_analyzer.analyze(peak_positions, peak_fwhms)
    
    print(f"      W-H crystallite size: {wh_result.crystallite_size_nm:.1f} nm")
    print(f"      Microstrain: {wh_result.microstrain:.2e}")
    print(f"      R²: {wh_result.r_squared:.4f}")
    
    # =========================================================================
    # Step 6: Texture Analysis
    # =========================================================================
    print("\n[6/6] Texture analysis...")
    
    texture_analyzer = TextureAnalyzer()
    
    # Match peaks to hkl
    intensities = {}
    for peak in fit_results:
        hkl = texture_analyzer.get_hkl_for_angle(peak.center, tolerance=1.0)
        if hkl:
            intensities[hkl] = peak.amplitude
    
    if len(intensities) >= 2:
        texture_result = texture_analyzer.analyze(intensities)
        print(f"      Degree of texture: {texture_result.degree_of_texture:.3f}")
        if texture_result.is_random:
            print("      Texture: Random (powder average)")
        else:
            pref = texture_result.preferred_orientation
            print(f"      Preferred orientation: ({pref[0]}{pref[1]}{pref[2]})")
    else:
        texture_result = None
        print("      Insufficient peaks for texture analysis")
    
    # =========================================================================
    # Save Results
    # =========================================================================
    print("\n" + "-" * 60)
    print("Saving results...")
    
    # Peak data CSV
    peak_df = pd.DataFrame(size_results)
    peak_csv = output_path / f"{sample_name}_peaks.csv"
    peak_df.to_csv(peak_csv, index=False)
    print(f"  ✓ Peak data: {peak_csv}")
    
    # Summary results
    summary = {
        'sample': sample_name,
        'timestamp': datetime.now().isoformat(),
        'peaks_detected': len(fit_results),
        'scherrer_size_nm': {
            'average': avg_size,
            'values': [r['size_nm'] for r in size_results]
        },
        'williamson_hall': {
            'size_nm': wh_result.crystallite_size_nm,
            'microstrain': wh_result.microstrain,
            'r_squared': wh_result.r_squared
        }
    }
    
    if texture_result:
        summary['texture'] = {
            'degree': texture_result.degree_of_texture,
            'is_random': texture_result.is_random,
            'tc_values': {
                f"({h}{k}{l})": v 
                for (h, k, l), v in texture_result.tc_values.items()
            }
        }
    
    summary_yaml = output_path / f"{sample_name}_summary.yaml"
    with open(summary_yaml, 'w') as f:
        yaml.dump(summary, f, default_flow_style=False)
    print(f"  ✓ Summary: {summary_yaml}")
    
    # =========================================================================
    # Generate Plots (if requested)
    # =========================================================================
    if generate_plots:
        print("\nGenerating plots...")
        
        # Fitted spectrum plot
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(two_theta, intensity_clean, 'b-', label='Data', alpha=0.7)
        
        # Plot fitted peaks
        fitted = np.zeros_like(two_theta)
        for peak in fit_results:
            peak_curve = PseudoVoigt.profile(
                two_theta, peak.center, peak.amplitude, peak.fwhm, peak.eta
            )
            fitted += peak_curve
        
        ax.plot(two_theta, fitted, 'r-', label='Fitted', linewidth=2)
        ax.set_xlabel('2θ (degrees)')
        ax.set_ylabel('Intensity (a.u.)')
        ax.set_title(f'XRD Pattern: {sample_name}')
        ax.legend()
        
        plot_path = output_path / f"{sample_name}_fit.png"
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  ✓ Fit plot: {plot_path}")
    
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Analyze a single XRD sample"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to XRD data file"
    )
    parser.add_argument(
        "--output", "-o",
        default="outputs/results",
        help="Output directory (default: outputs/results)"
    )
    parser.add_argument(
        "--config", "-c",
        default="config.yaml",
        help="Configuration file (default: config.yaml)"
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation"
    )
    
    args = parser.parse_args()
    
    config = load_config(args.config)
    analyze_sample(
        args.input,
        args.output,
        config,
        generate_plots=not args.no_plots
    )


if __name__ == "__main__":
    main()
