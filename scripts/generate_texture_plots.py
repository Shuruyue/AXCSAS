#!/usr/bin/env python3
"""
Texture Analysis Plot Generation Script (Type C)
================================================

Generates:
1. Polar Plots (Radar Charts) for Texture Coefficients
2. Texture Evolution Plots (TC vs Time/Concentration)
3. Texture Fraction Evolution (f = TC/3)

Methods:
- Harris Texture Coefficient (TC)
- Integrated Area Calculation (from Pseudo-Voigt params)
- Automatic cleanup of old inputs
"""

import sys
from pathlib import Path
import numpy as np
import warnings
import matplotlib.pyplot as plt

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

from axcsas.analysis.pipeline import load_bruker_txt, parse_filename
from axcsas.fitting.peak_fitter import fit_peak_with_diagnosis
from axcsas.methods.texture import analyze_texture
from axcsas.visualization.texture_plots import plot_texture_polar, plot_tc_evolution, plot_texture_fraction_single

# Suppress warnings
warnings.filterwarnings('ignore')

def calculate_area(amplitude, fwhm, eta):
    """
    Calculate integrated area of a Pseudo-Voigt peak.
    Area = Amp * FWHM * [eta * (pi/2) + (1-eta) * sqrt(pi/(4*ln2))]
    """
    # Constants
    CONST_LORENTZ = 1.570796327  # pi/2
    CONST_GAUSS = 1.064467       # sqrt(pi/(4*ln2))
    
    shape_factor = eta * CONST_LORENTZ + (1 - eta) * CONST_GAUSS
    return amplitude * fwhm * shape_factor

def process_all_samples():
    print("="*60)
    print("Generating Texture Analysis Plots (Type C)")
    print("="*60)
    
    # 1. Setup Output Directory
    output_dir = project_root / "outputs" / "plots" / "texture_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Cleanup old plots
    print(f"Cleaning output directory: {output_dir}")
    for f in output_dir.glob("*.png"):
        try:
            f.unlink()
        except:
            pass
            
    # Subdirectories
    dir_polar = output_dir / "polar_plots"
    dir_polar.mkdir(exist_ok=True)
    
    # 2. Load Data
    data_dir = project_root / "data" / "raw" / "202511"
    files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(files)} samples in {data_dir}")
    
    # Define Target Peaks for Texture Analysis
    # (111), (200), (220) are the primary face-centered cubic peaks
    target_peaks = {
        (1, 1, 1): 43.316,
        (2, 0, 0): 50.448,
        (2, 2, 0): 74.124
    }
    
    results_summary = []
    
    print("\nProcessing samples...")
    for i, filepath in enumerate(files):
        file_info = parse_filename(str(filepath))
        conc = file_info.get('concentration_ml', 0)
        time_h = file_info.get('time_hours', 0)
        sample_name = filepath.stem
        
        try:
            two_theta, intensity = load_bruker_txt(str(filepath))
        except:
            continue
            
        # Collect integrated intensities
        intensities_area = {} # Area-based (Correct for texture)
        
        for hkl, center_approx in target_peaks.items():
            try:
                # Use 100k iterations
                fit_res = fit_peak_with_diagnosis(
                    two_theta, intensity, center_approx, window=2.5, use_doublet=True
                )
                
                if fit_res['success'] and fit_res['r_squared'] > 0.8:
                    amp = fit_res['amplitude']
                    fwhm = fit_res['fwhm']
                    eta = fit_res['eta']
                    
                    # Calculate Area of Kα1
                    area_ka1 = calculate_area(amp, fwhm, eta)
                    
                    # Total Area = Area(Kα1) + Area(Kα2)
                    # Area(Kα2) scales with Amp(Kα2) = 0.5 * Amp(Kα1)
                    # So Total Area ≈ 1.5 * Area(Kα1)
                    total_area = area_ka1 * 1.5
                    
                    intensities_area[hkl] = total_area
            except Exception as e:
                pass
                
        # Perform Texture Analysis
        if len(intensities_area) >= 2:
            # We need at least 2 peaks to be meaningful, standard is 3
            
            # Analyze using Area
            tc_result = analyze_texture(intensities_area, use_area=True)
            
            # Store results
            result_entry = {
                'name': sample_name,
                'concentration': conc,
                'time': time_h,
                'tc_values': {hkl_label(k): v for k, v in tc_result.tc_values.items()},
                'dominant': hkl_label(tc_result.dominant_hkl) if tc_result.dominant_hkl else "None",
                'is_random': tc_result.is_random,
                'sigma': tc_result.degree_of_texture,
                'high_quality': True # Can add filter based on R² later
            }
            results_summary.append(result_entry)
            
            # Generate Individual Polar Plot
            plot_path = dir_polar / f"{sample_name}_TC.png"
            plot_texture_polar(
                result_entry['tc_values'],
                output_path=str(plot_path),
                sample_name=sample_name,
                dpi=300,
                show=False
            )
            plt.close('all')
            
        else:
            print(f"  [SKIP] {sample_name}: Not enough peaks found ({len(intensities_area)}/3)")

        if (i+1) % 10 == 0:
            print(f"  Processed {i+1}/{len(files)}...")

    print(f"Successfully processed {len(results_summary)} samples.")
    
    # 3. Generate Summary Plots
    print("\nGenerating Evolution Plots...")
    
    # (A) Evolution of TC (Standard)
    # (A) Evolution of TC (Standard) - Panel: Direction
    print("  - Generating TC Evolution (by Direction)...")
    plot_tc_evolution(
        results_summary,
        x_param='time',
        output_path=str(output_dir / "texture_evolution_TC_by_direction.png"),
        normalize=False,
        dpi=1200,
        show=False
    )
    
    # (B) Evolution of Fraction - Panel: Direction
    print("  - Generating Fraction Evolution (by Direction)...")
    plot_tc_evolution(
        results_summary,
        x_param='time',
        output_path=str(output_dir / "texture_evolution_Fraction_by_direction.png"),
        normalize=True,
        dpi=1200,
        show=False
    )
    
    # (C) Fraction Comparison - Panel: Concentration
    print("  - Generating Fraction Comparison (by Concentration)...")
    plot_texture_fraction_single(
        results_summary,
        x_param='time',
        output_path=str(output_dir / "texture_fraction_by_concentration.png"),
        metric="fraction",
        dpi=1200,
        show=False
    )
    
    # (D) TC Comparison - Panel: Concentration
    print("  - Generating TC Comparison (by Concentration)...")
    plot_texture_fraction_single(
        results_summary,
        x_param='time',
        output_path=str(output_dir / "texture_TC_by_concentration.png"),
        metric="tc",
        dpi=1200,
        show=False
    )
    
    print(f"Done! Plots saved to: {output_dir}")

def hkl_label(hkl):
    return f"({hkl[0]}{hkl[1]}{hkl[2]})"

if __name__ == "__main__":
    process_all_samples()
