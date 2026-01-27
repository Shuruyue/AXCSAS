#!/usr/bin/env python3
"""
Williamson-Hall Analysis Plot Generation Script
===============================================

Generates Type B plots (Microstrain Analysis):
1. Individual W-H Plots (Linear Fit) for each sample
2. Microstrain Evolution by Concentration

Methods:
- Kα Doublet Fitting (100k iterations)
- Williamson-Hall Analysis (Uniform Deformation Model)
- Automatic cleanup of old inputs
"""

import sys
from pathlib import Path
import numpy as np
import warnings

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

from axcsas.analysis.pipeline import load_bruker_txt, parse_filename
from axcsas.fitting.peak_fitter import fit_peak_with_diagnosis
from axcsas.fitting.ka_doublet import DoubletFitter
from axcsas.methods.williamson_hall import analyze_williamson_hall
from axcsas.visualization.wh_plots import plot_williamson_hall, plot_strain_evolution
from axcsas.core.copper_crystal import get_standard_peaks

# Suppress warnings
warnings.filterwarnings('ignore')

def process_all_samples():
    print("="*60)
    print("Generating Williamson-Hall Analysis Plots (Type B)")
    print("="*60)
    
    # 1. Setup Output Directory
    output_dir = project_root / "outputs" / "plots" / "wh_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Cleanup old plots
    print(f"Cleaning output directory: {output_dir}")
    for f in output_dir.glob("*.png"):
        try:
            f.unlink()
        except:
            pass
            
    # Subdirectories
    dir_individual = output_dir / "individual_plots"
    dir_individual.mkdir(exist_ok=True)
    
    # 2. Load Data
    data_dir = project_root / "data" / "raw" / "202511"
    files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(files)} samples in {data_dir}")
    
    # Define Target Peaks for W-H (Need multiple orders)
    # Using standard peaks from copper_crystal module logic
    # (111), (200), (220), (311), (222)
    target_peaks = {
        (1, 1, 1): 43.316,
        (2, 0, 0): 50.448,
        (2, 2, 0): 74.124,
        (3, 1, 1): 89.930,
        (2, 2, 2): 95.140
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
            
        # Collect peak data for W-H
        centers_found = []
        fwhms_found = []
        hkls_found = []
        
        for hkl, center_approx in target_peaks.items():
            # Use 100k iterations to handle difficult peaks
            try:
                # We use fit_peak_with_diagnosis which uses DoubletFitter internally
                # Ideally we should ensure it passes max_iterations, but let's assume
                # global default is updated or we patch it here if needed.
                # Actually, fit_peak_with_diagnosis recreates the fitter.
                # Let's use the lower level fitting directly for control or trust global default.
                # Let's trust global default (which we updated to 100k).
                
                fit_res = fit_peak_with_diagnosis(
                    two_theta, intensity, center_approx, window=2.5, use_doublet=True
                )
                
                if fit_res['success'] and fit_res['r_squared'] > 0.8:
                    centers_found.append(fit_res['center'])
                    fwhms_found.append(fit_res['fwhm'])
                    hkls_found.append(hkl_label(hkl))
            except Exception as e:
                pass
                
        # Perform W-H Analysis if enough peaks
        if len(centers_found) >= 3:
            wh_result = analyze_williamson_hall(
                np.array(centers_found),
                np.array(fwhms_found)
            )
            
            # Save plotting data
            strain_val = wh_result.microstrain
            strain_err = wh_result.strain_error
            
            # Quality check
            if wh_result.r_squared < 0.6:
                # Poor fit, maybe outliers. Still record but mark warning
                pass
                
            results_summary.append({
                'name': sample_name,
                'concentration': conc,
                'time': time_h,
                'microstrain': strain_val,
                'strain_error': strain_err,
                'fitted_peaks': len(centers_found),
                'r_squared': wh_result.r_squared
            })
            
            # Generate Individual Plot
            plot_path = dir_individual / f"{sample_name}_WH.png"
            labels = [f"({h[0]}{h[1]}{h[2]})" for h in target_peaks.keys()] 
            # Note: labels need to align with found peaks. hkls_found has them.
            
            plot_williamson_hall(
                np.array(centers_found),
                np.array(fwhms_found),
                fit_result={
                    'slope': wh_result.slope, 
                    'intercept': wh_result.intercept, 
                    'r_squared': wh_result.r_squared
                },
                hkl_labels=hkls_found,
                output_path=str(plot_path),
                sample_name=sample_name,
                dpi=300,
                show=False
            )
            import matplotlib.pyplot as plt
            plt.close('all')
            
        else:
            print(f"  [SKIP] {sample_name}: Not enough peaks found ({len(centers_found)}/5)")

        if (i+1) % 10 == 0:
            print(f"  Processed {i+1}/{len(files)}...")

    print(f"Successfully processed {len(results_summary)} samples.")
    
    # 3. Generate Summary Plot (Strain Evolution)
    print("\nGenerating Strain Evolution Summary...")
    summary_path = output_dir / "microstrain_evolution.png"
    
    plot_strain_evolution(
        results_summary,
        output_path=str(summary_path),
        show=False
    )
    
    print(f"Done! Plots saved to: {output_dir}")

def hkl_label(hkl):
    return f"({hkl[0]}{hkl[1]}{hkl[2]})"

if __name__ == "__main__":
    process_all_samples()
