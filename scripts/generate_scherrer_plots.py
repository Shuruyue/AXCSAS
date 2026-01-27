#!/usr/bin/env python3
"""
Scherrer Analysis Plot Generation Script
========================================

Generates Crystallite Size (Type A) plots:
1. Size Evolution by Peak (Size vs Time)
2. Size Evolution by Concentration (Size vs Time, grouped by Conc)

Methods:
- Peak Fitting: Kα Doublet Fitting (High Accuracy)
- Size Calculation: Scherrer Equation (FWHM-based)
- K-Factor: Langford & Wilson 1978 (Cubic Habit Directional)
  (111): 0.855, (200): 0.886, (220): 0.834
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
from axcsas.methods.scherrer import ScherrerCalculator
from axcsas.core.copper_crystal import get_standard_peaks
from axcsas.visualization.scherrer_plots import (
    plot_scherrer_evolution_by_peak,
    plot_scherrer_by_concentration
)
from axcsas.core.config import ParameterConfig

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

def process_all_samples():
    print("="*60)
    print("Generating Scherrer Crystallite Size Plots")
    print("="*60)
    
    # 1. Setup Configuration
    # Use default config (Instrument Broadening W=0.003 -> FWHM=0.055)
    config = ParameterConfig()
    
    # Scherrer Calculator with Cubic Habit (L&W 1978 constants)
    # This automatically selects 0.855/0.886/0.834 based on hkl
    calculator = ScherrerCalculator(
        use_cubic_habit=True,
        caglioti_params=(
            config.instrument.caglioti_u,
            config.instrument.caglioti_v,
            config.instrument.caglioti_w
        )
    )
    
    # 2. Load Data
    data_dir = project_root / "data" / "raw" / "202511"
    files = sorted(data_dir.glob("*.txt"))
    
    print(f"Found {len(files)} samples in {data_dir}")
    
    # Peaks to analyze (Low angle peaks are most reliable for size)
    # (111), (200), (220) are standard for Cu size analysis
    target_peaks = {
        (1, 1, 1): 43.316,
        (2, 0, 0): 50.448,
        (2, 2, 0): 74.124
    }
    
    results = []
    
    print("Processing samples...")
    for i, filepath in enumerate(files):
        # Parse basic info
        file_info = parse_filename(str(filepath))
        conc = file_info.get('concentration_ml', 0)
        time_h = file_info.get('time_hours', 0)
        
        # Load data
        try:
            two_theta, intensity = load_bruker_txt(str(filepath))
        except Exception as e:
            print(f"Skipping {filepath.name}: {e}")
            continue
            
        sample_data = {
            'filename': filepath.name,
            'concentration': conc,
            'time': time_h,
            'peaks': []
        }
        
        # Analyze each peak
        for hkl, center in target_peaks.items():
            # 1. Fit Peak (Method 3: Doublet Fitting for maximum accuracy)
            fit_res = fit_peak_with_diagnosis(
                two_theta, intensity, center, window=2.0, use_doublet=True
            )
            
            if fit_res['success'] and fit_res['r_squared'] > 0.8:
                # 2. Calculate Size
                # fit_res['fwhm'] is the FWHM of the Kα1 component (Observed)
                # ScherrerCalculator will handle instrumental subtraction
                
                scherrer_res = calculator.calculate(
                    two_theta=fit_res['center'],
                    fwhm_observed=fit_res['fwhm'],
                    hkl=hkl,
                    # We pass None for fwhm_instrumental to let calculator use Caglioti params
                    # derived from config above
                )
                
                if scherrer_res.size_nm > 0: # REMOVED UPPER LIMIT < 300 per user request
                    sample_data['peaks'].append({
                        'hkl': hkl,
                        'two_theta': fit_res['center'],
                        'size_nm': scherrer_res.size_nm,
                        'size_err': 0.0,
                        'k_factor': scherrer_res.k_factor
                    })
                else:
                    print(f"  [SKIP] {filepath.name} {hkl}: Size {scherrer_res.size_nm:.1f} nm <= 0")
            else:
                reason = "Fit Failed" if not fit_res['success'] else f"Low R2 ({fit_res.get('r_squared', 0):.4f})"
                if fit_res.get('amplitude', 0) < 50: reason = "Low Intensity"
                print(f"  [SKIP] {filepath.name} {hkl}: {reason}")
        
        results.append(sample_data)
        if (i+1) % 10 == 0:
            print(f"  Processed {i+1}/{len(files)}...")

    print(f"Successfully processed {len(results)} samples.")
    
    # 3. Generate Plots
    output_dir = project_root / "outputs" / "plots" / "fwhm_analysis" # Or create new specific folder?
    # Let's create a dedicated folder for Size Analysis
    size_output_dir = project_root / "outputs" / "plots" / "scherrer_size"
    size_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Cleanup old plots
    for f in size_output_dir.glob("*.png"):
        try:
            f.unlink()
        except Exception:
            pass
    print(f"Cleaned up old files in {size_output_dir}")
    
    print(f"\nGenerating plots in: {size_output_dir}")
    
    # Plot 1: Evolution by Peak (Separate plots for 111, 200, 220)
    print("  - Generating Size Evolution by Peak...")
    plot_scherrer_evolution_by_peak(
        results, 
        output_path=str(size_output_dir / "size_evolution_by_peak.png"),
        show=False
    )
    
    # Plot 2: Evolution by Concentration (Subplots for each concentration)
    print("  - Generating Size Evolution by Concentration...")
    plot_scherrer_by_concentration(
        results,
        output_path=str(size_output_dir / "size_evolution_by_concentration.png"),
        show=False
    )
    
    print("Done!")

if __name__ == "__main__":
    process_all_samples()
