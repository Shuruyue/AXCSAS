"""
Generate Defect Analysis Plots (Type D)
=======================================

Generates high-precision plots for:
1. Residual Stress Evolution (Anisotropic Model)
   - Uses direction-dependent Young's Modulus and Poisson Ratio.
2. Stacking Fault Probability Evolution (Warren Method)
   - Uses peak shift correlations.

This script implements the complete pipeline: Reading Raw Data -> Fitting -> Analysis -> Plotting.
"""

import sys
import yaml
import json
import logging
import numpy as np
from pathlib import Path
from typing import List, Dict, Any

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from axcsas.core.copper_crystal import (
    get_youngs_modulus, 
    get_poisson_ratio, 
    get_standard_peaks
)
from axcsas.methods.defect_analysis import (
    LatticeMonitor, 
    analyze_stacking_faults,
    StressType
)
from axcsas.visualization.defect_plots import (
    plot_stress_evolution,
    plot_stacking_fault_evolution
)
from axcsas.analysis.pipeline import load_bruker_txt, parse_filename
from axcsas.fitting.peak_fitter import fit_peak_with_diagnosis

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def load_config(config_path: Path):
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def process_raw_files(root_dir: Path) -> List[Dict]:
    """
    Process raw text files to calculate Stress and Stacking Faults.
    """
    monitor = LatticeMonitor()
    processed_count = 0
    samples = []
    
    # Locate data directory
    # Try both absolute and relative paths
    data_dir = root_dir / "data" / "raw" / "202511"
    if not data_dir.exists():
        # Try finding it elsewhere
        data_dir = Path("d:/Shuru/Git Project/AXCSAS/data/raw/202511")
        
    if not data_dir.exists():
        logger.error(f"Data directory not found: {data_dir}")
        return []
        
    files = sorted(list(data_dir.glob("*.txt")))
    logger.info(f"Found {len(files)} raw files in {data_dir}")
    
    # Targets for fitting
    # (111) approx 43.3, (200) approx 50.4
    target_peaks = {
        (1, 1, 1): 43.316,
        (2, 0, 0): 50.448,
        (2, 2, 0): 74.124
    }
    
    print("Processing samples (Fitting & Defect Analysis)...")
    
    for i, filepath in enumerate(files):
        try:
            # Parse metadata
            file_info = parse_filename(str(filepath))
            conc = file_info.get('concentration_ml', 0)
            time_h = file_info.get('time_hours', 0)
            
            # Load Data
            two_theta, intensity = load_bruker_txt(str(filepath))
            
            # Fit Peaks
            fitted_peaks = {}
            for hkl, center_approx in target_peaks.items():
                res = fit_peak_with_diagnosis(
                    two_theta, intensity, center_approx, window=2.5, use_doublet=True
                )
                if res['success']:
                    fitted_peaks[hkl] = res
            
            # Use results only if we have at least one peak, preferably both
            if not fitted_peaks:
                continue
                
            sample_entry = {
                'concentration': conc,
                'time': time_h,
                'sample_name': filepath.stem,
                'high_quality': True, 
                'stress_results': {},
                'stacking_fault': None
            }
            
            # 1. Calculate Residual Stress (Anisotropic)
            for hkl, peak_res in fitted_peaks.items():
                pos = peak_res['center']
                h, k, l = hkl
                hkl_str = f"({h}{k}{l})"
                
                # Anisotropic constants
                E = get_youngs_modulus(h, k, l)
                # Note: defect_analysis uses get_poisson_ratio internally with use_directional=True
                # We just pass E.
                
                stress_res = monitor.estimate_residual_stress(pos, hkl, E)
                
                sample_entry['stress_results'][hkl_str] = {
                    'stress_mpa': stress_res.stress_mpa,
                    'stress_type': stress_res.stress_type.value,
                    'd_measured': stress_res.d_measured,
                    'd_standard': stress_res.d_standard
                }
                
            # 2. Calculate Stacking Faults (Needs both peaks)
            if (1,1,1) in fitted_peaks and (2,0,0) in fitted_peaks:
                pos_111 = fitted_peaks[(1,1,1)]['center']
                pos_200 = fitted_peaks[(2,0,0)]['center']
                
                sf_res = analyze_stacking_faults(pos_111, pos_200)
                
                sample_entry['stacking_fault'] = {
                    'alpha_probability': sf_res.alpha_probability,
                    'alpha_percent': sf_res.alpha_percent,
                    'severity': sf_res.severity.value,
                    'peak_separation': sf_res.peak_separation_deg
                }
            
            samples.append(sample_entry)
            processed_count += 1
            
            if (i+1) % 10 == 0:
                print(f"  Processed {i+1}/{len(files)}...")

        except Exception as e:
            # logger.warning(f"Failed to process {filepath.name}: {e}")
            pass
            
    logger.info(f"Successfully processed {processed_count} samples.")
    return samples

def main():
    print("="*60)
    print("Generating Defect Analysis Plots (Type D)")
    print("="*60)
    
    root_dir = Path(__file__).parent.parent
    config = load_config(root_dir / "config.yaml")
    
    # Paths
    base_output_dir = Path(config['output']['paths']['root'])
    plots_dir = base_output_dir / "plots" / "defect_analysis"
    
    # Create directory if not exists
    if plots_dir.exists():
        print(f"Cleaning output directory: {plots_dir}")
        for f in plots_dir.glob("*.png"):
            f.unlink()
    else:
        plots_dir.mkdir(parents=True, exist_ok=True)
        
    # Process Data
    print("\nStarting Analysis Pipeline...")
    print("  - Using Anisotropic Young's Modulus (E_111=191, E_200=67, E_220=130 GPa)")
    print("  - Using Anisotropic Poisson's Ratio (ν_111=0.27, ν_200=0.42, ν_220=0.34)")
    print("  - Calculating Warren Stacking Fault Probability")
    
    samples = process_raw_files(root_dir)
    
    if not samples:
        logger.error("No samples processed. Exiting.")
        return
    
    # Generate Plots
    print("\nGenerating Plots...")
    
    # 1. Stress Evolution
    print("  - Generating Residual Stress Evolution...")
    plot_stress_evolution(
        samples,
        output_path=str(plots_dir / "stress_evolution.png"),
        dpi=1200,
        show=False
    )
    
    # 2. Stacking Fault Evolution
    print("  - Generating Stacking Fault Probability...")
    plot_stacking_fault_evolution(
        samples,
        output_path=str(plots_dir / "stacking_fault_evolution.png"),
        dpi=1200,
        show=False
    )
    
    print(f"\nDone! Plots saved to: {plots_dir}")

if __name__ == "__main__":
    main()
