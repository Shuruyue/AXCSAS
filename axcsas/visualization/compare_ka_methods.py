#!/usr/bin/env python3
"""
Fitting Method Comparison Script
=================================

Compare three fitting methods:
1. Simple Pseudo-Voigt (ignores Kα₂ splitting)
2. Kα Doublet fitting (shows both peaks)
3. Rigorous Pseudo-Voigt (primary method, same as fitting_diagnosis)

Ensures consistency with the main analysis pipeline.
確保與主分析管道的一致性。
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from axcsas.analysis.pipeline import load_bruker_txt, parse_filename
from axcsas.fitting.pseudo_voigt import PseudoVoigt
from axcsas.fitting.lm_optimizer import LMOptimizer
from axcsas.fitting.ka_doublet import DoubletFitter

# Import fitting function directly from diagnosis script for consistency
# (fitting_api.py was removed; it caused R² regression)
from axcsas.visualization.generate_fitting_diagnosis import fit_peak_with_diagnosis

# Cu peak positions - use unified function
from axcsas.core.copper_crystal import get_standard_peaks

# Peak configuration
PEAKS = get_standard_peaks()  # Returns (111), (200), (220)
PEAK_LABELS = [
    f"({h[0]}{h[1]}{h[2]})" for h in PEAKS.keys()
]


def compare_methods_for_sample(filepath: Path, output_dir: Path):
    """Generate comparison plot for one sample."""
    
    # Load data
    try:
        two_theta, intensity = load_bruker_txt(str(filepath))
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None
    
    file_info = parse_filename(str(filepath))
    sample_name = filepath.stem
    
    # Create figure
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    conc = file_info.get('concentration_ml', 0)
    time_h = file_info.get('time_hours', 0)
    fig.suptitle(
        f'Fitting Method Comparison: {sample_name}\n'
        f'(Leveler: {conc:.1f} mL/1.5L, Annealing Time: {time_h:.0f}h)', 
        fontsize=14, fontweight='bold'
    )
    
    results = []
    
    for row, ((hkl, expected_pos), label) in enumerate(zip(PEAKS.items(), PEAK_LABELS)):
        window = 2.5  # Consistent with pipeline default
        mask = (two_theta >= expected_pos - window) & (two_theta <= expected_pos + window)
        theta_region = two_theta[mask]
        int_region = intensity[mask]
        
        if len(theta_region) < 10:
            continue
        
        # ==== Method 1: Simple Pseudo-Voigt (Single Peak) ====
        ax1 = axes[row, 0]
        optimizer = LMOptimizer()
        result1 = optimizer.fit_single_peak(theta_region, int_region)
        
        ax1.scatter(theta_region, int_region, s=8, alpha=0.5, c='gray', label='Data')
        if result1.success:
            bg = np.min(int_region)
            fitted1 = PseudoVoigt.profile(theta_region, result1.params.center, 
                                          result1.params.amplitude, result1.params.fwhm, 
                                          result1.params.eta) + bg
            ax1.plot(theta_region, fitted1, 'b-', lw=2, label='Fit')
            ax1.set_title(f'{label} - Single Peak\nFWHM={result1.params.fwhm:.4f}°, R²={result1.r_squared:.4f}')
        else:
            ax1.set_title(f'{label} - Single Peak\n(Fit failed)')
        ax1.set_xlabel('2θ (°)')
        ax1.set_ylabel('Intensity')
        ax1.legend(fontsize=7)
        ax1.grid(True, alpha=0.3)
        
        # ==== Method 2: Kα Doublet (shows both peaks) ====
        ax2 = axes[row, 1]
        
        # Use simple DoubletFitter here just for visualization comparison
        # (This is distinct from the "Enhanced" method used in Method 3)
        doublet_fitter = DoubletFitter()
        result2 = doublet_fitter.fit(theta_region, int_region, expected_pos)
        
        ax2.scatter(theta_region, int_region, s=8, alpha=0.5, c='gray', label='Data')
        if result2.success and result2.fitted_curve is not None:
            ax2.plot(theta_region, result2.fitted_curve, 'g-', lw=2, label='Doublet fit')
            # Draw individual Kα₁ and Kα₂ peaks
            bg2 = np.min(result2.fitted_curve)
            ka1_peak = PseudoVoigt.profile(theta_region, result2.center_ka1, 
                                           result2.amplitude_ka1, result2.fwhm, result2.eta) + bg2
            ka2_peak = PseudoVoigt.profile(theta_region, result2.center_ka2,
                                           result2.amplitude_ka1 * 0.5, result2.fwhm, result2.eta) + bg2
            ax2.plot(theta_region, ka1_peak, 'b--', lw=1.5, alpha=0.7, label=f'Kα₁={result2.center_ka1:.3f}°')
            ax2.plot(theta_region, ka2_peak, 'r--', lw=1.5, alpha=0.7, label=f'Kα₂={result2.center_ka2:.3f}°')
            ax2.set_title(f'{label} - Doublet (Kα₁+Kα₂)\nFWHM={result2.fwhm:.4f}°, R²={result2.r_squared:.4f}')
        else:
            ax2.set_title(f'{label} - Doublet\n(Fit failed)')
        ax2.set_xlabel('2θ (°)')
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3)
        
        # ==== Method 3: Rigorous Pseudo-Voigt (SAME AS FITTING DIAGNOSIS) ====
        ax3 = axes[row, 2]
        
        # CALL THE EXACT FUNCTION FROM generate_fitting_diagnosis.py
        result3 = fit_peak_with_diagnosis(
            two_theta, intensity, expected_pos, use_doublet=True  # Same as fitting_diagnosis
        )
        
        ax3.scatter(theta_region, int_region, s=8, alpha=0.5, c='gray', label='Data')
        if result3['success'] and result3['fitted_curve'] is not None:
            curve3 = result3['fitted_curve']
            # fit_peak_with_diagnosis returns curve for its internal theta_range
            # we need to make sure we plot against the correct x-axis
            x3 = result3['theta_range']
            
            ax3.plot(x3, curve3, 'r-', lw=2, label='Enhanced PV')
            ax3.set_title(f'{label} - Enhanced PV\nFWHM={result3["fwhm"]:.4f}°, R²={result3["r_squared"]:.4f}')
        else:
             ax3.set_title(f'{label} - Enhanced PV\n(Fit Failed)')

        ax3.set_xlabel('2θ (°)')
        ax3.legend(fontsize=7)
        ax3.grid(True, alpha=0.3)
        
        results.append({
            'peak': label,
            'single_fwhm': result1.params.fwhm if result1.success else None,
            'single_r2': result1.r_squared if result1.success else None,
            'doublet_fwhm': result2.fwhm if result2.success else None,
            'doublet_r2': result2.r_squared if result2.success else None,
            'enhanced_fwhm': result3['fwhm'] if result3['success'] else None,
            'enhanced_r2': result3['r_squared'] if result3['success'] else None,
        })
    
    plt.tight_layout()
    output_path = output_dir / f'{sample_name}_method_comparison.png'
    plt.savefig(output_path, bbox_inches='tight', dpi=2400)
    plt.close()
    
    return results


def main():
    print("=" * 60)
    print("Fitting Method Comparison (using shared fitting_api)")
    print("1. Simple Pseudo-Voigt (ignores Kα₂)")
    print("2. Kα Doublet (Kα₁ + Kα₂ separate peaks)")
    print("3. Rigorous PV (SAME as fitting_diagnosis)")
    print("=" * 60)
    
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data" / "raw" / "202511" # Update this path as needed
    if not data_dir.exists():
        data_dir = project_root / "data" / "raw"
    
    output_dir = project_root / "outputs" / "plots" / "ka_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process ALL samples
    txt_files = sorted(data_dir.glob("*.txt"))
    
    print(f"\nProcessing {len(txt_files)} samples...")
    print(f"Output: {output_dir}\n")
    
    all_results = []
    for filepath in txt_files:
        print(f"Processing: {filepath.name}")
        result = compare_methods_for_sample(filepath, output_dir)
        if result:
            all_results.append({'file': filepath.name, 'results': result})
            print(f"  ✓ Saved")
    
    # Print summary
    print("\n" + "=" * 60)
    print(f"Complete! Generated {len(all_results)} comparison plots")
    print("All methods use shared fitting_api for data consistency.")
    print("=" * 60)


if __name__ == "__main__":
    main()
