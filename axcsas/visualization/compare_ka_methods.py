#!/usr/bin/env python3
"""
Kα Doublet Fitting Comparison Script
=====================================

Compare three fitting methods:
1. Standard Pseudo-Voigt (ignores Kα₂)
2. Kα₂ stripping + Pseudo-Voigt  
3. Kα₁/Kα₂ doublet fitting

Generates comparison plots for each sample.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from axcsas.analysis.pipeline import load_bruker_txt, parse_filename
from axcsas.fitting.lm_optimizer import LMOptimizer
from axcsas.fitting.pseudo_voigt import PseudoVoigt
from axcsas.fitting.ka_doublet import Ka2Stripper, DoubletFitter, calculate_ka2_position

# Cu peak positions
PEAKS = {
    (1, 1, 1): 43.3,
    (2, 0, 0): 50.4,
    (2, 2, 0): 74.1,
}
PEAK_LABELS = ['(111)', '(200)', '(220)']


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
    fig.suptitle(f'Kα Doublet Fitting Comparison: {sample_name}', fontsize=14, fontweight='bold')
    
    results = []
    
    for row, ((hkl, expected_pos), label) in enumerate(zip(PEAKS.items(), PEAK_LABELS)):
        window = 2.5
        mask = (two_theta >= expected_pos - window) & (two_theta <= expected_pos + window)
        theta_region = two_theta[mask]
        int_region = intensity[mask]
        
        if len(theta_region) < 10:
            continue
        
        # Method 1: Standard Pseudo-Voigt
        ax1 = axes[row, 0]
        optimizer = LMOptimizer()
        result1 = optimizer.fit_single_peak(theta_region, int_region)
        
        ax1.scatter(theta_region, int_region, s=8, alpha=0.5, c='gray', label='Data')
        if result1.success:
            fitted1 = PseudoVoigt.profile(theta_region, result1.params.center, 
                                          result1.params.amplitude, result1.params.fwhm, 
                                          result1.params.eta)
            # Add estimated background
            bg = np.min(int_region)
            ax1.plot(theta_region, fitted1 + bg, 'b-', lw=2, label='Pseudo-Voigt')
            ax1.set_title(f'{label} - Standard PV\nFWHM={result1.params.fwhm:.4f}°, R²={result1.r_squared:.4f}')
        ax1.set_xlabel('2θ (°)')
        ax1.set_ylabel('Intensity')
        ax1.legend(fontsize=7)
        ax1.grid(True, alpha=0.3)
        
        # Method 2: Kα₂ Stripping
        ax2 = axes[row, 1]
        stripper = Ka2Stripper()
        _, stripped = stripper.strip(theta_region, int_region)
        result2 = optimizer.fit_single_peak(theta_region, stripped)
        
        ax2.scatter(theta_region, int_region, s=8, alpha=0.3, c='gray', label='Original')
        ax2.scatter(theta_region, stripped, s=8, alpha=0.6, c='blue', label='Kα₂ stripped')
        if result2.success:
            fitted2 = PseudoVoigt.profile(theta_region, result2.params.center,
                                          result2.params.amplitude, result2.params.fwhm,
                                          result2.params.eta)
            bg = np.min(stripped)
            ax2.plot(theta_region, fitted2 + bg, 'r-', lw=2, label='Fit (stripped)')
            ax2.set_title(f'{label} - Kα₂ Stripped\nFWHM={result2.params.fwhm:.4f}°, R²={result2.r_squared:.4f}')
        ax2.set_xlabel('2θ (°)')
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3)
        
        # Method 3: Doublet Fitting
        ax3 = axes[row, 2]
        doublet_fitter = DoubletFitter()
        result3 = doublet_fitter.fit(theta_region, int_region, expected_pos)
        
        ax3.scatter(theta_region, int_region, s=8, alpha=0.5, c='gray', label='Data')
        if result3.success and result3.fitted_curve is not None:
            ax3.plot(theta_region, result3.fitted_curve, 'g-', lw=2, label='Doublet fit')
            # Mark Kα₁ and Kα₂ positions
            ax3.axvline(result3.center_ka1, color='blue', ls='--', alpha=0.5, label=f'Kα₁={result3.center_ka1:.3f}°')
            ax3.axvline(result3.center_ka2, color='orange', ls='--', alpha=0.5, label=f'Kα₂={result3.center_ka2:.3f}°')
            ax3.set_title(f'{label} - Doublet Fit\nFWHM={result3.fwhm:.4f}°, R²={result3.r_squared:.4f}')
        ax3.set_xlabel('2θ (°)')
        ax3.legend(fontsize=7)
        ax3.grid(True, alpha=0.3)
        
        results.append({
            'peak': label,
            'pv_fwhm': result1.params.fwhm if result1.success else None,
            'pv_r2': result1.r_squared if result1.success else None,
            'stripped_fwhm': result2.params.fwhm if result2.success else None,
            'stripped_r2': result2.r_squared if result2.success else None,
            'doublet_fwhm': result3.fwhm if result3.success else None,
            'doublet_r2': result3.r_squared if result3.success else None,
        })
    
    plt.tight_layout()
    output_path = output_dir / f'{sample_name}_ka_comparison.png'
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close()
    
    return results


def main():
    print("=" * 60)
    print("Kα Doublet Fitting Comparison")
    print("=" * 60)
    
    project_root = Path(__file__).parent
    data_dir = project_root / "data" / "raw" / "202511"
    output_dir = project_root / "outputs" / "plots" / "ka_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Only process a few samples for comparison
    txt_files = sorted(data_dir.glob("*.txt"))[:10]  # First 10 samples
    
    print(f"\nProcessing {len(txt_files)} samples...")
    print(f"Output: {output_dir}\n")
    
    all_results = []
    for filepath in txt_files:
        print(f"Processing: {filepath.name}")
        result = compare_methods_for_sample(filepath, output_dir)
        if result:
            all_results.append({'file': filepath.name, 'results': result})
            print(f"  ✓ Saved")
    
    # Print summary table
    print("\n" + "=" * 80)
    print("Summary: R² Comparison")
    print("=" * 80)
    print(f"{'File':<30} {'Peak':<8} {'PV R²':<10} {'Strip R²':<10} {'Doublet R²':<10}")
    print("-" * 80)
    for item in all_results:
        for r in item['results']:
            pv_r2 = f"{r['pv_r2']:.4f}" if r['pv_r2'] is not None else "N/A"
            strip_r2 = f"{r['stripped_r2']:.4f}" if r['stripped_r2'] is not None else "N/A"
            doublet_r2 = f"{r['doublet_r2']:.4f}" if r['doublet_r2'] is not None else "N/A"
            print(f"{item['file'][:28]:<30} {r['peak']:<8} {pv_r2:<10} {strip_r2:<10} {doublet_r2:<10}")
    
    print("\n" + "=" * 60)
    print(f"Complete! Generated {len(all_results)} comparison plots")
    print("=" * 60)


if __name__ == "__main__":
    main()
