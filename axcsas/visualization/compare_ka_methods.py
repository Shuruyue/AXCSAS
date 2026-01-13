#!/usr/bin/env python3
"""
Kα Fitting Method Comparison Script
====================================

Compare three fitting methods:
1. Single Pseudo-Voigt (ignores Kα₂ splitting)
2. Kα₁/Kα₂ Doublet fitting (shows both peaks)
3. Enhanced Pseudo-Voigt with polynomial background (current fitting method)

Generates comparison plots for ALL samples.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import least_squares

from axcsas.analysis.pipeline import load_bruker_txt, parse_filename
from axcsas.fitting.lm_optimizer import LMOptimizer
from axcsas.fitting.pseudo_voigt import PseudoVoigt
from axcsas.fitting.ka_doublet import DoubletFitter

# Cu peak positions
PEAKS = {
    (1, 1, 1): 43.3,
    (2, 0, 0): 50.4,
    (2, 2, 0): 74.1,
}
PEAK_LABELS = ['(111)', '(200)', '(220)']


def enhanced_pv_fit(theta_range, int_range, expected_center):
    """
    Enhanced Pseudo-Voigt fitting with polynomial background.
    Same method as fitting_diagnosis.py.
    """
    # Find peak position
    peak_idx = np.argmax(int_range)
    peak_theta = theta_range[peak_idx]
    peak_int = int_range[peak_idx]
    
    # Estimate FWHM from half-max
    half_max = peak_int / 2
    above_half = int_range > half_max
    if np.sum(above_half) > 2:
        indices = np.where(above_half)[0]
        initial_fwhm = theta_range[indices[-1]] - theta_range[indices[0]]
        initial_fwhm = max(0.1, min(initial_fwhm, 1.0))
    else:
        initial_fwhm = 0.3
    
    # Estimate polynomial background
    n_edge = max(5, len(theta_range) // 8)
    bg_x = np.concatenate([theta_range[:n_edge], theta_range[-n_edge:]])
    bg_y = np.concatenate([int_range[:n_edge], int_range[-n_edge:]])
    bg_coeffs = np.polyfit(bg_x, bg_y, 2)
    
    # Model: Pseudo-Voigt + quadratic background
    def enhanced_pv_model(x, center, amplitude, fwhm, eta, a2, a1, a0):
        pv = PseudoVoigt.profile(x, center, amplitude, fwhm, np.clip(eta, 0, 1))
        background = a2 * x**2 + a1 * x + a0
        return pv + background
    
    def residual(params):
        return int_range - enhanced_pv_model(theta_range, *params)
    
    # Multi-start optimization
    best_r2 = -1
    best_params = None
    best_fitted = None
    
    for eta_init in [0.3, 0.5, 0.7, 0.9]:
        x0 = np.array([peak_theta, peak_int * 0.9, initial_fwhm, eta_init, 
                      bg_coeffs[0], bg_coeffs[1], bg_coeffs[2]])
        
        lb = np.array([theta_range.min(), 10, 0.03, 0.0, -np.inf, -np.inf, -np.inf])
        ub = np.array([theta_range.max(), peak_int * 1.5, 3.0, 1.0, np.inf, np.inf, np.inf])
        
        try:
            opt_result = least_squares(residual, x0, bounds=(lb, ub), 
                                       max_nfev=5000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
            
            if opt_result.success or opt_result.status >= 1:
                fitted = enhanced_pv_model(theta_range, *opt_result.x)
                ss_res = np.sum((int_range - fitted)**2)
                ss_tot = np.sum((int_range - np.mean(int_range))**2)
                r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                
                if r2 > best_r2:
                    best_r2 = r2
                    best_params = opt_result.x
                    best_fitted = fitted
        except Exception:
            continue
    
    if best_params is not None:
        return {
            'success': True,
            'center': best_params[0],
            'fwhm': best_params[2],
            'eta': best_params[3],
            'r_squared': best_r2,
            'fitted_curve': best_fitted,
        }
    return {'success': False}


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
        window = 2.5
        mask = (two_theta >= expected_pos - window) & (two_theta <= expected_pos + window)
        theta_region = two_theta[mask]
        int_region = intensity[mask]
        
        if len(theta_region) < 10:
            continue
        
        # ==== Method 1: Simple Pseudo-Voigt (no Kα splitting) ====
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
        ax1.set_xlabel('2θ (°)')
        ax1.set_ylabel('Intensity')
        ax1.legend(fontsize=7)
        ax1.grid(True, alpha=0.3)
        
        # ==== Method 2: Kα Doublet (shows both peaks) ====
        ax2 = axes[row, 1]
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
        ax2.set_xlabel('2θ (°)')
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3)
        
        # ==== Method 3: Enhanced Pseudo-Voigt (polynomial bg) ====
        ax3 = axes[row, 2]
        result3 = enhanced_pv_fit(theta_region, int_region, expected_pos)
        
        ax3.scatter(theta_region, int_region, s=8, alpha=0.5, c='gray', label='Data')
        if result3['success']:
            ax3.plot(theta_region, result3['fitted_curve'], 'r-', lw=2, label='Enhanced PV')
            ax3.set_title(f'{label} - Enhanced PV\nFWHM={result3["fwhm"]:.4f}°, R²={result3["r_squared"]:.4f}')
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
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close()
    
    return results


def main():
    print("=" * 60)
    print("Fitting Method Comparison")
    print("1. Single Pseudo-Voigt (ignores Kα₂)")
    print("2. Kα Doublet (Kα₁ + Kα₂ separate peaks)")
    print("3. Enhanced PV (polynomial background)")
    print("=" * 60)
    
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data" / "raw" / "202511"
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
    
    # Print summary table
    print("\n" + "=" * 90)
    print("Summary: R² Comparison")
    print("=" * 90)
    print(f"{'File':<30} {'Peak':<8} {'Single R²':<12} {'Doublet R²':<12} {'Enhanced R²':<12}")
    print("-" * 90)
    for item in all_results:
        for r in item['results']:
            single_r2 = f"{r['single_r2']:.4f}" if r['single_r2'] is not None else "N/A"
            doublet_r2 = f"{r['doublet_r2']:.4f}" if r['doublet_r2'] is not None else "N/A"
            enhanced_r2 = f"{r['enhanced_r2']:.4f}" if r['enhanced_r2'] is not None else "N/A"
            print(f"{item['file'][:28]:<30} {r['peak']:<8} {single_r2:<12} {doublet_r2:<12} {enhanced_r2:<12}")
    
    print("\n" + "=" * 60)
    print(f"Complete! Generated {len(all_results)} comparison plots")
    print("=" * 60)


if __name__ == "__main__":
    main()
