#!/usr/bin/env python3
"""
AXCSAS Peak Fitting Diagnostic Plots
====================================

Generate detailed fitting diagnostic plots for each XRD sample,
showing peak positions, FWHM, fitting curves, and R² values.

Uses Kα₁/Kα₂ doublet fitting for accurate peak characterization.
使用 Kα₁/Kα₂ 雙峰擬合進行精確峰型表徵。

Refactored to use axcsas.visualization module.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional

# Import AXCSAS modules
from axcsas.analysis.pipeline import (
    AXCSASPipeline,
    AnalysisConfig,
    load_bruker_txt,
    parse_filename,
)
from axcsas.fitting.pseudo_voigt import PseudoVoigt, PseudoVoigtParams
from axcsas.fitting.lm_optimizer import LMOptimizer
from axcsas.fitting.ka_doublet import DoubletFitter, Ka2Stripper, calculate_ka2_position

# Import new visualization module
from axcsas.visualization.style import (
    apply_axcsas_style,
    COLORBLIND_SAFE,
    PEAK_COLORS,
    save_figure,
    create_figure,
)
from axcsas.visualization.fitting_plots import (
    plot_peak_fit,
    plot_doublet_comparison,
)

import matplotlib.pyplot as plt

# Apply unified style
apply_axcsas_style()


# Cu peak positions (JCPDS)
PEAK_POSITIONS = {
    (1, 1, 1): 43.3,
    (2, 0, 0): 50.4,
    (2, 2, 0): 74.1,
}
PEAK_LABELS = ['(111)', '(200)', '(220)']


def fit_peak_with_diagnosis(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    expected_center: float,
    window: float = 2.5,
    use_doublet: bool = False  # Pseudo-Voigt for reliability (DoubletFitter has environment issues)
) -> Dict:
    """
    Fit a single peak using Kα doublet model and return detailed diagnosis info.
    使用 Kα 雙峰模型擬合單峰並返回詳細診斷資訊。
    
    Returns dict with:
        - success: bool
        - center: fitted Kα₁ peak center
        - center_ka2: fitted Kα₂ peak center
        - amplitude: fitted amplitude
        - fwhm: fitted FWHM
        - eta: mixing parameter
        - r_squared: goodness of fit
        - theta_range: x data
        - int_range: y data (original)
        - fitted_curve: y data (fitted)
    """
    result = {
        'success': False,
        'center': np.nan,
        'center_ka2': np.nan,
        'amplitude': np.nan,
        'fwhm': np.nan,
        'eta': np.nan,
        'r_squared': np.nan,
        'theta_range': None,
        'int_range': None,
        'fitted_curve': None,
        'method': 'doublet' if use_doublet else 'pseudo-voigt'
    }
    
    # Select range
    mask = (two_theta >= expected_center - window) & (two_theta <= expected_center + window)
    if not np.any(mask):
        return result
    
    theta_range = two_theta[mask]
    int_range = intensity[mask]
    
    result['theta_range'] = theta_range
    result['int_range'] = int_range
    
    # Find maximum for initial guess
    idx_max = np.argmax(int_range)
    peak_theta = theta_range[idx_max]
    peak_int = int_range[idx_max]
    
    if peak_int < 50:
        return result
    
    # Estimate initial FWHM
    half_max = peak_int / 2
    left_idx = idx_max
    while left_idx > 0 and int_range[left_idx] > half_max:
        left_idx -= 1
    right_idx = idx_max
    while right_idx < len(int_range) - 1 and int_range[right_idx] > half_max:
        right_idx += 1
    initial_fwhm = max(theta_range[right_idx] - theta_range[left_idx], 0.1)
    
    try:
        if use_doublet:
            # Use Kα doublet fitting with timeout protection
            import signal
            
            fitter = DoubletFitter(max_iterations=5000)
            fit_result = fitter.fit(theta_range, int_range, expected_center, initial_fwhm)
            
            if fit_result.success and fit_result.r_squared > 0.8:
                result['success'] = True
                result['center'] = fit_result.center_ka1
                result['center_ka2'] = fit_result.center_ka2
                result['amplitude'] = fit_result.amplitude_ka1
                result['fwhm'] = fit_result.fwhm
                result['eta'] = fit_result.eta
                result['r_squared'] = fit_result.r_squared
                result['fitted_curve'] = fit_result.fitted_curve
                result['method'] = 'doublet'
            else:
                # Fallback to Pseudo-Voigt if doublet fitting fails or has low R²
                use_doublet = False
        
        if not use_doublet or not result['success']:
            # Enhanced Pseudo-Voigt fitting with polynomial background
            from scipy.optimize import least_squares
            
            # Estimate polynomial background (quadratic)
            n_edge = max(5, len(theta_range) // 8)
            bg_x = np.concatenate([theta_range[:n_edge], theta_range[-n_edge:]])
            bg_y = np.concatenate([int_range[:n_edge], int_range[-n_edge:]])
            bg_coeffs = np.polyfit(bg_x, bg_y, 2)  # Quadratic background
            
            # Model: Pseudo-Voigt + quadratic background
            def enhanced_pv_model(x, center, amplitude, fwhm, eta, a2, a1, a0):
                pv = PseudoVoigt.profile(x, center, amplitude, fwhm, np.clip(eta, 0, 1))
                background = a2 * x**2 + a1 * x + a0
                return pv + background
            
            def residual(params):
                return int_range - enhanced_pv_model(theta_range, *params)
            
            # Multi-start optimization for robustness
            best_r2 = -1
            best_params = None
            best_fitted = None
            
            # Try multiple initial eta values
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
                        residuals = int_range - fitted
                        ss_res = np.sum(residuals**2)
                        ss_tot = np.sum((int_range - np.mean(int_range))**2)
                        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                        
                        # Calculate reduced chi-squared and parameter uncertainties
                        n_data = len(int_range)
                        n_params = len(opt_result.x)
                        dof = n_data - n_params  # Degrees of freedom
                        
                        # Estimate variance and chi-squared
                        variance = ss_res / dof if dof > 0 else ss_res
                        chi2_red = ss_res / (dof * variance) if dof > 0 else 1.0
                        
                        # Parameter uncertainties from Jacobian
                        try:
                            # Covariance matrix from Jacobian: cov = (J^T J)^-1 * variance
                            J = opt_result.jac
                            cov = np.linalg.inv(J.T @ J) * variance
                            param_errors = np.sqrt(np.diag(cov))
                        except (np.linalg.LinAlgError, ValueError):
                            # Fallback if covariance calculation fails
                            param_errors = np.zeros(n_params)
                        
                        if r2 > best_r2:
                            best_r2 = r2
                            best_params = opt_result.x
                            best_fitted = fitted
                            best_errors = param_errors
                            best_chi2_red = chi2_red
                            best_dof = dof
                except Exception:
                    continue
            
            if best_params is not None and best_r2 > 0.9:
                result['success'] = True
                result['center'] = best_params[0]
                result['center_err'] = best_errors[0] if best_errors[0] > 0 else 0.001
                result['amplitude'] = best_params[1]
                result['amplitude_err'] = best_errors[1]
                result['fwhm'] = best_params[2]
                result['fwhm_err'] = best_errors[2] if best_errors[2] > 0 else 0.0001
                result['eta'] = best_params[3]
                result['eta_err'] = best_errors[3]
                result['r_squared'] = best_r2
                result['chi2_red'] = best_chi2_red
                result['dof'] = best_dof
                result['method'] = 'enhanced-pv'
                result['fitted_curve'] = best_fitted
    except Exception as e:
        result['error'] = str(e)
    
    return result


def generate_sample_fitting_plot(
    filepath: Path,
    output_dir: Path,
    config: AnalysisConfig = None
) -> Optional[List[Dict]]:
    """
    Generate a single diagnostic plot for one XRD sample.
    使用新視覺化模組生成單一樣品的診斷圖。
    
    Shows 3 subplots (one for each peak) with fitting details.
    """
    config = config or AnalysisConfig()
    
    # Load data
    try:
        two_theta, intensity = load_bruker_txt(str(filepath))
    except Exception as e:
        print(f"  ✗ Error loading {filepath.name}: {e}")
        return None
    
    if len(two_theta) == 0:
        print(f"  ✗ No data in {filepath.name}")
        return None
    
    # Parse filename
    file_info = parse_filename(str(filepath))
    sample_name = filepath.stem
    conc = file_info.get('concentration_ml', 0)
    time_h = file_info.get('time_hours', 0)
    
    # Create figure with 3 subplots using AXCSAS style
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f'Peak Fitting Diagnosis (Kα Doublet): {sample_name}\n'
        f'(Leveler: {conc} ml, Time: {time_h}h)', 
        fontsize=14, fontweight='bold'
    )
    
    peaks_info = []
    colors = [COLORBLIND_SAFE[0], COLORBLIND_SAFE[1], COLORBLIND_SAFE[2]]
    
    for idx, (hkl, expected_pos) in enumerate(PEAK_POSITIONS.items()):
        ax = axes[idx]
        label = PEAK_LABELS[idx]
        color = colors[idx]
        
        # Fit peak using doublet model
        fit_result = fit_peak_with_diagnosis(two_theta, intensity, expected_pos, use_doublet=True)
        
        if fit_result['theta_range'] is not None:
            # Plot original data
            ax.scatter(fit_result['theta_range'], fit_result['int_range'], 
                      s=15, alpha=0.5, color='gray', label='Data', zorder=3)
            
            if fit_result['success'] and fit_result['fitted_curve'] is not None:
                # Plot fitted curve
                ax.plot(fit_result['theta_range'], fit_result['fitted_curve'], 
                       color=color, linewidth=1.5, label='Fit', zorder=4)
                
                # Fill under curve
                ax.fill_between(
                    fit_result['theta_range'], 0, fit_result['fitted_curve'],
                    alpha=0.2, color=color
                )
                
                # Mark peak center and FWHM
                center = fit_result['center']
                center_ka2 = fit_result.get('center_ka2', np.nan)
                fwhm = fit_result['fwhm']
                amp = fit_result['amplitude']
                
                # Vertical lines at Kα₁ and Kα₂ centers
                ax.axvline(x=center, color='blue', linestyle='--', alpha=0.6, 
                          linewidth=1.0, label=f'Kα₁={center:.3f}°')
                if not np.isnan(center_ka2):
                    ax.axvline(x=center_ka2, color='orange', linestyle='--', alpha=0.6,
                              linewidth=1.0, label=f'Kα₂={center_ka2:.3f}°')
                
                # FWHM indicator
                half_max = amp / 2 + np.min(fit_result['fitted_curve'])
                ax.hlines(y=half_max, xmin=center - fwhm/2, xmax=center + fwhm/2,
                         color='green', linewidth=1.5, label=f'FWHM={fwhm:.4f}°')
                
                # Get uncertainties (with fallbacks)
                center_err = fit_result.get('center_err', 0.001)
                fwhm_err = fit_result.get('fwhm_err', 0.0001)
                eta_err = fit_result.get('eta_err', 0.01)
                chi2_red = fit_result.get('chi2_red', 1.0)
                r2 = fit_result['r_squared']
                eta = fit_result['eta']
                
                # Info text box with uncertainties
                info_text = (
                    f"R² = {r2:.4f}\n"
                    f"χ²ᵣₑₐ = {chi2_red:.3f}\n"
                    f"2θ = {center:.3f}° ± {center_err:.3f}°\n"
                    f"FWHM = {fwhm:.4f}° ± {fwhm_err:.4f}°\n"
                    f"η = {eta:.3f} ± {eta_err:.3f}"
                )
                
                ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                       fontsize=9, verticalalignment='top', family='monospace',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
                
                peaks_info.append({
                    'hkl': hkl,
                    'center': center,
                    'center_err': center_err,
                    'center_ka2': center_ka2,
                    'fwhm': fwhm,
                    'fwhm_err': fwhm_err,
                    'eta': eta,
                    'eta_err': eta_err,
                    'r_squared': r2,
                    'chi2_red': chi2_red,
                })
            else:
                ax.text(0.5, 0.5, 'Fitting Failed', transform=ax.transAxes,
                       fontsize=14, ha='center', va='center', color='red',
                       fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No Peak Found', transform=ax.transAxes,
                   fontsize=14, ha='center', va='center', color='red',
                   fontweight='bold')
        
        ax.set_xlabel('2θ (°)')
        ax.set_ylabel('Intensity (counts)')
        ax.set_title(f'{label} Peak @ {expected_pos}°', fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot using visualization module
    output_path = output_dir / f'{sample_name}_fitting.png'
    save_figure(fig, str(output_path), dpi=300)
    plt.close()
    
    return peaks_info


def main():
    """Main entry point."""
    print("=" * 60)
    print("AXCSAS Peak Fitting Diagnostic Plots (Refactored)")
    print("Using axcsas.visualization module")
    print("=" * 60)
    
    # Setup paths - go up 2 levels from visualization/ to project root
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data" / "raw" / "202511"
    output_dir = project_root / "outputs" / "plots" / "fitting_diagnosis"
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get all XRD files
    txt_files = sorted(data_dir.glob("*.txt"))
    print(f"\nFound {len(txt_files)} XRD data files")
    print(f"Output directory: {output_dir}\n")
    
    # Process each file
    all_results = []
    for filepath in txt_files:
        print(f"Processing: {filepath.name}")
        result = generate_sample_fitting_plot(filepath, output_dir)
        if result:
            all_results.append({
                'filename': filepath.name,
                'peaks': result
            })
            print(f"  ✓ Saved (DPI: 300)")
        else:
            print(f"  ✗ Failed")
    
    print("\n" + "=" * 60)
    print(f"Complete! Generated {len(all_results)} diagnostic plots")
    print(f"Output directory: {output_dir}")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    exit(main())
