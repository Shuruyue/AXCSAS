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


# Cu peak positions (JCPDS) - use unified function
from axcsas.core.copper_crystal import get_standard_peaks

PEAK_POSITIONS = get_standard_peaks()  # Returns (111), (200), (220)
PEAK_LABELS = [f"({h[0]}{h[1]}{h[2]})" for h in PEAK_POSITIONS.keys()]


from axcsas.fitting.peak_fitter import fit_peak_with_diagnosis


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
        f'Peak Fitting Diagnosis: {sample_name}\n'
        f'(Leveler: {conc:.1f} mL/1.5L, Annealing Time: {time_h:.0f}h)', 
        fontsize=14, fontweight='bold'
    )
    
    peaks_info = []
    colors = [COLORBLIND_SAFE[0], COLORBLIND_SAFE[1], COLORBLIND_SAFE[2]]
    
    for idx, (hkl, expected_pos) in enumerate(PEAK_POSITIONS.items()):
        ax = axes[idx]
        label = PEAK_LABELS[idx]
        color = colors[idx]
        
        # Fit peak using doublet model
        # Use DoubletFitter first (higher R²), with Enhanced PV as fallback
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
                
                # Vertical lines at Kα₁ and Kα₂ centers (gray, different dash styles)
                ax.axvline(x=center, color='gray', linestyle='--', alpha=0.7, 
                          linewidth=1.0, label=f'Kα₁={center:.3f}°')
                if not np.isnan(center_ka2):
                    ax.axvline(x=center_ka2, color='gray', linestyle=':', alpha=0.7,
                              linewidth=1.0, label=f'Kα₂={center_ka2:.3f}°')
                    
                    # Plot separate Kα₁ and Kα₂ components if doublet method used
                    if fit_result.get('method') in ['doublet', 'doublet-true-voigt']:
                         from axcsas.fitting.pseudo_voigt import TrueVoigt
                         
                         eta_val = fit_result['eta']
                         amp_ka1 = fit_result['amplitude']
                         amp_ka2 = amp_ka1 * 0.5 # Approximation
                         fwhm_val = fit_result['fwhm']
                         
                         # Get X range
                         x_vals = fit_result['theta_range']
                         
                         # Recover True Voigt params from FWHM/Eta
                         # This allows plotting consistent with the reported metrics
                         # Note: params_from_fwhm is an estimation but consistent with the fitting pipeline
                         sigma, gamma = TrueVoigt.params_from_fwhm(fwhm_val, eta_val)
                         
                         # Calculate profiles using True Voigt
                         y_ka1 = TrueVoigt.profile(x_vals, center, amp_ka1, sigma, gamma)
                         y_ka2 = TrueVoigt.profile(x_vals, center_ka2, amp_ka2, sigma, gamma)
                         
                         # Estimate background from fitted_curve - (ka1 + ka2)
                         y_total = fit_result['fitted_curve']
                         y_bg = y_total - (y_ka1 + y_ka2)
                         
                         # Plot components with background added (stacked-ish)
                         # User requested uniform gray style, and implicit legend association
                         # Both set to nolegend to avoid redundancy with the vertical line labels
                         # Correction: Match linestyle to vertical lines (Ka1=dashed, Ka2=dotted)
                         ax.plot(x_vals, y_ka1 + y_bg, color='gray', linestyle='--', linewidth=1.0, alpha=0.5, label='_nolegend_')
                         ax.plot(x_vals, y_ka2 + y_bg, color='gray', linestyle=':', linewidth=1.0, alpha=0.5, label='_nolegend_')
                
                # FWHM indicator (muted red)
                intensity_at_half_max = amp / 2 + np.min(fit_result['fitted_curve'])
                ax.hlines(y=intensity_at_half_max, xmin=center - fwhm/2, xmax=center + fwhm/2,
                          color='#C44E52', linewidth=1.5, label=f'FWHM={fwhm:.4f}°')
                # Get uncertainties - NO misleading defaults! Use NaN if calculation failed
                center_err = fit_result.get('center_err', np.nan)
                fwhm_err = fit_result.get('fwhm_err', np.nan)
                eta_err = fit_result.get('eta_err', np.nan)
                chi2_red = fit_result.get('chi2_red', np.nan)  # ← Changed from 1.0
                r2 = fit_result['r_squared']
                eta = fit_result['eta']
                low_quality = fit_result.get('low_quality', False)
                
                # Info text box with uncertainties
                # Add prominent warning if low quality fit
                if low_quality:
                    quality_warning = "⚠️ R² < 0.995\n"
                else:
                    quality_warning = ""
                
                info_text = (
                    f"{quality_warning}"
                    f"R² = {r2:.4f}\n"
                    f"2θ = {center:.3f}° ± {center_err:.3f}°\n"
                    f"FWHM = {fwhm:.4f}° ± {fwhm_err:.4f}°\n"
                    f"η = {eta:.3f} ± {eta_err:.3f}"
                )
                
                # Use orange box if low quality
                box_color = 'orange' if low_quality else 'white'
                edge_color = 'red' if low_quality else 'gray'
                ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                       fontsize=9, verticalalignment='top', family='serif',
                       bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.9, edgecolor=edge_color))
                
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
    save_figure(fig, str(output_path), dpi=1200)
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
            print(f"  [OK] Saved (DPI: 1000)")
        else:
            print(f"  ✗ Failed")
    
    print("\n" + "=" * 60)
    print(f"Complete! Generated {len(all_results)} diagnostic plots")
    print(f"Output directory: {output_dir}")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    exit(main())
