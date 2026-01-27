#!/usr/bin/env python3
"""
Fitting Method Comparison Script
=================================

Compare three XRD peak fitting methods (from simple to accurate):
1. Single Pseudo-Voigt - Ignores Kα₁/Kα₂ splitting entirely
2. Kα₂ Stripping - Remove Kα₂ first (Rachinger method), then fit single peak
3. Kα Doublet Fitting - Fit Kα₁ + Kα₂ simultaneously (most accurate)

Style unified with fitting_diagnosis: Times New Roman, gray Kα lines, FWHM markers.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from axcsas.analysis.pipeline import load_bruker_txt, parse_filename
from axcsas.fitting.pseudo_voigt import PseudoVoigt
from axcsas.fitting.lm_optimizer import LMOptimizer
from axcsas.fitting.ka_doublet import DoubletFitter, Ka2Stripper
from axcsas.core.copper_crystal import get_standard_peaks

# Import AXCSAS style (Times New Roman)
from axcsas.visualization.style import (
    apply_axcsas_style,
    COLORBLIND_SAFE,
    save_figure,
)

# Apply unified style
apply_axcsas_style()

# Peak configuration
PEAKS = get_standard_peaks()
PEAK_LABELS = [f"({h[0]}{h[1]}{h[2]})" for h in PEAKS.keys()]


def compare_methods_for_sample(filepath: Path, output_dir: Path):
    """Generate comparison plot for one sample (unified fitting_diagnosis style)."""
    
    # Load data
    try:
        two_theta, intensity = load_bruker_txt(str(filepath))
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None
    
    file_info = parse_filename(str(filepath))
    sample_name = filepath.stem
    conc = file_info.get('concentration_ml', 0)
    time_h = file_info.get('time_hours', 0)
    
    # Create figure: 3 rows (methods) × 3 columns (peaks)
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), facecolor='white')
    fig.suptitle(
        f'Fitting Method Comparison: {sample_name}\n'
        f'Leveler: {conc:.1f} mL/1.5L, Annealing Time: {time_h:.1f}h', 
        fontsize=12, fontweight='normal'
    )
    
    method_names = [
        "Method 1: Single Pseudo-Voigt",
        "Method 2: Kα₂ Stripping",
        "Method 3: Doublet Fitting"
    ]
    
    colors = [COLORBLIND_SAFE[0], COLORBLIND_SAFE[1], COLORBLIND_SAFE[2]]  # 藍、橙、青
    
    for col, ((hkl, expected_pos), label) in enumerate(zip(PEAKS.items(), PEAK_LABELS)):
        window = 2.5
        mask = (two_theta >= expected_pos - window) & (two_theta <= expected_pos + window)
        theta_region = two_theta[mask]
        int_region = intensity[mask]
        
        if len(theta_region) < 10:
            continue
        
        color = colors[col]
        
        # ==== Row 0: Method 1 - Single Pseudo-Voigt ====
        ax1 = axes[0, col]
        optimizer = LMOptimizer()
        result1 = optimizer.fit_single_peak(theta_region, int_region)
        
        ax1.scatter(theta_region, int_region, s=15, alpha=0.5, color='gray', label='Data', zorder=3)
        if result1.success:
            bg = np.min(int_region)
            fitted1 = PseudoVoigt.profile(
                theta_region, result1.params.center, 
                result1.params.amplitude, result1.params.fwhm, 
                result1.params.eta
            ) + bg
            ax1.plot(theta_region, fitted1, color=color, lw=1.5, label='Fit', zorder=4)
            ax1.fill_between(theta_region, 0, fitted1, alpha=0.15, color=color)
            
            # 2θ position line (gray dashed)
            fwhm = result1.params.fwhm
            center = result1.params.center
            ax1.axvline(center, color='gray', linestyle='--', alpha=0.7, lw=1,
                       label=f'2θ={center:.3f}°')
            
            # FWHM marker
            half_max = result1.params.amplitude / 2 + bg
            ax1.hlines(y=half_max, xmin=center - fwhm/2, xmax=center + fwhm/2,
                      color='#C44E52', linewidth=1.5, label=f'FWHM={fwhm:.4f}°')
            
            # Info text box - unified serif font
            info_text = (
                f"R² = {result1.r_squared:.4f}\n"
                f"2θ = {result1.params.center:.3f}°\n"
                f"FWHM = {result1.params.fwhm:.4f}°\n"
                f"η = {result1.params.eta:.3f}"
            )
            ax1.text(0.02, 0.98, info_text, transform=ax1.transAxes,
                    fontsize=9, verticalalignment='top', family='serif',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
        
        ax1.set_xlabel('2θ (°)')
        ax1.set_ylabel('Intensity (counts)')
        ax1.set_title(f'{label} Peak @ {expected_pos}°', fontsize=10, fontweight='bold')
        ax1.legend(loc='upper right', fontsize=7)
        ax1.grid(True, alpha=0.3, linestyle=':')
        
        # ==== Row 1: Method 2 - Kα₂ Stripping ====
        ax2 = axes[1, col]
        
        stripper = Ka2Stripper()
        theta_stripped, int_stripped = stripper.strip_peak_region(
            two_theta, intensity, expected_pos, window=window
        )
        
        result2 = optimizer.fit_single_peak(theta_stripped, int_stripped)
        
        # Gray dots for both original and stripped
        ax2.scatter(theta_region, int_region, s=10, alpha=0.3, color='lightgray', label='Original', zorder=2)
        ax2.scatter(theta_stripped, int_stripped, s=15, alpha=0.5, color='gray', label='Stripped', zorder=3)
        
        if result2.success:
            bg2 = np.min(int_stripped)
            fitted2 = PseudoVoigt.profile(
                theta_stripped, result2.params.center,
                result2.params.amplitude, result2.params.fwhm,
                result2.params.eta
            ) + bg2
            ax2.plot(theta_stripped, fitted2, color=color, lw=1.5, label='Fit', zorder=4)
            ax2.fill_between(theta_stripped, 0, fitted2, alpha=0.15, color=color)
            
            # Gray vertical line
            ax2.axvline(result2.params.center, color='gray', linestyle='--', alpha=0.7, lw=1)
            
            # FWHM marker
            fwhm = result2.params.fwhm
            center = result2.params.center
            half_max = result2.params.amplitude / 2 + bg2
            ax2.hlines(y=half_max, xmin=center - fwhm/2, xmax=center + fwhm/2,
                      color='#C44E52', linewidth=1.5, label=f'FWHM={fwhm:.4f}°')
            
            # Info text box
            info_text = (
                f"R² = {result2.r_squared:.4f}\n"
                f"2θ = {result2.params.center:.3f}°\n"
                f"FWHM = {result2.params.fwhm:.4f}°\n"
                f"η = {result2.params.eta:.3f}"
            )
            ax2.text(0.02, 0.98, info_text, transform=ax2.transAxes,
                    fontsize=9, verticalalignment='top', family='serif',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
        
        ax2.set_xlabel('2θ (°)')
        ax2.set_ylabel('Intensity (counts)')
        ax2.set_title(f'{label} Peak @ {expected_pos}°', fontsize=10, fontweight='bold')
        ax2.legend(loc='upper right', fontsize=7)
        ax2.grid(True, alpha=0.3, linestyle=':')
        
        # ==== Row 2: Method 3 - Kα Doublet Fitting ====
        ax3 = axes[2, col]
        
        doublet_fitter = DoubletFitter(max_iterations=100000)
        result3 = doublet_fitter.fit(theta_region, int_region, expected_pos)
        
        ax3.scatter(theta_region, int_region, s=15, alpha=0.5, color='gray', label='Data', zorder=3)
        if result3.success and result3.fitted_curve is not None:
            ax3.plot(theta_region, result3.fitted_curve, color=color, lw=1.5, label='Fit', zorder=4)
            ax3.fill_between(theta_region, 0, result3.fitted_curve, alpha=0.15, color=color)
            
            # Gray Kα₁ and Kα₂ individual peak curves (True Voigt)
            from axcsas.fitting.pseudo_voigt import TrueVoigt
            bg3 = np.min(result3.fitted_curve)
            
            # Recover True Voigt params
            sigma, gamma = TrueVoigt.params_from_fwhm(result3.fwhm, result3.eta)
            
            ka1_curve = TrueVoigt.profile(
                theta_region, result3.center_ka1,
                result3.amplitude_ka1, sigma, gamma
            ) + bg3
            
            ka2_curve = TrueVoigt.profile(
                theta_region, result3.center_ka2,
                result3.amplitude_ka1 * 0.5, sigma, gamma
            ) + bg3
            
            ax3.plot(theta_region, ka1_curve, color='gray', linestyle='--', lw=1, alpha=0.7, label='Kα₁ peak')
            ax3.plot(theta_region, ka2_curve, color='gray', linestyle=':', lw=1, alpha=0.7, label='Kα₂ peak')
            
            # Gray vertical lines for Kα₁ and Kα₂ positions
            ax3.axvline(result3.center_ka1, color='gray', linestyle='--', alpha=0.5, lw=0.8)
            ax3.axvline(result3.center_ka2, color='gray', linestyle=':', alpha=0.5, lw=0.8)
            
            # FWHM marker
            fwhm = result3.fwhm
            center = result3.center_ka1
            half_max = result3.amplitude_ka1 / 2 + bg3
            ax3.hlines(y=half_max, xmin=center - fwhm/2, xmax=center + fwhm/2,
                      color='#C44E52', linewidth=1.5, label=f'FWHM={fwhm:.4f}°')
            
            # Info text box
            info_text = (
                f"R² = {result3.r_squared:.4f}\n"
                f"Kα₁ = {result3.center_ka1:.3f}°\n"
                f"Kα₂ = {result3.center_ka2:.3f}°\n"
                f"FWHM = {result3.fwhm:.4f}°\n"
                f"η = {result3.eta:.3f}"
            )
            ax3.text(0.02, 0.98, info_text, transform=ax3.transAxes,
                    fontsize=9, verticalalignment='top', family='serif',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
        
        ax3.set_xlabel('2θ (°)')
        ax3.set_ylabel('Intensity (counts)')
        ax3.set_title(f'{label} Peak @ {expected_pos}°', fontsize=10, fontweight='bold')
        ax3.legend(loc='upper right', fontsize=7)
        ax3.grid(True, alpha=0.3, linestyle=':')
    
    # Add row labels
    for row, method in enumerate(method_names):
        axes[row, 0].text(-0.35, 0.5, method, transform=axes[row, 0].transAxes,
                          fontsize=11, fontweight='bold', rotation=90,
                          ha='center', va='center', family='serif')
    
    plt.tight_layout(rect=[0.05, 0, 1, 0.97])
    
    # Save using AXCSAS save_figure
    output_path = output_dir / f'{sample_name}_method_comparison.png'
    save_figure(fig, str(output_path), dpi=1200)
    plt.close()
    
    return True


def main():
    print("=" * 60)
    print("XRD Peak Fitting Method Comparison")
    print("=" * 60)
    print("Method 1: Single Pseudo-Voigt (ignores Kα splitting)")
    print("Method 2: Kα₂ Stripping (Rachinger method)")
    print("Method 3: Kα Doublet Fitting (simultaneous Kα₁ + Kα₂)")
    print("=" * 60)
    
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data" / "raw" / "202511"
    if not data_dir.exists():
        data_dir = project_root / "data" / "raw"
    
    output_dir = project_root / "outputs" / "plots" / "compare_integral_methods"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Cleanup old plots
    for f in output_dir.glob("*.png"):
        f.unlink()
    print(f"Cleaned up old comparison plots in {output_dir}")
    
    txt_files = sorted(data_dir.glob("*.txt"))
    
    print(f"\nProcessing {len(txt_files)} samples...")
    print(f"Output: {output_dir}\n")
    
    success_count = 0
    for filepath in txt_files:
        print(f"Processing: {filepath.name}")
        result = compare_methods_for_sample(filepath, output_dir)
        if result:
            success_count += 1
            print(f"  ✓ Saved")
    
    print("\n" + "=" * 60)
    print(f"Complete! Generated {success_count} comparison plots")
    print("=" * 60)


if __name__ == "__main__":
    main()
