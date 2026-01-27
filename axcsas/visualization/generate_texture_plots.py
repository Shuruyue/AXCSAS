"""
Texture Analysis Plot Generator
================================

Generates Harris Texture Coefficient (TC) analysis plots for all samples.
使用 Harris 紋理係數分析所有樣品並生成 TC 演化圖。
"""

import os
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from axcsas.analysis.pipeline import load_bruker_txt
from axcsas.fitting.ka_doublet import DoubletFitter
from axcsas.fitting.pv_area import calculate_pv_area  # Correct PV area calculation
from axcsas.methods.texture import TextureAnalyzer, analyze_texture
from axcsas.visualization.texture_plots import plot_texture_polar, plot_tc_evolution
from axcsas.visualization.style import save_figure
from axcsas.visualization.generate_fitting_diagnosis import fit_peak_with_diagnosis

import matplotlib.pyplot as plt
import numpy as np
import re


# Peak positions for Cu (111), (200), (220)
# Use unified functions from core and analysis
from axcsas.core.copper_crystal import get_standard_peaks
from axcsas.analysis.pipeline import parse_filename

# Peak positions for Cu (111), (200), (220) - from unified source
PEAK_POSITIONS = get_standard_peaks()

# Quality threshold for texture analysis
# R² threshold for high-quality texture coefficient calculation
# 織構係數計算的高品質 R² 閾值
R2_THRESHOLD = 0.95  # User-defined quality standard / 使用者自定義品質標準


# parse_sample_info removed - use parse_filename() from pipeline
# parse_filename returns {'concentration_ml': float, 'time_hours': float, 'name': str}
# Usage: info = parse_filename(filename); info['concentration_ml'], info['time_hours']



def fit_peaks_and_get_intensities(two_theta: np.ndarray, intensity: np.ndarray, verbose: bool = False) -> tuple:
    """Fit ALL peaks using DoubletFitter and calculate accurate integrated areas.
    
    Returns:
        (intensities_dict, all_passed, r2_values, area_values)
        - intensities_dict: Dict of (hkl) -> integrated intensity (only R² ≥ 0.95)
        - all_passed: True if all 3 peaks have R² ≥ 0.95
        - r2_values: Dict of (hkl) -> R² value
        - area_values: Dict of (hkl) -> fitted area (all peaks)
    
    Note:
        Uses calculate_pv_area() with actual η from fitting for accurate area calculation.
        This is critical for texture coefficient accuracy!
    """
    intensities = {}
    r2_values = {}
    area_values = {}
    
    for hkl, expected_pos in PEAK_POSITIONS.items():
        try:
            # Use DoubletFitter (same as fitting_diagnosis)
            result = fit_peak_with_diagnosis(
                two_theta, 
                intensity, 
                expected_pos, 
                window=2.5, 
                use_doublet=True
            )
            
            if result['success']:
                amplitude = result['amplitude']
                fwhm = result['fwhm']
                eta = result['eta']  # ← CRITICAL: Use actual fitted η
                r2 = result['r_squared']
                
                # ✓ CORRECT: Use dynamic η-based area calculation
                peak_area = calculate_pv_area(amplitude, fwhm, eta)
                
                r2_values[hkl] = r2
                area_values[hkl] = peak_area
                
                if r2 >= R2_THRESHOLD and peak_area > 0:
                    intensities[hkl] = peak_area
                        
        except Exception as e:
            if verbose:
                print(f"  Error fitting {hkl}: {e}")
            continue
    
    # Check if ALL 3 peaks passed
    all_passed = len(intensities) == 3 and all(r2 >= R2_THRESHOLD for r2 in r2_values.values())
    
    return intensities, all_passed, r2_values, area_values


def plot_texture_diagnosis(
    two_theta: np.ndarray, 
    intensity: np.ndarray, 
    sample_name: str,
    output_path: str
) -> None:
    """Generate texture fitting diagnosis plot matching fitting_diagnosis format.
    
    Creates a 3-panel figure (same format as fitting_diagnosis) with:
    - Gray scatter data points
    - Colored fit line with fill
    - Peak center vertical lines
    - FWHM indicator
    - Info box with R², Amp, FWHM, η, Area (with correct calculation)
    """
    from axcsas.visualization.style import apply_axcsas_style, COLORBLIND_SAFE
    apply_axcsas_style()
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    PEAK_LABELS = ['(111)', '(200)', '(220)']
    colors = [COLORBLIND_SAFE[0], COLORBLIND_SAFE[1], COLORBLIND_SAFE[2]]
    
    for idx, (hkl, expected_pos) in enumerate(PEAK_POSITIONS.items()):
        ax = axes[idx]
        hkl_str = PEAK_LABELS[idx]
        color = colors[idx]
        
        try:
            result = fit_peak_with_diagnosis(
                two_theta, intensity, expected_pos, window=2.5, use_doublet=True
            )
            
            theta_range = result.get('theta_range')
            int_range = result.get('int_range')
            fitted_curve = result.get('fitted_curve')
            
            if theta_range is not None:
                # Plot original data (gray scatter like fitting_diagnosis)
                ax.scatter(theta_range, int_range, s=15, alpha=0.5, 
                          color='gray', label='Data', zorder=3)
                
                if result['success'] and fitted_curve is not None:
                    # Plot fitted curve (colored line)
                    ax.plot(theta_range, fitted_curve, color=color, 
                           linewidth=1.5, label='Fit', zorder=4)
                    
                    # Fill under curve
                    ax.fill_between(theta_range, 0, fitted_curve,
                                   alpha=0.2, color=color)
                    
                    r2 = result['r_squared']
                    amp = result['amplitude']
                    fwhm = result['fwhm']
                    center = result['center']
                    center_ka2 = result.get('center_ka2', np.nan)
                    eta = result.get('eta', np.nan)
                    
                    # ✓ CORRECT: Use dynamic η-based area calculation
                    area = calculate_pv_area(amp, fwhm, eta)
                    
                    # Get uncertainties (matching fitting_diagnosis)
                    center_err = result.get('center_err', np.nan)
                    fwhm_err = result.get('fwhm_err', np.nan)
                    eta_err = result.get('eta_err', np.nan)
                    
                    # Vertical lines at Kα₁ and Kα₂ centers (gray, matching fitting_diagnosis)
                    ax.axvline(x=center, color='gray', linestyle='--', alpha=0.7,
                              linewidth=1.0, label=f'Kα₁={center:.3f}°')
                    if not np.isnan(center_ka2):
                        ax.axvline(x=center_ka2, color='gray', linestyle=':', alpha=0.7,
                                  linewidth=1.0, label=f'Kα₂={center_ka2:.3f}°')
                    
                    # FWHM indicator (muted red, matching fitting_diagnosis)
                    half_max = amp / 2 + np.min(fitted_curve)
                    ax.hlines(y=half_max, xmin=center - fwhm/2, xmax=center + fwhm/2,
                             color='#C44E52', linewidth=1.5, label=f'FWHM={fwhm:.4f}°')
                    
                    # Info text box (matching fitting_diagnosis format exactly)
                    info_text = (
                        f"R² = {r2:.4f}\n"
                        f"2θ = {center:.3f}° ± {center_err:.3f}°\n"
                        f"FWHM = {fwhm:.4f}° ± {fwhm_err:.4f}°\n"
                        f"η = {eta:.3f} ± {eta_err:.3f}\n"
                        f"Area = {area:.1f}"
                    )
                    
                    # Always use white box (no warning for texture plots)
                    ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                           fontsize=9, verticalalignment='top', family='serif',
                           bbox=dict(boxstyle='round', facecolor='white', 
                                    edgecolor='gray', alpha=0.9))
                
        except Exception as e:
            ax.text(0.5, 0.5, f'Error: {e}', transform=ax.transAxes,
                   ha='center', va='center', color='red')
        
        ax.set_xlabel('2θ (degrees)')
        ax.set_ylabel('Intensity (counts)')
        ax.set_title(f'{hkl_str} Peak', fontweight='bold')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    fig.suptitle(f'Texture Fitting Diagnosis: {sample_name}', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Use save_figure for consistent DPI (2400)
    from axcsas.visualization.style import save_figure
    save_figure(plt.gcf(), output_path)
    plt.close(fig)


def main():
    print("=" * 60)
    print("AXCSAS Texture Analysis Plot Generator")
    print("=" * 60)
    
    # Setup paths
    data_dir = project_root / "data" / "raw" / "202511"
    output_dir = project_root / "outputs" / "plots" / "texture"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Diagnosis plots output directory
    diagnosis_dir = project_root / "outputs" / "plots" / "texture_diagnosis"
    diagnosis_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nData directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Diagnosis directory: {diagnosis_dir}")
    
    # Find all data files
    data_files = sorted(data_dir.glob("*.txt"))
    print(f"\nFound {len(data_files)} data files")
    
    if not data_files:
        print("ERROR: No data files found!")
        return
    
    # Initialize texture analyzer
    analyzer = TextureAnalyzer(use_area=True)
    
    # Process all samples
    all_tc_data = []
    
    for filepath in data_files:
        filename = filepath.name
        print(f"Processing: {filename}")
        
        try:
            two_theta, intensity = load_bruker_txt(str(filepath))
            sample_info = parse_filename(str(filepath))
            # Map keys: parse_filename uses 'concentration_ml' and 'time_hours'
            # Convert to match existing usage
            sample_info = {
                'concentration': sample_info.get('concentration_ml', 0),
                'time': sample_info.get('time_hours', 0),
            }
            
            # Generate texture diagnosis plot
            sample_basename = filename.replace('.txt', '')
            diagnosis_path = diagnosis_dir / f"tc_diagnosis_{sample_basename}.png"
            plot_texture_diagnosis(two_theta, intensity, sample_basename, str(diagnosis_path))
            
            # Fit peaks and get intensities (requires R² ≥ 0.995 for all 3 peaks)
            intensities, all_passed, r2_values, area_values = fit_peaks_and_get_intensities(two_theta, intensity)
            
            # Display R² and fitted areas for diagnostics
            r2_strs = []
            area_strs = []
            for hkl in PEAK_POSITIONS.keys():
                hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
                if hkl in r2_values:
                    r2 = r2_values[hkl]
                    status = "✓" if r2 >= 0.995 else "✗"
                    r2_strs.append(f"{hkl_str}:{r2:.4f}{status}")
                if hkl in area_values:
                    area_strs.append(f"{hkl_str}:{area_values[hkl]:.1f}")
            
            print(f"  R²: {', '.join(r2_strs)}")
            print(f"  Area: {', '.join(area_strs)}")
            
            # For low-quality samples, use all fitted peaks for gray points
            if not all_passed:
                if len(area_values) >= 2:
                    # Use areas with R² >= 0.9 for gray points
                    intensities_low = {k: v for k, v in area_values.items() 
                                       if k in r2_values and r2_values[k] >= 0.9}
                    
                    if len(intensities_low) >= 2:
                        tc_result = analyzer.analyze(intensities_low)
                        tc_values_str = {f"({h[0]}{h[1]}{h[2]})": v for h, v in tc_result.tc_values.items()}
                        
                        all_tc_data.append({
                            'name': filename,
                            'concentration': sample_info['concentration'],
                            'time': sample_info['time'],
                            'tc_values': tc_values_str,
                            'tc_result': tc_result,
                            'high_quality': False,  # Gray points
                        })
                        print(f"  ○ TC (R²<0.995): {', '.join([f'{k}={v:.2f}' for k, v in tc_values_str.items()])}")
                    else:
                        print(f"  ⚠ Insufficient peaks (R² < 0.9)")
                else:
                    print(f"  ⚠ Insufficient peaks for TC analysis")
                continue
            
            if len(intensities) < 3:
                print(f"  ⚠ Insufficient peaks for TC analysis")
                continue
            
            # Calculate TC values
            tc_result = analyzer.analyze(intensities)
            
            # Store for evolution plot
            tc_values_str = {f"({h[0]}{h[1]}{h[2]})": v for h, v in tc_result.tc_values.items()}
            
            all_tc_data.append({
                'name': filename,
                'concentration': sample_info['concentration'],
                'time': sample_info['time'],
                'tc_values': tc_values_str,
                'tc_result': tc_result,
                'high_quality': True,  # Colored points
            })
            
            print(f"  ✓ TC: {', '.join([f'{k}={v:.2f}' for k, v in tc_values_str.items()])}")
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    high_q = sum(1 for d in all_tc_data if d.get('high_quality', True))
    low_q = len(all_tc_data) - high_q
    print(f"\nSuccessfully analyzed {len(all_tc_data)} samples ({high_q} high R², {low_q} low R²)")
    
    if not all_tc_data:
        print("ERROR: No texture data collected!")
        return
    
    # Generate plots
    print("\n" + "-" * 40)
    print("Generating Texture Plots...")
    print("-" * 40)
    
    # 1. TC Evolution by Time (grouped by concentration)
    try:
        fig = plot_tc_evolution(
            all_tc_data,
            x_param='time',
            output_path=str(output_dir / "tc_evolution_by_time.png"),
            show=False,
            dpi=1200,
            normalize=False,  # Standard TC (sum=3)
        )
        plt.close(fig)
        print("  ✓ tc_evolution_by_time.png")
    except Exception as e:
        print(f"  ✗ tc_evolution_by_time.png: {e}")
    
    
    # 1b. Texture Fraction Evolution (normalized, sum=1) - DISABLED per user request
    # User feedback: Not needed when we have fraction_single plot
    # try:
    #     fig = plot_tc_evolution(
    #         all_tc_data,
    #         x_param='time',
    #         output_path=str(output_dir / "texture_fraction_by_time.png"),
    #         show=False,
    #         dpi=1200,
    #         normalize=True,  # Fraction mode (sum=1)
    #     )
    #     plt.close(fig)
    #     print("  ✓ texture_fraction_by_time.png (3 subplots)")
    # except Exception as e:
    #     print(f("  ✗ texture_fraction_by_time.png: {e}")
    
    
    # 1c. Texture Fraction - Single Plot (all peaks on one chart)
    try:
        from axcsas.visualization.texture_plots import plot_texture_fraction_single
        fig = plot_texture_fraction_single(
            all_tc_data,
            x_param='time',
            output_path=str(output_dir / "texture_fraction_single.png"),
            show=False,
            dpi=1200,
        )
        plt.close(fig)
        print("  ✓ texture_fraction_single.png (combined plot)")
    except Exception as e:
        print(f"  ✗ texture_fraction_single.png: {e}")
    
    # 2. TC Evolution by Concentration (grouped by time) - DISABLED per user request
    # User feedback: This plot is confusing and not useful
    # try:
    #     fig = plot_tc_evolution(
    #         all_tc_data,
    #         x_param='concentration',
    #         output_path=str(output_dir / "tc_evolution_by_concentration.png"),
    #         show=False,
    #         dpi=1200,
    #     )
    #     plt.close(fig)
    #     print("  ✓ tc_evolution_by_concentration.png")
    # except Exception as e:
    #     print(f"  ✗ tc_evolution_by_concentration.png: {e}")
    
    
    # 3. Generate polar plots for selected samples (initial and final for each concentration)
    concentrations = sorted(set(d['concentration'] for d in all_tc_data))
    
    for conc in concentrations:
        conc_samples = [d for d in all_tc_data if d['concentration'] == conc]
        if not conc_samples:
            continue
        
        # Sort by time
        conc_samples.sort(key=lambda x: x['time'])
        
        # Plot first and last sample
        for sample in [conc_samples[0], conc_samples[-1]]:
            try:
                # Use integer hours for cleaner filenames
                time_hr = int(round(sample['time']))
                sample_name = f"{sample['concentration']}ml_{time_hr}hr"
                fig = plot_texture_polar(
                    sample['tc_values'],
                    output_path=str(output_dir / f"tc_polar_{sample_name}.png"),
                    show=False,
                    sample_name=sample_name,
                    dpi=1200,
                )
                plt.close(fig)
                print(f"  ✓ tc_polar_{sample_name}.png")
            except Exception as e:
                print(f"  ✗ tc_polar_{sample_name}.png: {e}")
    
    print("\n" + "=" * 60)
    print(f"Complete! Output directory: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
