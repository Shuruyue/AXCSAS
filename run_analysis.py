#!/usr/bin/env python3
"""
XRD Analysis with Direct Data Parsing
Analyzes Bruker TXT format XRD data and generates plots.
"""

import sys
from pathlib import Path
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import re

# Constants
CU_KA1 = 1.54056  # Angstrom
SCHERRER_K = 0.89

# Cu peak positions and standard intensities
CU_PEAKS = {
    (1, 1, 1): {'two_theta': 43.3, 'std_intensity': 100},
    (2, 0, 0): {'two_theta': 50.4, 'std_intensity': 46},
    (2, 2, 0): {'two_theta': 74.1, 'std_intensity': 20},
    (3, 1, 1): {'two_theta': 89.9, 'std_intensity': 17},
}


def load_bruker_txt(filepath: str):
    """Load Bruker TXT format XRD data."""
    two_theta = []
    intensity = []
    
    in_data_section = False
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('[Data]'):
                in_data_section = True
                continue
            
            if in_data_section:
                # Skip header line
                if 'Angle' in line or 'PSD' in line:
                    continue
                
                # Parse data line
                parts = line.replace(',', '').split()
                if len(parts) >= 2:
                    try:
                        angle = float(parts[0])
                        inten = float(parts[1])
                        two_theta.append(angle)
                        intensity.append(inten)
                    except ValueError:
                        continue
    
    return np.array(two_theta), np.array(intensity)


def find_peak(two_theta, intensity, center, window=2.0):
    """Find peak near center position."""
    mask = (two_theta > center - window) & (two_theta < center + window)
    if not np.any(mask):
        return None
    
    idx_window = np.where(mask)[0]
    local_intensity = intensity[idx_window]
    peak_idx = idx_window[np.argmax(local_intensity)]
    
    peak_pos = two_theta[peak_idx]
    peak_height = intensity[peak_idx]
    
    # Estimate FWHM by finding half-max points
    half_max = (peak_height + intensity[idx_window[0]]) / 2
    
    left_idx = peak_idx
    while left_idx > 0 and intensity[left_idx] > half_max:
        left_idx -= 1
    
    right_idx = peak_idx
    while right_idx < len(intensity) - 1 and intensity[right_idx] > half_max:
        right_idx += 1
    
    fwhm = two_theta[right_idx] - two_theta[left_idx]
    
    return {
        'position': peak_pos,
        'height': peak_height,
        'fwhm': max(fwhm, 0.1),  # Minimum FWHM
        'idx': peak_idx
    }


def scherrer_size(two_theta_deg, fwhm_deg, wavelength=CU_KA1, k=SCHERRER_K):
    """Calculate crystallite size using Scherrer equation."""
    theta_rad = np.radians(two_theta_deg / 2)
    fwhm_rad = np.radians(fwhm_deg)
    
    size_angstrom = (k * wavelength) / (fwhm_rad * np.cos(theta_rad))
    size_nm = size_angstrom / 10
    
    return size_nm


def texture_coefficient(intensities: dict):
    """Calculate Harris texture coefficients."""
    tc = {}
    
    # Calculate sum of I/I0
    sum_ratio = 0
    for hkl, obs_int in intensities.items():
        std_int = CU_PEAKS[hkl]['std_intensity']
        sum_ratio += obs_int / std_int
    
    n = len(intensities)
    
    for hkl, obs_int in intensities.items():
        std_int = CU_PEAKS[hkl]['std_intensity']
        tc[hkl] = (obs_int / std_int) / (sum_ratio / n)
    
    return tc


def analyze_sample(filepath: str, output_dir: Path):
    """Analyze a single XRD file."""
    sample_name = Path(filepath).stem
    
    print(f"\n{'='*50}")
    print(f"Sample: {sample_name}")
    print('='*50)
    
    # Load data
    two_theta, intensity = load_bruker_txt(filepath)
    
    if len(two_theta) == 0:
        print("  Error: No data loaded")
        return None
    
    print(f"  Data range: {two_theta.min():.1f}° - {two_theta.max():.1f}°")
    print(f"  Data points: {len(two_theta)}")
    
    # Smooth data
    intensity_smooth = savgol_filter(intensity, 11, 3)
    
    # Find peaks
    peaks_found = {}
    intensities_for_tc = {}
    
    for hkl, info in CU_PEAKS.items():
        if two_theta.min() <= info['two_theta'] <= two_theta.max():
            peak = find_peak(two_theta, intensity_smooth, info['two_theta'], window=2.0)
            if peak and peak['height'] > np.median(intensity_smooth) * 1.3:
                peaks_found[hkl] = peak
                intensities_for_tc[hkl] = peak['height']
                print(f"  ({hkl[0]}{hkl[1]}{hkl[2]}): 2θ={peak['position']:.2f}°, FWHM={peak['fwhm']:.3f}°")
    
    if not peaks_found:
        print("  No peaks found")
        return None
    
    # Calculate crystallite sizes
    sizes = {}
    print(f"\n  Crystallite Size (Scherrer):")
    for hkl, peak in peaks_found.items():
        size = scherrer_size(peak['position'], peak['fwhm'])
        sizes[hkl] = size
        status = "✓" if 5 < size < 200 else "⚠"
        print(f"    ({hkl[0]}{hkl[1]}{hkl[2]}): D = {size:.1f} nm {status}")
    
    avg_size = np.mean(list(sizes.values()))
    print(f"    Average: {avg_size:.1f} nm")
    
    # Texture analysis
    if len(intensities_for_tc) >= 2:
        tc = texture_coefficient(intensities_for_tc)
        print(f"\n  Texture Coefficients:")
        for hkl, val in tc.items():
            marker = "★" if val > 1.5 else ""
            print(f"    TC({hkl[0]}{hkl[1]}{hkl[2]}) = {val:.2f} {marker}")
    else:
        tc = {}
    
    # Create plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Full pattern
    ax = axes[0, 0]
    ax.plot(two_theta, intensity, 'b-', alpha=0.3, linewidth=0.5, label='Raw')
    ax.plot(two_theta, intensity_smooth, 'r-', linewidth=1, label='Smoothed')
    
    for hkl, peak in peaks_found.items():
        ax.axvline(peak['position'], color='green', linestyle='--', alpha=0.5)
        ax.annotate(f"({hkl[0]}{hkl[1]}{hkl[2]})", 
                   xy=(peak['position'], peak['height']),
                   xytext=(peak['position']+1, peak['height']*1.05),
                   fontsize=9, fontweight='bold')
    
    ax.set_xlabel('2θ (degrees)')
    ax.set_ylabel('Intensity (counts)')
    ax.set_title(f'XRD Pattern: {sample_name}', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # (111) peak detail
    ax = axes[0, 1]
    if (1, 1, 1) in peaks_found:
        peak = peaks_found[(1, 1, 1)]
        mask = (two_theta > peak['position'] - 4) & (two_theta < peak['position'] + 4)
        ax.plot(two_theta[mask], intensity_smooth[mask], 'b-', linewidth=2)
        ax.axvline(peak['position'], color='red', linestyle='--', 
                  label=f"2θ = {peak['position']:.2f}°")
        ax.fill_between(two_theta[mask], intensity_smooth[mask], alpha=0.3)
        
        # Mark FWHM
        fwhm_y = peak['height'] / 2
        ax.axhline(fwhm_y, color='green', linestyle=':', alpha=0.7,
                  label=f"FWHM = {peak['fwhm']:.3f}°")
        
        ax.set_xlabel('2θ (degrees)')
        ax.set_ylabel('Intensity')
        ax.set_title(f'Cu (111) Peak - D = {sizes[(1,1,1)]:.1f} nm', fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
    
    # Size bar chart
    ax = axes[1, 0]
    if sizes:
        labels = [f"({h}{k}{l})" for h, k, l in sizes.keys()]
        values = list(sizes.values())
        colors = ['green' if 5 < v < 200 else 'orange' for v in values]
        
        bars = ax.bar(labels, values, color=colors, edgecolor='black')
        ax.axhline(avg_size, color='red', linestyle='--', linewidth=2,
                  label=f'Average: {avg_size:.1f} nm')
        ax.set_xlabel('(hkl) Reflection')
        ax.set_ylabel('Crystallite Size (nm)')
        ax.set_title('Scherrer Crystallite Size', fontweight='bold')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
    
    # Texture coefficients
    ax = axes[1, 1]
    if tc:
        labels = [f"({h}{k}{l})" for h, k, l in tc.keys()]
        values = list(tc.values())
        colors = ['coral' if v > 1.5 else 'steelblue' for v in values]
        
        ax.bar(labels, values, color=colors, edgecolor='black')
        ax.axhline(1.0, color='red', linestyle='--', linewidth=2,
                  label='Random texture (TC=1)')
        ax.set_xlabel('(hkl) Reflection')
        ax.set_ylabel('Texture Coefficient')
        ax.set_title('Harris Texture Analysis', fontweight='bold')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    plot_file = output_dir / f"{sample_name}_analysis.png"
    fig.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  ✓ Saved: {plot_file.name}")
    
    return {
        'sample': sample_name,
        'avg_size': avg_size,
        'sizes': sizes,
        'tc': tc,
        'plot': str(plot_file)
    }


def main():
    output_dir = Path("outputs/plots")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find XRD files
    data_dir = Path("XRD data/202511")
    files = list(data_dir.glob("*.txt"))
    
    if not files:
        data_dir = Path("XRD data/old")
        files = list(data_dir.glob("*.txt"))
    
    print(f"\n{'='*60}")
    print("AXCSAS XRD Analysis")
    print('='*60)
    print(f"Found {len(files)} XRD files")
    
    # Process a subset for demo
    sample_files = files[:8]
    
    results = []
    for f in sample_files:
        try:
            result = analyze_sample(str(f), output_dir)
            if result:
                results.append(result)
        except Exception as e:
            print(f"  Error: {e}")
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print('='*60)
    
    if results:
        for r in results:
            tc_111 = r['tc'].get((1,1,1), 0)
            print(f"{r['sample']}: D = {r['avg_size']:.1f} nm, TC(111) = {tc_111:.2f}")
        
        # Summary plot
        fig, ax = plt.subplots(figsize=(12, 6))
        
        names = [r['sample'].replace('20251125_', '').replace('20251126_', '') for r in results]
        sizes = [r['avg_size'] for r in results]
        
        ax.bar(names, sizes, color='steelblue', edgecolor='black')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Crystallite Size (nm)')
        ax.set_title('Crystallite Size Comparison', fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        summary_plot = output_dir / "size_comparison.png"
        fig.savefig(summary_plot, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"\n✓ Summary plot: {summary_plot}")
    
    print(f"\nPlots saved to: {output_dir}")


if __name__ == "__main__":
    main()
