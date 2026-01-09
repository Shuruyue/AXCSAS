#!/usr/bin/env python3
"""
FWHM Analysis and Plotting Script
==================================

Analyzes all XRD data files and generates FWHM plots.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from integration.pipeline import (
    AXCSASPipeline,
    load_bruker_txt,
    find_peak_in_range,
    parse_filename,
)

# Peak positions
PEAK_POSITIONS = {
    (1, 1, 1): 43.3,
    (2, 0, 0): 50.4,
    (2, 2, 0): 74.1,
    (3, 1, 1): 89.9,
}

def analyze_all_files(data_dir: Path):
    """Analyze all XRD files and extract FWHM data."""
    
    results = []
    
    for filepath in sorted(data_dir.glob("*.txt")):
        # Parse filename
        info = parse_filename(str(filepath))
        
        # Load data
        try:
            two_theta, intensity = load_bruker_txt(str(filepath))
        except Exception as e:
            print(f"Error loading {filepath.name}: {e}")
            continue
        
        if len(two_theta) == 0:
            continue
        
        # Find peaks and FWHM
        sample_data = {
            'name': info['name'],
            'concentration': info['concentration_ml'],
            'time': info['time_hours'],
            'peaks': {}
        }
        
        for hkl, expected_pos in PEAK_POSITIONS.items():
            peak = find_peak_in_range(two_theta, intensity, expected_pos, window=2.0)
            if peak:
                sample_data['peaks'][hkl] = {
                    'two_theta': peak.two_theta,
                    'intensity': peak.intensity,
                    'fwhm': peak.fwhm,
                }
        
        results.append(sample_data)
        print(f"Analyzed: {filepath.name}")
    
    return results


def plot_fwhm_by_peak(results, output_path):
    """Create FWHM plot organized by peak."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('FWHM vs. Plating Time by Peak', fontsize=14, fontweight='bold')
    
    concentrations = [0, 4.5, 9, 18]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    markers = ['o', 's', '^', 'D']
    
    peaks = [(1,1,1), (2,0,0), (2,2,0), (3,1,1)]
    
    # Boundary conditions
    # Caglioti W parameter (from config, placeholder value)
    # FWHM_inst = sqrt(W) at low angles
    CAGLIOTI_W = 0.01  # deg² (typical value, should be calibrated)
    INSTRUMENT_LIMIT = np.sqrt(CAGLIOTI_W)  # ~0.1° from Caglioti
    
    Y_MAX = 1.6  # Based on actual data max ~1.43°
    Y_MIN = 0.0
    
    for idx, hkl in enumerate(peaks):
        ax = axes[idx // 2, idx % 2]
        hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
        
        # Add gray shaded region for model limitation (below instrument limit)
        ax.axhspan(Y_MIN, INSTRUMENT_LIMIT, color='gray', alpha=0.3, 
                   label='Model Limit (FWHM < Inst.)')
        
        for conc, color, marker in zip(concentrations, colors, markers):
            times = []
            fwhms = []
            
            for r in results:
                if r['concentration'] == conc and hkl in r['peaks']:
                    times.append(r['time'])
                    fwhms.append(r['peaks'][hkl]['fwhm'])
            
            if times:
                # Sort by time
                sorted_data = sorted(zip(times, fwhms))
                times, fwhms = zip(*sorted_data)
                
                ax.plot(times, fwhms, marker=marker, color=color, 
                       label=f'{conc} ml', linewidth=1.5, markersize=6)
        
        # Instrument limit line
        ax.axhline(y=INSTRUMENT_LIMIT, color='red', linestyle='--', 
                   linewidth=2, label=f'Inst. Limit ({INSTRUMENT_LIMIT}°)')
        
        ax.set_xlabel('Plating Time (hours)')
        ax.set_ylabel('FWHM (°)')
        ax.set_title(f'{hkl_str} Peak')
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(Y_MIN, Y_MAX)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_fwhm_by_concentration(results, output_path):
    """Create FWHM plot organized by concentration."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('FWHM vs. Plating Time by Leveler Concentration', 
                 fontsize=14, fontweight='bold')
    
    concentrations = [0, 4.5, 9, 18]
    peaks = [(1,1,1), (2,0,0), (2,2,0)]
    colors = ['#e41a1c', '#377eb8', '#4daf4a']
    markers = ['o', 's', '^']
    
    # Boundary conditions (same as plot_fwhm_by_peak)
    CAGLIOTI_W = 0.01  # deg²
    INSTRUMENT_LIMIT = np.sqrt(CAGLIOTI_W)  # ~0.1° from Caglioti
    
    Y_MAX = 1.6  # Based on actual data max
    Y_MIN = 0.0
    
    for idx, conc in enumerate(concentrations):
        ax = axes[idx // 2, idx % 2]
        
        # Add gray shaded region for model limitation
        ax.axhspan(Y_MIN, INSTRUMENT_LIMIT, color='gray', alpha=0.3,
                   label='Model Limit')
        
        for hkl, color, marker in zip(peaks, colors, markers):
            times = []
            fwhms = []
            
            for r in results:
                if r['concentration'] == conc and hkl in r['peaks']:
                    times.append(r['time'])
                    fwhms.append(r['peaks'][hkl]['fwhm'])
            
            if times:
                sorted_data = sorted(zip(times, fwhms))
                times, fwhms = zip(*sorted_data)
                
                hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
                ax.plot(times, fwhms, marker=marker, color=color,
                       label=hkl_str, linewidth=1.5, markersize=6)
        
        # Instrument limit
        ax.axhline(y=INSTRUMENT_LIMIT, color='red', linestyle='--',
                   linewidth=2, label=f'Inst. ({INSTRUMENT_LIMIT}°)')
        
        ax.set_xlabel('Plating Time (hours)')
        ax.set_ylabel('FWHM (°)')
        ax.set_title(f'Leveler: {conc} ml')
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(Y_MIN, Y_MAX)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def plot_crystallite_size(results, output_path):
    """Plot estimated crystallite size from FWHM."""
    
    # Scherrer constants
    WAVELENGTH = 1.54056  # Å
    K = 0.9
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    concentrations = [0, 4.5, 9, 18]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    markers = ['o', 's', '^', 'D']
    
    for conc, color, marker in zip(concentrations, colors, markers):
        times = []
        sizes = []
        
        for r in results:
            if r['concentration'] == conc and (1,1,1) in r['peaks']:
                fwhm = r['peaks'][(1,1,1)]['fwhm']
                two_theta = r['peaks'][(1,1,1)]['two_theta']
                
                # Scherrer calculation
                fwhm_rad = np.radians(fwhm)
                theta_rad = np.radians(two_theta / 2)
                size_A = K * WAVELENGTH / (fwhm_rad * np.cos(theta_rad))
                size_nm = size_A / 10
                
                times.append(r['time'])
                sizes.append(size_nm)
        
        if times:
            sorted_data = sorted(zip(times, sizes))
            times, sizes = zip(*sorted_data)
            
            ax.plot(times, sizes, marker=marker, color=color,
                   label=f'{conc} ml', linewidth=2, markersize=8)
    
    ax.set_xlabel('Plating Time (hours)', fontsize=12)
    ax.set_ylabel('Crystallite Size (nm)', fontsize=12)
    ax.set_title('Crystallite Size Evolution from (111) Peak', fontsize=14, fontweight='bold')
    ax.legend(title='Leveler', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 50)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def main():
    """Main function."""
    
    # Data directory
    data_dir = Path(__file__).parent / 'data' / 'raw' / '202511'
    output_dir = Path(__file__).parent / 'outputs' / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("AXCSAS FWHM Analysis")
    print("=" * 60)
    
    # Analyze all files
    print("\nAnalyzing XRD files...")
    results = analyze_all_files(data_dir)
    print(f"\nTotal samples analyzed: {len(results)}")
    
    # Generate plots
    print("\nGenerating plots...")
    
    plot_fwhm_by_peak(results, output_dir / 'fwhm_by_peak.png')
    plot_fwhm_by_concentration(results, output_dir / 'fwhm_by_concentration.png')
    
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print(f"Plots saved to: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
