#!/usr/bin/env python3
"""
TC Texture Coefficient Plots
=============================

Generates texture coefficient (TC) evolution plots for each crystallographic direction.

Output: 3 plots (TC_111.png, TC_200.png, TC_220.png)
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from integration.pipeline import (
    load_bruker_txt,
    find_peak_in_range,
    parse_filename,
)
from physics.texture_enhanced import (
    TextureAnalyzerEnhanced,
    JCPDS_STANDARD_INTENSITY,
)

# Peak positions
PEAK_POSITIONS = {
    (1, 1, 1): 43.3,
    (2, 0, 0): 50.4,
    (2, 2, 0): 74.1,
}


def analyze_all_tc(data_dir: Path):
    """Analyze all XRD files and extract TC data."""
    
    analyzer = TextureAnalyzerEnhanced()
    results = []
    
    for filepath in sorted(data_dir.glob("*.txt")):
        info = parse_filename(str(filepath))
        
        try:
            two_theta, intensity = load_bruker_txt(str(filepath))
        except Exception as e:
            print(f"Error loading {filepath.name}: {e}")
            continue
        
        if len(two_theta) == 0:
            continue
        
        # Find peaks and get intensities
        intensities = {}
        for hkl, expected_pos in PEAK_POSITIONS.items():
            peak = find_peak_in_range(two_theta, intensity, expected_pos, window=2.0)
            if peak:
                intensities[hkl] = peak.intensity
        
        if len(intensities) < 2:
            continue
        
        # Calculate TC
        tc_result = analyzer.analyze(intensities)
        
        sample_data = {
            'name': info['name'],
            'concentration': info['concentration_ml'],
            'time': info['time_hours'],
            'tc_values': tc_result.tc_values,
            'dominant_hkl': tc_result.dominant_hkl,
            'is_random': tc_result.is_random,
        }
        
        results.append(sample_data)
        print(f"Analyzed: {filepath.name}")
    
    return results


def plot_tc_by_direction(results, output_dir: Path):
    """Generate TC plots for each crystallographic direction."""
    
    directions = [(1, 1, 1), (2, 0, 0), (2, 2, 0)]
    concentrations = [0, 4.5, 9, 18]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    markers = ['o', 's', '^', 'D']
    
    for hkl in directions:
        hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        for conc, color, marker in zip(concentrations, colors, markers):
            times = []
            tc_values = []
            
            for r in results:
                if r['concentration'] == conc and hkl in r['tc_values']:
                    times.append(r['time'])
                    tc_values.append(r['tc_values'][hkl])
            
            if times:
                sorted_data = sorted(zip(times, tc_values))
                times, tc_values = zip(*sorted_data)
                
                ax.plot(times, tc_values, marker=marker, color=color,
                       label=f'{conc} ml', linewidth=2, markersize=8)
        
        # Reference line at TC = 1 (random)
        ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=2, 
                   label='Random (TC=1)')
        
        # Random zone shading (0.9-1.1)
        ax.axhspan(0.9, 1.1, color='gray', alpha=0.2, label='Random Zone')
        
        ax.set_xlabel('Plating Time (hours)', fontsize=12)
        ax.set_ylabel(f'TC{hkl_str}', fontsize=12)
        ax.set_title(f'Texture Coefficient {hkl_str} Evolution', 
                     fontsize=14, fontweight='bold')
        ax.legend(title='Leveler', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 2.5)
        
        # Add physical interpretation note
        ax.text(0.98, 0.02, 
                'TC > 1: Preferred\nTC = 1: Random\nTC < 1: Suppressed',
                transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        filename = f'TC_{hkl[0]}{hkl[1]}{hkl[2]}.png'
        plt.savefig(output_dir / filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")


def plot_defect_data(results, output_dir: Path):
    """Generate peak separation and lattice constant plots."""
    
    # Need to re-analyze for defect data
    from physics.defect_analysis import analyze_stacking_faults, analyze_lattice
    
    concentrations = [0, 4.5, 9, 18]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    markers = ['o', 's', '^', 'D']
    
    # Load and analyze defect data
    data_dir = Path(__file__).parent / 'data' / 'raw' / '202511'
    defect_results = []
    
    for filepath in sorted(data_dir.glob("*.txt")):
        info = parse_filename(str(filepath))
        
        try:
            two_theta, intensity = load_bruker_txt(str(filepath))
        except:
            continue
        
        # Find (111) and (200) peaks
        peak_111 = find_peak_in_range(two_theta, intensity, 43.3)
        peak_200 = find_peak_in_range(two_theta, intensity, 50.4)
        peak_220 = find_peak_in_range(two_theta, intensity, 74.1)
        
        if peak_111 and peak_200:
            sf_result = analyze_stacking_faults(peak_111.two_theta, peak_200.two_theta)
            lattice = analyze_lattice(peak_220.two_theta if peak_220 else peak_200.two_theta, 
                                      (2, 2, 0) if peak_220 else (2, 0, 0))
            
            defect_results.append({
                'concentration': info['concentration_ml'],
                'time': info['time_hours'],
                'peak_separation': sf_result.peak_separation_deg,
                'lattice_constant': lattice.lattice_constant,
            })
    
    # Plot 1: Peak Separation
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for conc, color, marker in zip(concentrations, colors, markers):
        times = []
        sep_values = []
        
        for r in defect_results:
            if r['concentration'] == conc:
                times.append(r['time'])
                sep_values.append(r['peak_separation'])
        
        if times:
            sorted_data = sorted(zip(times, sep_values))
            times, sep_values = zip(*sorted_data)
            ax.plot(times, sep_values, marker=marker, color=color,
                   label=f'{conc} ml', linewidth=2, markersize=8)
    
    # Standard peak separation line
    ax.axhline(y=7.136, color='red', linestyle='--', linewidth=2, 
               label='Standard (7.136°)')
    
    ax.set_xlabel('Plating Time (hours)', fontsize=12)
    ax.set_ylabel('Peak Separation (111)-(200) (°)', fontsize=12)
    ax.set_title('Peak Separation Evolution', fontsize=14, fontweight='bold')
    ax.legend(title='Leveler', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'peak_separation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: peak_separation.png")
    
    # Plot 2: Lattice Constant
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for conc, color, marker in zip(concentrations, colors, markers):
        times = []
        a_values = []
        
        for r in defect_results:
            if r['concentration'] == conc:
                times.append(r['time'])
                a_values.append(r['lattice_constant'])
        
        if times:
            sorted_data = sorted(zip(times, a_values))
            times, a_values = zip(*sorted_data)
            ax.plot(times, a_values, marker=marker, color=color,
                   label=f'{conc} ml', linewidth=2, markersize=8)
    
    # Standard lattice constant
    ax.axhline(y=3.6150, color='green', linestyle='--', linewidth=2, 
               label='Standard (3.6150 Å)')
    
    # Warning threshold
    ax.axhspan(3.618, 3.63, color='red', alpha=0.2, label='Impurity Warning')
    
    ax.set_xlabel('Plating Time (hours)', fontsize=12)
    ax.set_ylabel('Lattice Constant a (Å)', fontsize=12)
    ax.set_title('Lattice Constant Evolution', fontsize=14, fontweight='bold')
    ax.legend(title='Leveler', fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'lattice_constant.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: lattice_constant.png")


def main():
    """Main function."""
    
    data_dir = Path(__file__).parent / 'data' / 'raw' / '202511'
    output_dir = Path(__file__).parent / 'outputs' / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("AXCSAS TC & Defect Analysis")
    print("=" * 60)
    
    # Analyze TC
    print("\nAnalyzing TC data...")
    results = analyze_all_tc(data_dir)
    print(f"\nTotal samples: {len(results)}")
    
    # Generate TC plots
    print("\nGenerating TC plots...")
    plot_tc_by_direction(results, output_dir)
    
    # Generate defect plots
    print("\nGenerating defect plots...")
    plot_defect_data(results, output_dir)
    
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print(f"Plots saved to: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
