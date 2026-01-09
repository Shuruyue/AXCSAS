#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FWHM Visualization for Cu XRD Peaks

Generates plots showing FWHM evolution over time for each Cu peak (111, 200, 220).
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

# Style settings for publication-quality figures
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'axes.grid': True,
    'grid.alpha': 0.3
})

# Color scheme for concentrations
COLORS = {
    '0ml': '#1f77b4',
    '4.5ml': '#ff7f0e', 
    '9ml': '#2ca02c',
    '18ml': '#d62728'
}

MARKERS = {
    '0ml': 'o',
    '4.5ml': 's',
    '9ml': '^',
    '18ml': 'D'
}

# Instrument detection limit (from 計畫.txt)
# When FWHM < 0.05°, crystallite size > 200nm exceeds detection limit
INSTRUMENT_LIMIT = 0.05  # degrees


def load_fwhm_data(output_dir: Path) -> dict:
    """Load all FWHM CSV files."""
    data = {}
    for conc in ['0ml', '4.5ml', '9ml', '18ml']:
        filepath = output_dir / f"fwhm_{conc}.csv"
        if filepath.exists():
            df = pd.read_csv(filepath)
            # Ensure Peak column is string type
            df['Peak'] = df['Peak'].astype(str)
            df = df.set_index('Peak')
            data[conc] = df
    return data


def plot_fwhm_by_peak(data: dict, output_dir: Path):
    """
    Create 3 subplots, one for each peak direction.
    Each subplot shows all 4 concentrations.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    peaks = ['111', '200', '220']
    
    # Calculate global Y max for consistent axis limits
    y_max = 0
    for conc, df in data.items():
        for peak in peaks:
            if peak in df.index:
                y_max = max(y_max, df.loc[peak].values.max())
    y_max = y_max * 1.1  # Add 10% margin
    
    for ax, peak in zip(axes, peaks):
        for conc, df in data.items():
            if peak in df.index:
                # Get time points (column names) and values
                time_labels = df.columns.tolist()
                time_values = [int(t.replace('h', '')) for t in time_labels]
                fwhm_values = df.loc[peak].values
                
                ax.plot(time_values, fwhm_values, 
                       marker=MARKERS[conc], 
                       color=COLORS[conc],
                       label=conc,
                       linewidth=2,
                       markersize=6)
        
        ax.set_xlabel('Plating Time (h)')
        ax.set_ylabel('FWHM (° 2θ)')
        ax.set_title(f'Cu ({peak})')
        ax.legend(title='Leveler')
        ax.set_xlim(-1, 25)
        ax.set_ylim(0, y_max)  # Consistent Y-axis for all subplots
        # Add undetectable region (FWHM < 0.05°)
        ax.axhspan(0, INSTRUMENT_LIMIT, facecolor='gray', alpha=0.25, hatch='///', edgecolor='gray', linewidth=0)
    
    plt.tight_layout()
    output_path = output_dir / "fwhm_by_peak.png"
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    return output_path


def plot_fwhm_by_concentration(data: dict, output_dir: Path):
    """
    Create 4 subplots, one for each concentration.
    Each subplot shows all 3 peaks.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    peak_colors = {'111': '#e41a1c', '200': '#377eb8', '220': '#4daf4a'}
    peak_markers = {'111': 'o', '200': 's', '220': '^'}
    peaks = ['111', '200', '220']
    
    # Calculate global Y max for consistent axis limits
    y_max = 0
    for conc, df in data.items():
        for peak in peaks:
            if peak in df.index:
                y_max = max(y_max, df.loc[peak].values.max())
    y_max = y_max * 1.1  # Add 10% margin
    
    for ax, (conc, df) in zip(axes, data.items()):
        time_labels = df.columns.tolist()
        time_values = [int(t.replace('h', '')) for t in time_labels]
        
        for peak in peaks:
            if peak in df.index:
                fwhm_values = df.loc[peak].values
                ax.plot(time_values, fwhm_values,
                       marker=peak_markers[peak],
                       color=peak_colors[peak],
                       label=f'({peak})',
                       linewidth=2,
                       markersize=6)
        
        ax.set_xlabel('Plating Time (h)')
        ax.set_ylabel('FWHM (° 2θ)')
        ax.set_title(f'Leveler: {conc}')
        ax.legend(title='Cu Peak')
        ax.set_xlim(-1, 25)
        ax.set_ylim(0, y_max)  # Consistent Y-axis for all subplots
        # Add undetectable region (FWHM < 0.05°)
        ax.axhspan(0, INSTRUMENT_LIMIT, facecolor='gray', alpha=0.25, hatch='///', edgecolor='gray', linewidth=0)
    
    plt.tight_layout()
    output_path = output_dir / "fwhm_by_concentration.png"
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    return output_path


def plot_fwhm_combined(data: dict, output_dir: Path):
    """
    Create a single combined figure with all data.
    3 rows (peaks) × 1 column, with all concentrations overlaid.
    """
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    peaks = ['111', '200', '220']
    
    # Calculate global Y max for consistent axis limits
    y_max = 0
    for conc, df in data.items():
        for peak in peaks:
            if peak in df.index:
                y_max = max(y_max, df.loc[peak].values.max())
    y_max = y_max * 1.1  # Add 10% margin
    
    for ax, peak in zip(axes, peaks):
        for conc, df in data.items():
            if peak in df.index:
                time_labels = df.columns.tolist()
                time_values = [int(t.replace('h', '')) for t in time_labels]
                fwhm_values = df.loc[peak].values
                
                ax.plot(time_values, fwhm_values,
                       marker=MARKERS[conc],
                       color=COLORS[conc],
                       label=conc,
                       linewidth=2.5,
                       markersize=8)
        
        ax.set_xlabel('Plating Time (h)')
        ax.set_ylabel('FWHM (° 2θ)')
        ax.set_title(f'FWHM Evolution - Cu ({peak})')
        ax.legend(title='Leveler Conc.', loc='upper right')
        ax.set_xlim(-1, 25)
        ax.set_ylim(0, y_max)  # Consistent Y-axis for all subplots
    
    plt.tight_layout()
    output_path = output_dir / "fwhm_evolution.png"
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    return output_path


def main():
    script_path = Path(__file__).resolve()
    project_dir = script_path.parent.parent
    output_dir = project_dir / "outputs"
    
    print("=" * 60)
    print("FWHM Visualization")
    print("=" * 60)
    
    # Load data
    data = load_fwhm_data(output_dir)
    print(f"Loaded {len(data)} concentration datasets")
    
    # Generate plots
    print("\nGenerating plots...")
    plot_fwhm_by_peak(data, output_dir)
    plot_fwhm_by_concentration(data, output_dir)
    
    print("\n" + "=" * 60)
    print("Visualization complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
