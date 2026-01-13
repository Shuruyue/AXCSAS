#!/usr/bin/env python3
"""
AXCSAS Plot Generator
=====================

Generate FWHM evolution and Scherrer Size plots from XRD data.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import re

# Configure matplotlib for better output
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150
plt.rcParams['figure.figsize'] = (12, 8)

# Import AXCSAS modules
from axcsas.analysis.pipeline import AXCSASPipeline, AnalysisConfig, PipelineResult


@dataclass
class SampleData:
    """Parsed sample data."""
    filepath: str
    filename: str
    concentration_ml: float
    time_hours: float
    result: Optional[PipelineResult] = None


def parse_sample_info(filename: str) -> Tuple[float, float]:
    """
    Parse leveler concentration and plating time from filename.
    
    Format: YYYYMMDD_Xml_Xh.txt or YYYYMMDD_Xml_Xh_Xmin.txt
    """
    # Extract concentration (e.g., "4.5ml" -> 4.5)
    conc_match = re.search(r'_(\d+\.?\d*)ml_', filename)
    concentration = float(conc_match.group(1)) if conc_match else 0.0
    
    # Extract time in hours
    hours_match = re.search(r'_(\d+)h', filename)
    hours = float(hours_match.group(1)) if hours_match else 0.0
    
    # Check for minutes
    min_match = re.search(r'_(\d+)min', filename)
    if min_match:
        hours += float(min_match.group(1)) / 60.0
    
    return concentration, hours


def analyze_all_samples(data_dir: Path) -> List[SampleData]:
    """Analyze all XRD samples in directory."""
    samples = []
    pipeline = AXCSASPipeline()
    
    txt_files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(txt_files)} XRD data files")
    
    for filepath in txt_files:
        try:
            concentration, time_hours = parse_sample_info(filepath.name)
            
            # Run AXCSAS analysis
            result = pipeline.analyze(str(filepath))
            
            sample = SampleData(
                filepath=str(filepath),
                filename=filepath.name,
                concentration_ml=concentration,
                time_hours=time_hours,
                result=result
            )
            samples.append(sample)
            print(f"  ✓ {filepath.name}")
            
        except Exception as e:
            print(f"  ✗ {filepath.name}: {e}")
    
    return samples


def plot_fwhm_evolution(samples: List[SampleData], output_dir: Path):
    """
    Generate FWHM evolution plots.
    
    Creates:
    1. FWHM vs Time for each concentration (by peak)
    2. FWHM heatmap
    """
    # Group by concentration
    concentrations = sorted(set(s.concentration_ml for s in samples))
    peaks = [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1)]
    peak_labels = ['(111)', '(200)', '(220)', '(311)']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    # =========================================================================
    # Plot 1: FWHM Evolution by Concentration
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    # Collect all FWHM values for consistent Y-axis
    all_fwhms = []
    for s in samples:
        if s.result and s.result.peaks:
            for p in s.result.peaks:
                all_fwhms.append(p.fwhm)
    y_max_fwhm = max(all_fwhms) * 1.1 if all_fwhms else 1.5
    
    for idx, (peak, label) in enumerate(zip(peaks, peak_labels)):
        ax = axes[idx]
        
        for conc in concentrations:
            conc_samples = [s for s in samples if s.concentration_ml == conc and s.result]
            conc_samples.sort(key=lambda x: x.time_hours)
            
            times = []
            fwhms = []
            
            for s in conc_samples:
                if s.result and s.result.peaks:
                    for p in s.result.peaks:
                        if p.hkl == peak:
                            times.append(s.time_hours)
                            fwhms.append(p.fwhm)
                            break
            
            if times and fwhms:
                ax.plot(times, fwhms, 'o-', label=f'{conc} ml', markersize=6)
        
        ax.set_xlabel('Plating Time (hours)')
        ax.set_ylabel('FWHM (°)')
        ax.set_title(f'{label} Peak FWHM Evolution')
        ax.legend(title='Leveler Conc.')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0.05, color='red', linestyle='--', alpha=0.5, label='Instrument Limit')
        ax.set_ylim(0, y_max_fwhm)  # Consistent Y-axis
    
    plt.tight_layout()
    output_path = output_dir / 'fwhm_evolution_by_peak.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    
    # =========================================================================
    # Plot 2: FWHM by Concentration (all peaks together)
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for idx, conc in enumerate(concentrations):
        if idx >= 4:
            break
        ax = axes[idx]
        
        conc_samples = [s for s in samples if s.concentration_ml == conc and s.result]
        conc_samples.sort(key=lambda x: x.time_hours)
        
        for peak, label, color in zip(peaks, peak_labels, colors):
            times = []
            fwhms = []
            
            for s in conc_samples:
                if s.result and s.result.peaks:
                    for p in s.result.peaks:
                        if p.hkl == peak:
                            times.append(s.time_hours)
                            fwhms.append(p.fwhm)
                            break
            
            if times and fwhms:
                ax.plot(times, fwhms, 'o-', label=label, color=color, markersize=6)
        
        ax.set_xlabel('Plating Time (hours)')
        ax.set_ylabel('FWHM (°)')
        ax.set_title(f'FWHM Evolution - {conc} ml Leveler')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0.05, color='red', linestyle='--', alpha=0.5)
        ax.set_ylim(0, y_max_fwhm)  # Consistent Y-axis
    
    plt.tight_layout()
    output_path = output_dir / 'fwhm_by_concentration.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    
    return 2  # Number of plots generated


def plot_scherrer_size(samples: List[SampleData], output_dir: Path):
    """
    Generate Scherrer crystallite size plots.
    
    Creates:
    1. Size vs Time for each concentration (by peak)
    2. Size comparison across concentrations
    """
    # Group by concentration
    concentrations = sorted(set(s.concentration_ml for s in samples))
    peaks = [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1)]
    peak_labels = ['(111)', '(200)', '(220)', '(311)']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    # Helper function to get size for a specific peak
    def get_size_for_peak(result, target_hkl):
        if result and result.scherrer_results:
            for sr in result.scherrer_results:
                if sr.hkl == target_hkl:
                    return sr.size_nm
        return None
    
    # =========================================================================
    # Plot 1: Scherrer Size Evolution by Peak
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    # Collect all sizes for consistent Y-axis
    all_sizes = []
    for s in samples:
        if s.result and s.result.scherrer_results:
            for sr in s.result.scherrer_results:
                if sr.size_nm and sr.size_nm > 0:
                    all_sizes.append(sr.size_nm)
    y_max_size = max(all_sizes) * 1.1 if all_sizes else 500
    
    for idx, (peak, label) in enumerate(zip(peaks, peak_labels)):
        ax = axes[idx]
        
        for conc in concentrations:
            conc_samples = [s for s in samples if s.concentration_ml == conc and s.result]
            conc_samples.sort(key=lambda x: x.time_hours)
            
            times = []
            sizes = []
            
            for s in conc_samples:
                size = get_size_for_peak(s.result, peak)
                if size and size > 0:
                    times.append(s.time_hours)
                    sizes.append(size)
            
            if times and sizes:
                ax.plot(times, sizes, 'o-', label=f'{conc} ml', markersize=6)
        
        ax.set_xlabel('Plating Time (hours)')
        ax.set_ylabel('Crystallite Size (nm)')
        ax.set_title(f'{label} Scherrer Crystallite Size')
        ax.legend(title='Leveler Conc.')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, y_max_size)  # Consistent Y-axis
    
    plt.tight_layout()
    output_path = output_dir / 'scherrer_size_by_peak.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    
    # =========================================================================
    # Plot 2: Average Size by Concentration
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for idx, conc in enumerate(concentrations):
        if idx >= 4:
            break
        ax = axes[idx]
        
        conc_samples = [s for s in samples if s.concentration_ml == conc and s.result]
        conc_samples.sort(key=lambda x: x.time_hours)
        
        for peak, label, color in zip(peaks, peak_labels, colors):
            times = []
            sizes = []
            
            for s in conc_samples:
                size = get_size_for_peak(s.result, peak)
                if size and size > 0:
                    times.append(s.time_hours)
                    sizes.append(size)
            
            if times and sizes:
                ax.plot(times, sizes, 'o-', label=label, color=color, markersize=6)
        
        ax.set_xlabel('Plating Time (hours)')
        ax.set_ylabel('Crystallite Size (nm)')
        ax.set_title(f'Scherrer Size Evolution - {conc} ml Leveler')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, y_max_size)  # Consistent Y-axis
    
    plt.tight_layout()
    output_path = output_dir / 'scherrer_size_by_concentration.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    
    # =========================================================================
    # Plot 3: Size Comparison Heatmap
    # =========================================================================
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get unique times
    all_times = sorted(set(s.time_hours for s in samples if s.result))
    
    # Create data matrix for (111) peak
    data = np.full((len(concentrations), len(all_times)), np.nan)
    
    for i, conc in enumerate(concentrations):
        for s in samples:
            if s.concentration_ml == conc and s.result:
                size = get_size_for_peak(s.result, (1, 1, 1))
                if size and s.time_hours in all_times:
                    j = all_times.index(s.time_hours)
                    data[i, j] = size
    
    im = ax.imshow(data, aspect='auto', cmap='viridis')
    ax.set_xticks(range(len(all_times)))
    ax.set_xticklabels([f'{t:.1f}h' for t in all_times], rotation=45, ha='right')
    ax.set_yticks(range(len(concentrations)))
    ax.set_yticklabels([f'{c} ml' for c in concentrations])
    ax.set_xlabel('Plating Time')
    ax.set_ylabel('Leveler Concentration')
    ax.set_title('(111) Crystallite Size Heatmap (nm)')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Size (nm)')
    
    plt.tight_layout()
    output_path = output_dir / 'scherrer_size_heatmap_111.png'
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    
    return 3  # Number of plots generated


def main():
    """Main entry point."""
    print("=" * 60)
    print("AXCSAS Plot Generator")
    print("=" * 60)
    
    # Setup paths
    project_root = Path(__file__).parent
    data_dir = project_root / "data" / "raw" / "202511"
    output_dir = project_root / "outputs" / "plots"
    
    # Create output directories
    (output_dir / "fwhm").mkdir(parents=True, exist_ok=True)
    (output_dir / "scherrer").mkdir(parents=True, exist_ok=True)
    
    # Analyze all samples
    print(f"\nAnalyzing samples from: {data_dir}")
    samples = analyze_all_samples(data_dir)
    
    if not samples:
        print("No samples analyzed!")
        return 1
    
    print(f"\nSuccessfully analyzed {len(samples)} samples")
    
    # Generate FWHM plots
    print("\n" + "-" * 40)
    print("Generating FWHM Evolution Plots...")
    print("-" * 40)
    fwhm_count = plot_fwhm_evolution(samples, output_dir / "fwhm")
    
    # Generate Scherrer plots
    print("\n" + "-" * 40)
    print("Generating Scherrer Size Plots...")
    print("-" * 40)
    scherrer_count = plot_scherrer_size(samples, output_dir / "scherrer")
    
    print("\n" + "=" * 60)
    print(f"Complete! Generated {fwhm_count + scherrer_count} plots")
    print(f"Output directory: {output_dir}")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    exit(main())
