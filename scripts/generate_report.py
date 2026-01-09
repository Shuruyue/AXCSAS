#!/usr/bin/env python3
"""
Report Generation Script
Generates comprehensive analysis reports in Markdown and PDF formats.

Usage:
    python generate_report.py --input outputs/results/batch_summary.csv --output outputs/reports/
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_batch_results(results_dir: str) -> dict:
    """
    Load batch analysis results.
    
    Args:
        results_dir: Directory containing analysis results
        
    Returns:
        Dictionary with summary data and individual samples
    """
    results_path = Path(results_dir)
    
    data = {
        'summary': None,
        'samples': []
    }
    
    # Load batch summary if exists
    summary_csv = results_path / "batch_summary.csv"
    if summary_csv.exists():
        data['summary'] = pd.read_csv(summary_csv)
    
    # Load individual sample results
    for sample_dir in results_path.iterdir():
        if sample_dir.is_dir():
            summary_file = sample_dir / f"{sample_dir.name}_summary.yaml"
            if summary_file.exists():
                with open(summary_file, 'r') as f:
                    sample_data = yaml.safe_load(f)
                    data['samples'].append(sample_data)
    
    return data


def generate_markdown_report(
    data: dict,
    output_path: Path,
    title: str = "XRD Analysis Report"
) -> str:
    """
    Generate comprehensive Markdown report.
    
    Args:
        data: Analysis data dictionary
        output_path: Output file path
        title: Report title
        
    Returns:
        Path to generated report
    """
    report_file = output_path / "analysis_report.md"
    
    with open(report_file, 'w', encoding='utf-8') as f:
        # Header
        f.write(f"# {title}\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**Analysis System:** AXCSAS v0.1.0\n\n")
        f.write("---\n\n")
        
        # Executive Summary
        f.write("## Executive Summary\n\n")
        
        if data['summary'] is not None:
            df = data['summary']
            n_samples = len(df)
            
            f.write(f"This report presents XRD crystallite size analysis for **{n_samples} samples**.\n\n")
            
            # Key findings
            f.write("### Key Findings\n\n")
            
            avg_scherrer = df['scherrer_size_nm'].mean()
            avg_wh = df['wh_size_nm'].mean()
            avg_strain = df['microstrain'].mean()
            
            f.write(f"- **Average Crystallite Size (Scherrer):** {avg_scherrer:.1f} ± {df['scherrer_size_nm'].std():.1f} nm\n")
            f.write(f"- **Average Crystallite Size (W-H):** {avg_wh:.1f} ± {df['wh_size_nm'].std():.1f} nm\n")
            f.write(f"- **Average Microstrain:** {avg_strain:.2e}\n\n")
        
        f.write("---\n\n")
        
        # Methodology
        f.write("## Methodology\n\n")
        f.write("### Data Processing Pipeline\n\n")
        f.write("1. **Data Loading:** Bruker TXT/XY format\n")
        f.write("2. **Preprocessing:** Savitzky-Golay smoothing, Chebyshev background subtraction\n")
        f.write("3. **Peak Fitting:** Pseudo-Voigt profile, Levenberg-Marquardt optimization\n")
        f.write("4. **Size Calculation:** Scherrer equation (K=0.89)\n")
        f.write("5. **Strain Analysis:** Williamson-Hall method\n\n")
        
        # Equations
        f.write("### Equations\n\n")
        f.write("**Scherrer Equation:**\n")
        f.write("$$D = \\frac{K\\lambda}{\\beta\\cos\\theta}$$\n\n")
        f.write("**Williamson-Hall:**\n")
        f.write("$$\\beta\\cos\\theta = \\frac{K\\lambda}{D} + 4\\varepsilon\\sin\\theta$$\n\n")
        
        f.write("---\n\n")
        
        # Results Table
        f.write("## Results\n\n")
        
        if data['summary'] is not None:
            df = data['summary']
            
            f.write("### Summary Table\n\n")
            f.write("| Sample | Scherrer (nm) | W-H Size (nm) | Microstrain | W-H R² |\n")
            f.write("|--------|---------------|---------------|-------------|--------|\n")
            
            for _, row in df.iterrows():
                f.write(f"| {row['sample']} | {row['scherrer_size_nm']:.1f} | "
                       f"{row['wh_size_nm']:.1f} | {row['microstrain']:.2e} | "
                       f"{row['wh_r_squared']:.3f} |\n")
            
            f.write("\n")
            
            # Statistics
            f.write("### Statistical Summary\n\n")
            f.write("| Metric | Mean | Std Dev | Min | Max |\n")
            f.write("|--------|------|---------|-----|-----|\n")
            
            for col in ['scherrer_size_nm', 'wh_size_nm', 'microstrain']:
                col_label = col.replace('_', ' ').title()
                f.write(f"| {col_label} | {df[col].mean():.2f} | {df[col].std():.2f} | "
                       f"{df[col].min():.2f} | {df[col].max():.2f} |\n")
            
            f.write("\n")
        
        f.write("---\n\n")
        
        # Discussion
        f.write("## Discussion\n\n")
        
        if data['summary'] is not None:
            df = data['summary']
            
            # Size comparison
            scherrer_mean = df['scherrer_size_nm'].mean()
            wh_mean = df['wh_size_nm'].mean()
            diff_pct = abs(scherrer_mean - wh_mean) / scherrer_mean * 100
            
            f.write("### Crystallite Size Analysis\n\n")
            f.write(f"The Scherrer and Williamson-Hall methods yield crystallite sizes "
                   f"differing by **{diff_pct:.1f}%**. ")
            
            if wh_mean > scherrer_mean:
                f.write("The W-H method gives larger sizes, indicating significant microstrain "
                       "contribution to peak broadening.\n\n")
            else:
                f.write("The methods show good agreement, suggesting minimal strain broadening.\n\n")
            
            # Quality assessment
            avg_r2 = df['wh_r_squared'].mean()
            f.write("### Fit Quality\n\n")
            if avg_r2 > 0.95:
                f.write(f"The average W-H fit R² = {avg_r2:.3f} indicates **excellent** linear correlation.\n\n")
            elif avg_r2 > 0.85:
                f.write(f"The average W-H fit R² = {avg_r2:.3f} indicates **good** linear correlation.\n\n")
            else:
                f.write(f"The average W-H fit R² = {avg_r2:.3f} suggests potential non-uniform strain "
                       "or measurement issues.\n\n")
        
        f.write("---\n\n")
        
        # Appendix
        f.write("## Appendix\n\n")
        f.write("### Analysis Parameters\n\n")
        f.write("| Parameter | Value | Reference |\n")
        f.write("|-----------|-------|----------|\n")
        f.write("| X-ray Wavelength (λ) | 1.54056 Å | Cu Kα₁ |\n")
        f.write("| Scherrer Constant (K) | 0.89 | Spherical grains |\n")
        f.write("| Smoothing | Savitzky-Golay | Window=11, Order=3 |\n")
        f.write("| Background | Chebyshev | Degree=5 |\n\n")
        
        # References
        f.write("### References\n\n")
        f.write("1. Langford, J. I., & Wilson, A. J. C. (1978). *Scherrer after sixty years*. "
               "J. Appl. Cryst., 11, 102-113.\n")
        f.write("2. Williamson, G. K., & Hall, W. H. (1953). *X-ray line broadening from filed "
               "aluminium and wolfram*. Acta Metallurgica, 1(1), 22-31.\n")
    
    return str(report_file)


def generate_plots(
    data: dict,
    output_path: Path
) -> list:
    """
    Generate analysis plots.
    
    Args:
        data: Analysis data
        output_path: Output directory
        
    Returns:
        List of generated plot paths
    """
    plots = []
    
    if data['summary'] is None:
        return plots
    
    df = data['summary']
    
    # 1. Size comparison plot
    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(df))
    width = 0.35
    
    bars1 = ax.bar([i - width/2 for i in x], df['scherrer_size_nm'], width,
                   label='Scherrer', color='steelblue', edgecolor='black')
    bars2 = ax.bar([i + width/2 for i in x], df['wh_size_nm'], width,
                   label='W-H', color='coral', edgecolor='black')
    
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('Crystallite Size (nm)', fontsize=12)
    ax.set_title('Crystallite Size Comparison: Scherrer vs Williamson-Hall', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(df['sample'], rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plot1 = output_path / "size_comparison.png"
    fig.savefig(plot1, dpi=300, bbox_inches='tight')
    plt.close(fig)
    plots.append(str(plot1))
    
    # 2. Size vs Strain correlation
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(df['wh_size_nm'], df['microstrain'] * 1e3,
                        c=df['wh_r_squared'], cmap='viridis',
                        s=100, edgecolors='black', alpha=0.8)
    
    ax.set_xlabel('Crystallite Size (nm)', fontsize=12)
    ax.set_ylabel('Microstrain (×10⁻³)', fontsize=12)
    ax.set_title('Size-Strain Relationship', fontsize=14)
    
    cbar = plt.colorbar(scatter)
    cbar.set_label('W-H R²')
    
    plt.tight_layout()
    plot2 = output_path / "size_strain_correlation.png"
    fig.savefig(plot2, dpi=300, bbox_inches='tight')
    plt.close(fig)
    plots.append(str(plot2))
    
    # 3. Distribution histograms
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].hist(df['scherrer_size_nm'], bins=10, color='steelblue',
                 edgecolor='black', alpha=0.7)
    axes[0].axvline(df['scherrer_size_nm'].mean(), color='red',
                    linestyle='--', linewidth=2, label='Mean')
    axes[0].set_xlabel('Crystallite Size (nm)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Scherrer Size Distribution')
    axes[0].legend()
    
    axes[1].hist(df['microstrain'] * 1e3, bins=10, color='green',
                 edgecolor='black', alpha=0.7)
    axes[1].axvline(df['microstrain'].mean() * 1e3, color='red',
                    linestyle='--', linewidth=2, label='Mean')
    axes[1].set_xlabel('Microstrain (×10⁻³)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Microstrain Distribution')
    axes[1].legend()
    
    plt.tight_layout()
    plot3 = output_path / "distributions.png"
    fig.savefig(plot3, dpi=300, bbox_inches='tight')
    plt.close(fig)
    plots.append(str(plot3))
    
    return plots


def generate_report(
    input_dir: str,
    output_dir: str,
    title: str = "XRD Crystallite Size Analysis Report"
):
    """
    Generate complete analysis report.
    
    Args:
        input_dir: Directory containing analysis results
        output_dir: Output directory for reports
        title: Report title
    """
    print("=" * 60)
    print("AXCSAS Report Generation")
    print("=" * 60)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("\n[1/3] Loading analysis results...")
    data = load_batch_results(input_dir)
    
    if data['summary'] is not None:
        print(f"      Found {len(data['summary'])} samples")
    else:
        print("      Warning: No batch summary found")
    
    # Generate plots
    print("\n[2/3] Generating plots...")
    plots = generate_plots(data, output_path)
    for p in plots:
        print(f"      ✓ {Path(p).name}")
    
    # Generate report
    print("\n[3/3] Generating report...")
    report = generate_markdown_report(data, output_path, title)
    print(f"      ✓ {Path(report).name}")
    
    print("\n" + "=" * 60)
    print("Report Generation Complete!")
    print(f"Output: {output_dir}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Generate XRD analysis report"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Directory containing analysis results"
    )
    parser.add_argument(
        "--output", "-o",
        default="outputs/reports",
        help="Output directory (default: outputs/reports)"
    )
    parser.add_argument(
        "--title", "-t",
        default="XRD Crystallite Size Analysis Report",
        help="Report title"
    )
    
    args = parser.parse_args()
    
    generate_report(args.input, args.output, args.title)


if __name__ == "__main__":
    main()
