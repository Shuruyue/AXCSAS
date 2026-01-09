#!/usr/bin/env python3
"""
Batch Analysis Script
Processes multiple XRD samples and generates comparative reports.

Usage:
    python batch_analysis.py --input-dir data/raw/202511/ --output-dir outputs/results/
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from analyze_sample import analyze_sample, load_config


def find_xrd_files(input_dir: str, extensions: list = None) -> list:
    """
    Find all XRD data files in directory.
    
    Args:
        input_dir: Directory to search
        extensions: File extensions to include
        
    Returns:
        List of file paths
    """
    extensions = extensions or ['.xy', '.csv', '.txt']
    input_path = Path(input_dir)
    
    files = []
    for ext in extensions:
        files.extend(input_path.glob(f'*{ext}'))
    
    return sorted(files)


def batch_analyze(
    input_dir: str,
    output_dir: str,
    config: dict = None,
    parallel: bool = True,
    max_workers: int = 4
):
    """
    Perform batch analysis on multiple samples.
    
    Args:
        input_dir: Directory containing XRD files
        output_dir: Output directory
        config: Configuration dictionary
        parallel: Use parallel processing
        max_workers: Maximum parallel workers
    """
    config = config or {}
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("AXCSAS Batch Analysis")
    print("=" * 60)
    
    # Find files
    files = find_xrd_files(input_dir)
    print(f"\nFound {len(files)} XRD files in {input_dir}")
    
    if not files:
        print("No XRD files found. Exiting.")
        return
    
    # Process files
    results = []
    errors = []
    
    print("\nProcessing samples...")
    print("-" * 60)
    
    for i, filepath in enumerate(files, 1):
        sample_name = filepath.stem
        print(f"\n[{i}/{len(files)}] {sample_name}")
        
        try:
            sample_output = output_path / sample_name
            result = analyze_sample(
                str(filepath),
                str(sample_output),
                config,
                generate_plots=True
            )
            results.append(result)
            print(f"      ✓ Completed")
        except Exception as e:
            errors.append({'sample': sample_name, 'error': str(e)})
            print(f"      ✗ Error: {e}")
    
    # =========================================================================
    # Generate Comparative Report
    # =========================================================================
    print("\n" + "=" * 60)
    print("Generating Comparative Report")
    print("=" * 60)
    
    # Compile results into DataFrame
    summary_data = []
    for r in results:
        row = {
            'sample': r['sample'],
            'peaks': r['peaks_detected'],
            'scherrer_size_nm': r['scherrer_size_nm']['average'],
            'wh_size_nm': r['williamson_hall']['size_nm'],
            'microstrain': r['williamson_hall']['microstrain'],
            'wh_r_squared': r['williamson_hall']['r_squared'],
        }
        
        if 'texture' in r:
            row['texture_degree'] = r['texture']['degree']
            row['is_random'] = r['texture']['is_random']
        
        summary_data.append(row)
    
    df = pd.DataFrame(summary_data)
    
    # Save summary CSV
    summary_csv = output_path / "batch_summary.csv"
    df.to_csv(summary_csv, index=False)
    print(f"✓ Summary CSV: {summary_csv}")
    
    # Generate comparison plots
    if len(results) > 1:
        # Crystallite size comparison
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Size comparison
        ax = axes[0, 0]
        x = range(len(df))
        width = 0.35
        ax.bar([i - width/2 for i in x], df['scherrer_size_nm'], width, 
               label='Scherrer', color='steelblue')
        ax.bar([i + width/2 for i in x], df['wh_size_nm'], width,
               label='W-H', color='coral')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Crystallite Size (nm)')
        ax.set_title('Crystallite Size Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels(df['sample'], rotation=45, ha='right')
        ax.legend()
        
        # Microstrain
        ax = axes[0, 1]
        ax.bar(df['sample'], df['microstrain'] * 100, color='green', alpha=0.7)
        ax.set_xlabel('Sample')
        ax.set_ylabel('Microstrain (%)')
        ax.set_title('Microstrain from W-H Analysis')
        ax.tick_params(axis='x', rotation=45)
        
        # Size trend
        ax = axes[1, 0]
        ax.plot(df['sample'], df['scherrer_size_nm'], 'o-', 
               label='Scherrer', color='steelblue')
        ax.plot(df['sample'], df['wh_size_nm'], 's-', 
               label='W-H', color='coral')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Crystallite Size (nm)')
        ax.set_title('Size Trend')
        ax.tick_params(axis='x', rotation=45)
        ax.legend()
        
        # W-H fit quality
        ax = axes[1, 1]
        colors = ['green' if r > 0.9 else 'orange' if r > 0.7 else 'red' 
                  for r in df['wh_r_squared']]
        ax.bar(df['sample'], df['wh_r_squared'], color=colors, alpha=0.7)
        ax.axhline(y=0.9, color='green', linestyle='--', alpha=0.5)
        ax.axhline(y=0.7, color='orange', linestyle='--', alpha=0.5)
        ax.set_xlabel('Sample')
        ax.set_ylabel('R²')
        ax.set_title('W-H Fit Quality')
        ax.set_ylim(0, 1.1)
        ax.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        comparison_plot = output_path / "batch_comparison.png"
        fig.savefig(comparison_plot, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"✓ Comparison plot: {comparison_plot}")
    
    # Error report
    if errors:
        error_csv = output_path / "batch_errors.csv"
        pd.DataFrame(errors).to_csv(error_csv, index=False)
        print(f"✓ Error log: {error_csv}")
    
    # Generate Markdown report
    report_md = output_path / "batch_report.md"
    with open(report_md, 'w', encoding='utf-8') as f:
        f.write(f"# AXCSAS Batch Analysis Report\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**Input Directory:** `{input_dir}`\n\n")
        f.write(f"---\n\n")
        
        f.write("## Summary Statistics\n\n")
        f.write(f"| Metric | Mean | Std | Min | Max |\n")
        f.write(f"|--------|------|-----|-----|-----|\n")
        f.write(f"| Scherrer Size (nm) | {df['scherrer_size_nm'].mean():.1f} | "
                f"{df['scherrer_size_nm'].std():.1f} | "
                f"{df['scherrer_size_nm'].min():.1f} | "
                f"{df['scherrer_size_nm'].max():.1f} |\n")
        f.write(f"| W-H Size (nm) | {df['wh_size_nm'].mean():.1f} | "
                f"{df['wh_size_nm'].std():.1f} | "
                f"{df['wh_size_nm'].min():.1f} | "
                f"{df['wh_size_nm'].max():.1f} |\n")
        f.write(f"| Microstrain | {df['microstrain'].mean():.2e} | "
                f"{df['microstrain'].std():.2e} | "
                f"{df['microstrain'].min():.2e} | "
                f"{df['microstrain'].max():.2e} |\n\n")
        
        f.write("## All Samples\n\n")
        f.write(df.to_markdown(index=False))
        f.write("\n\n")
        
        if errors:
            f.write("## Errors\n\n")
            for e in errors:
                f.write(f"- **{e['sample']}**: {e['error']}\n")
    
    print(f"✓ Report: {report_md}")
    
    print("\n" + "=" * 60)
    print(f"Batch Analysis Complete!")
    print(f"  Processed: {len(results)} / {len(files)} samples")
    if errors:
        print(f"  Errors: {len(errors)}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Batch analyze multiple XRD samples"
    )
    parser.add_argument(
        "--input-dir", "-i",
        required=True,
        help="Directory containing XRD data files"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="outputs/results",
        help="Output directory (default: outputs/results)"
    )
    parser.add_argument(
        "--config", "-c",
        default="config.yaml",
        help="Configuration file (default: config.yaml)"
    )
    
    args = parser.parse_args()
    
    config = load_config(args.config)
    batch_analyze(args.input_dir, args.output_dir, config)


if __name__ == "__main__":
    main()
