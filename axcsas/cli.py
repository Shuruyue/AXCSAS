#!/usr/bin/env python3
"""
AXCSAS Command Line Interface
=============================

Unified entry point for all analysis operations.
所有分析操作的統一入口點。
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List

from axcsas.__version__ import __version__


def main(argv: Optional[List[str]] = None) -> int:
    """
    Main CLI entry point.
    主 CLI 入口點。
    """
    parser = argparse.ArgumentParser(
        prog="axcsas",
        description="Advanced XRD Crystallite Size Analysis System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  axcsas analyze data/sample.txt -o outputs/
  axcsas calibrate data/LaB6_standard.txt
  axcsas report outputs/results.csv
        """
    )
    
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands 可用指令")
    
    # analyze command
    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Analyze XRD data 分析 XRD 資料"
    )
    analyze_parser.add_argument(
        "input",
        type=Path,
        help="Input file or directory 輸入檔案或目錄"
    )
    analyze_parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("outputs"),
        help="Output directory 輸出目錄"
    )
    analyze_parser.add_argument(
        "-c", "--config",
        type=Path,
        help="Config file path 配置檔案路徑"
    )
    analyze_parser.add_argument(
        "--batch",
        action="store_true",
        help="Batch processing mode 批次處理模式"
    )
    
    # calibrate command
    cal_parser = subparsers.add_parser(
        "calibrate",
        help="Calibrate instrument 校準儀器"
    )
    cal_parser.add_argument(
        "standard",
        type=Path,
        help="Standard material data file 標準樣品資料檔案"
    )
    cal_parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("calibration.yaml"),
        help="Output calibration file 輸出校正檔案"
    )
    
    # report command
    report_parser = subparsers.add_parser(
        "report",
        help="Generate analysis report 產生分析報告"
    )
    report_parser.add_argument(
        "results",
        type=Path,
        help="Results CSV file 結果 CSV 檔案"
    )
    report_parser.add_argument(
        "-f", "--format",
        choices=["markdown", "html", "pdf"],
        default="markdown",
        help="Report format 報告格式"
    )
    
    args = parser.parse_args(argv)
    
    if args.command is None:
        parser.print_help()
        return 0
    
    try:
        if args.command == "analyze":
            return _run_analyze(args)
        elif args.command == "calibrate":
            return _run_calibrate(args)
        elif args.command == "report":
            return _run_report(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    return 0


def _run_analyze(args) -> int:
    """
    Run analysis command.
    執行分析指令。
    """
    from axcsas.core.config_loader import load_config
    from axcsas.analysis.pipeline import AnalysisPipeline
    
    print(f"Loading configuration... 載入配置...")
    config = load_config(args.config)
    
    print(f"Analyzing: {args.input}")
    print(f"Output to: {args.output}")
    
    # Ensure output directory exists
    args.output.mkdir(parents=True, exist_ok=True)
    
    pipeline = AnalysisPipeline(config)
    
    if args.batch or args.input.is_dir():
        # Batch processing
        files = list(args.input.glob("*.txt")) if args.input.is_dir() else [args.input]
        for f in files:
            print(f"  Processing: {f.name}")
            pipeline.process_file(f, args.output)
    else:
        pipeline.process_file(args.input, args.output)
    
    print("Analysis complete. 分析完成。")
    return 0


def _run_calibrate(args) -> int:
    """
    Run calibration command.
    執行校正指令。
    
    Calibrates Caglioti parameters (U, V, W) using a standard material
    like LaB6 (NIST SRM 660c) or Si.
    使用標準樣品（如 LaB6 或 Si）校準 Caglioti 參數。
    """
    import numpy as np
    import yaml
    from axcsas.analysis.pipeline import load_bruker_txt, find_peak_in_range
    from axcsas.methods.caglioti import CagliotiCorrection
    
    print(f"Calibrating with standard: {args.standard}")
    
    # Load standard data
    try:
        two_theta, intensity = load_bruker_txt(str(args.standard))
    except Exception as e:
        print(f"Error loading file: {e}")
        return 1
    
    if len(two_theta) == 0:
        print("Error: No data found in file")
        return 1
    
    print(f"  Loaded {len(two_theta)} data points")
    print(f"  2θ range: {two_theta.min():.1f}° - {two_theta.max():.1f}°")
    
    # LaB6 standard peak positions (NIST SRM 660c)
    LAB6_PEAKS = {
        (1, 0, 0): 21.358,
        (1, 1, 0): 30.385,
        (1, 1, 1): 37.442,
        (2, 0, 0): 43.505,
        (2, 1, 0): 48.958,
        (2, 1, 1): 53.993,
        (2, 2, 0): 63.218,
        (3, 0, 0): 67.548,
        (3, 1, 0): 71.704,
        (3, 1, 1): 75.724,
    }
    
    # Find peaks and measure FWHM
    peaks_data = []
    print("\nFitting standard peaks...")
    
    for hkl, expected_pos in LAB6_PEAKS.items():
        if two_theta.min() <= expected_pos <= two_theta.max():
            peak = find_peak_in_range(two_theta, intensity, expected_pos, window=2.0)
            if peak is not None and peak.fwhm > 0:
                peaks_data.append({
                    'hkl': hkl,
                    'two_theta': peak.two_theta,
                    'fwhm': peak.fwhm,
                    'intensity': peak.intensity,
                })
                print(f"  ✓ {hkl}: 2θ={peak.two_theta:.3f}°, FWHM={peak.fwhm:.4f}°")
    
    if len(peaks_data) < 3:
        print(f"\nWarning: Only {len(peaks_data)} peaks found. Need at least 3 for calibration.")
        print("Using default Caglioti parameters.")
        U, V, W = 0.001, -0.001, 0.0025
    else:
        # Fit Caglioti: FWHM² = U·tan²θ + V·tanθ + W
        theta_rad = np.array([np.radians(p['two_theta'] / 2) for p in peaks_data])
        tan_theta = np.tan(theta_rad)
        fwhm_sq = np.array([p['fwhm'] ** 2 for p in peaks_data])
        
        # Linear regression: y = U·x² + V·x + W
        X = np.column_stack([tan_theta**2, tan_theta, np.ones_like(tan_theta)])
        coeffs, residuals, rank, s = np.linalg.lstsq(X, fwhm_sq, rcond=None)
        U, V, W = coeffs
        
        print(f"\nCaglioti parameters fitted from {len(peaks_data)} peaks:")
    
    print(f"  U = {U:.6f}")
    print(f"  V = {V:.6f}")
    print(f"  W = {W:.6f}")
    print(f"  FWHM at 43° ≈ {np.sqrt(U * np.tan(np.radians(21.5))**2 + V * np.tan(np.radians(21.5)) + W):.4f}°")
    
    # Save calibration
    calibration = {
        'caglioti': {
            'U': float(U),
            'V': float(V),
            'W': float(W),
        },
        'standard': str(args.standard),
        'n_peaks_used': len(peaks_data),
    }
    
    with open(args.output, 'w') as f:
        yaml.dump(calibration, f, default_flow_style=False)
    
    print(f"\nCalibration saved to: {args.output}")
    print("Calibration complete. 校正完成。")
    return 0


def _run_report(args) -> int:
    """
    Run report generation command.
    執行報告生成指令。
    
    Generates analysis report from results CSV file.
    從結果 CSV 檔案生成分析報告。
    """
    import pandas as pd
    from datetime import datetime
    from axcsas.analysis.report_generator import generate_comprehensive_report
    
    print(f"Generating {args.format} report from: {args.results}")
    
    # Load results
    if not args.results.exists():
        print(f"Error: Results file not found: {args.results}")
        return 1
    
    try:
        df = pd.read_csv(args.results)
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return 1
    
    print(f"  Loaded {len(df)} analysis records")
    
    # Generate report content
    report_lines = [
        "# AXCSAS Analysis Report",
        f"# AXCSAS 分析報告",
        "",
        f"**Generated / 生成時間**: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "---",
        "",
        "## Summary / 摘要",
        "",
    ]
    
    # Add summary statistics
    if 'size_nm' in df.columns:
        valid_sizes = df['size_nm'].dropna()
        if len(valid_sizes) > 0:
            report_lines.extend([
                f"- **Total samples / 總樣品數**: {len(df)}",
                f"- **Average crystallite size / 平均晶粒尺寸**: {valid_sizes.mean():.1f} nm",
                f"- **Size range / 尺寸範圍**: {valid_sizes.min():.1f} - {valid_sizes.max():.1f} nm",
                "",
            ])
    
    # Add data table
    report_lines.extend([
        "## Data Table / 資料表",
        "",
        df.to_markdown(index=False) if hasattr(df, 'to_markdown') else df.to_string(),
        "",
    ])
    
    report_content = "\n".join(report_lines)
    
    # Determine output path
    output_path = args.results.parent / f"report_{args.results.stem}.{args.format}"
    if args.format == "markdown":
        output_path = output_path.with_suffix(".md")
    
    # Save report
    if args.format == "markdown":
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report_content)
        print(f"\nMarkdown report saved to: {output_path}")
    elif args.format == "html":
        try:
            import markdown
            html_content = markdown.markdown(report_content, extensions=['tables'])
            html_output = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>AXCSAS Report</title>
    <style>
        body {{ font-family: 'Segoe UI', sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background: #4a90d9; color: white; }}
    </style>
</head>
<body>
{html_content}
</body>
</html>"""
            output_path = output_path.with_suffix(".html")
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html_output)
            print(f"\nHTML report saved to: {output_path}")
        except ImportError:
            print("Warning: 'markdown' package not installed. Saving as markdown instead.")
            output_path = output_path.with_suffix(".md")
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(report_content)
    else:
        print(f"Warning: {args.format} format not yet implemented. Saving as markdown.")
        output_path = output_path.with_suffix(".md")
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report_content)
    
    print("Report generated. 報告已生成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
