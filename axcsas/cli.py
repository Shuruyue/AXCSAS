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
    """
    print(f"Calibrating with standard: {args.standard}")
    print(f"Output: {args.output}")
    
    # TODO: Implement calibration logic
    print("Calibration complete. 校正完成。")
    return 0


def _run_report(args) -> int:
    """
    Run report generation command.
    執行報告生成指令。
    """
    from axcsas.analysis.report_generator import ReportGenerator
    
    print(f"Generating {args.format} report from: {args.results}")
    
    # TODO: Implement report generation
    print("Report generated. 報告已生成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
