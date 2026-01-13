"""
Analysis Pipeline Module
=========================

High-level analysis orchestration for XRD data.
XRD 資料的高階分析協調模組。
"""

from axcsas.analysis.pipeline import (
    AXCSASPipeline,
    AnalysisConfig,
    PipelineResult,
    run_full_analysis,
    batch_analyze,
)
from axcsas.analysis.report_generator import (
    ComprehensiveResult,
    generate_comprehensive_report,
    generate_csv_summary,
)

# Alias for backward compatibility  
AnalysisPipeline = AXCSASPipeline
ReportGenerator = None  # Module doesn't have this class but keep for compatibility

__all__ = [
    "AXCSASPipeline",
    "AnalysisPipeline",
    "AnalysisConfig",
    "PipelineResult",
    "run_full_analysis",
    "batch_analyze",
    "ComprehensiveResult",
    "generate_comprehensive_report",
    "generate_csv_summary",
]
