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
# report_generator.py has been removed - no longer available


# Alias for backward compatibility  
AnalysisPipeline = AXCSASPipeline


class _DeprecatedReportGenerator:
    """
    Deprecated placeholder for backward compatibility.
    已棄用的向後相容佔位類別。
    """
    def __init__(self, *args, **kwargs):
        import warnings
        warnings.warn(
            "ReportGenerator is deprecated. Use generate_comprehensive_report() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        raise NotImplementedError(
            "ReportGenerator has been replaced. "
            "Please use generate_comprehensive_report() function."
        )


ReportGenerator = _DeprecatedReportGenerator

__all__ = [
    "AXCSASPipeline",
    "AnalysisPipeline",
    "AnalysisConfig",
    "PipelineResult",
    "run_full_analysis",
    "batch_analyze",
]
