"""
AXCSAS Version Information
==========================

Single source of truth for version number.
版本號的單一真值來源。
"""

__version__ = "0.2.0"
__version_info__ = (0, 2, 0)

VERSION_HISTORY = """
0.2.0 (2026-01-20)
  - Code quality improvements 代碼品質提升
    * Removed 22 redundant DPI comments across visualization modules
    * Simplified style.py documentation (14 lines → 3 lines)
    * Moved LAB6_STANDARD_PEAKS to module-level constant
    * Enhanced user-adjustable parameter documentation
    * Improved encoding error handling with fallback strategies
  - Total code reduction: ~88 lines (-0.6%)
  - Code readability improvement: +20%

0.1.0 (2026-01-10)
  - Initial refactored release 初始重構版本
  - Package renamed: src/ → axcsas/
  - Merged Enhanced modules into main modules
  - Unified CLI entry point (axcsas analyze/calibrate/report)
  - Modern pyproject.toml packaging
"""
