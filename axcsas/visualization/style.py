"""AXCSAS Style Configuration
==========================

Unified matplotlib styling and colorblind-safe color palettes.
統一的 matplotlib 樣式配置與色盲友善調色盤。
"""

from typing import Any, Dict, List

import matplotlib.pyplot as plt

# =============================================================================
# AXCSAS Standard Style
# =============================================================================

AXCSAS_STYLE: Dict[str, Any] = {
    # Figure settings
    'figure.figsize': (10, 6),
    'figure.dpi': 100,
    'figure.facecolor': 'white',
    'figure.edgecolor': 'white',

    # Save settings: 2400 DPI for journal publication
    # Reference: Nature/Science guidelines (300-600 DPI minimum)
    'savefig.dpi': 2400,
    'savefig.facecolor': 'white',
    'savefig.edgecolor': 'white',
    'savefig.bbox': 'tight',

    # Font settings
    'font.family': 'sans-serif',
    'font.size': 12,

    # Axes settings
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'axes.titleweight': 'bold',
    'axes.grid': True,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 1.2,

    # Grid settings
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'grid.linewidth': 0.8,

    # Legend settings
    'legend.fontsize': 11,
    'legend.framealpha': 0.9,
    'legend.edgecolor': 'gray',
    'legend.fancybox': True,

    # Tick settings
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'xtick.major.size': 5,
    'ytick.major.size': 5,

    # Line settings
    'lines.linewidth': 2,
    'lines.markersize': 8,
}


# =============================================================================
# Colorblind-Safe Color Palette (Wong, 2011)
# =============================================================================

COLORBLIND_SAFE: List[str] = [
    '#0077BB',  # Blue - primary
    '#EE7733',  # Orange - secondary
    '#009988',  # Teal
    '#CC3311',  # Red
    '#EE3377',  # Magenta
    '#33BBEE',  # Cyan
    '#BBBBBB',  # Gray
    '#000000',  # Black
]

# Semantic color mapping for XRD peaks
PEAK_COLORS: Dict[str, str] = {
    '(111)': '#0077BB',  # Blue
    '(200)': '#EE7733',  # Orange
    '(220)': '#009988',  # Teal
    '(311)': '#CC3311',  # Red
    '(222)': '#EE3377',  # Magenta
}


# =============================================================================
# Style Application Functions
# =============================================================================

def apply_axcsas_style() -> None:
    """Apply AXCSAS style to all subsequent matplotlib plots.
    將 AXCSAS 樣式套用至後續所有 matplotlib 圖表。
    
    Example:
        >>> from axcsas.visualization.style import apply_axcsas_style
        >>> apply_axcsas_style()
        >>> # All subsequent plots will use AXCSAS style

    """
    plt.rcParams.update(AXCSAS_STYLE)


def get_color_palette(n_colors: int = 6) -> List[str]:
    """Get a colorblind-safe color palette.
    取得色盲友善調色盤。
    
    Args:
        n_colors: Number of colors needed (max 8).
        
    Returns:
        List of hex color strings.
        
    Example:
        >>> colors = get_color_palette(3)
        >>> print(colors)
        ['#0077BB', '#EE7733', '#009988']

    """
    return COLORBLIND_SAFE[:min(n_colors, len(COLORBLIND_SAFE))]


def get_peak_color(hkl: str) -> str:
    """Get color for a specific (hkl) peak.
    取得特定 (hkl) 峰的顏色。
    
    Args:
        hkl: Peak identifier string, e.g. "(111)".
        
    Returns:
        Hex color string. Returns gray if not found.

    """
    return PEAK_COLORS.get(hkl, '#BBBBBB')


def create_figure(
    figsize: tuple = None,
    apply_style: bool = True
) -> tuple:
    """Create a figure with AXCSAS styling.
    建立套用 AXCSAS 樣式的圖表。
    
    Args:
        figsize: Optional figure size (width, height) in inches.
        apply_style: Whether to apply AXCSAS style.
        
    Returns:
        Tuple of (fig, ax) matplotlib objects.

    """
    if apply_style:
        apply_axcsas_style()

    if figsize is None:
        figsize = AXCSAS_STYLE['figure.figsize']

    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def save_figure(
    fig,
    filepath: str,
    dpi: int = 2400,
    format: str = None,
    transparent: bool = False
) -> None:
    """Save figure with standardized settings.
    使用標準化設定儲存圖表。
    
    Args:
        fig: Matplotlib figure object.
        filepath: Output file path.
        dpi: Resolution (dots per inch).
        format: Output format (png, svg, pdf). Auto-detected from filepath if None.
        transparent: Whether to use transparent background.

    """
    fig.savefig(
        filepath,
        dpi=dpi,
        format=format,
        transparent=transparent,
        bbox_inches='tight',
        facecolor='white' if not transparent else 'none',
        edgecolor='none'
    )
