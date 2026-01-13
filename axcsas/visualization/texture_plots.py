"""Texture Coefficient Visualization Module
========================================

Plots for Harris Texture Coefficient (TC) analysis.
Harris 紋理係數（TC）分析繪圖模組。
"""

from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .style import (
    COLORBLIND_SAFE,
    apply_axcsas_style,
    get_color_palette,
    save_figure,
)


def plot_texture_polar(
    tc_values: Dict[str, float],
    output_path: Optional[str] = None,
    dpi: int = 300,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (8, 8),
    sample_name: str = "Sample",
    show_random_circle: bool = True,
) -> plt.Figure:
    """Plot texture coefficients as a polar (radar) chart.
    以極座標雷達圖繪製紋理係數。
    
    TC = 1 indicates random orientation (powder average).
    TC > 1 indicates preferred orientation.
    TC < 1 indicates under-represented direction.
    
    Args:
        tc_values: Dictionary mapping hkl strings to TC values.
                   e.g., {'(111)': 1.25, '(200)': 0.85, '(220)': 1.10}
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        sample_name: Name for title.
        show_random_circle: Whether to show TC=1 reference circle.
        
    Returns:
        Matplotlib Figure object.
        
    Example:
        >>> tc = {'(111)': 1.35, '(200)': 0.72, '(220)': 0.93}
        >>> fig = plot_texture_polar(tc)

    """
    apply_axcsas_style()

    hkl_labels = list(tc_values.keys())
    tc_vals = list(tc_values.values())
    n_peaks = len(hkl_labels)

    if n_peaks < 2:
        raise ValueError("Need at least 2 peaks for polar plot")

    # Compute angle for each peak
    angles = np.linspace(0, 2 * np.pi, n_peaks, endpoint=False).tolist()

    # Close the polygon
    tc_vals_closed = tc_vals + [tc_vals[0]]
    angles_closed = angles + [angles[0]]

    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))

    # Plot TC polygon
    ax.plot(angles_closed, tc_vals_closed, 'o-', linewidth=2.5,
           color=COLORBLIND_SAFE[0], markersize=12, label='TC values')
    ax.fill(angles_closed, tc_vals_closed, alpha=0.25, color=COLORBLIND_SAFE[0])

    # Random orientation reference circle (TC = 1)
    if show_random_circle:
        random_circle = [1.0] * (n_peaks + 1)
        ax.plot(angles_closed, random_circle, '--', linewidth=2,
               color='red', alpha=0.7, label='Random (TC=1)')

    # Configure labels
    ax.set_xticks(angles)
    ax.set_xticklabels(hkl_labels, fontsize=12, fontweight='bold')

    # Configure radial axis
    max_tc = max(tc_vals) * 1.2
    ax.set_ylim(0, max(max_tc, 1.5))
    ax.set_ylabel('TC', labelpad=30)

    ax.set_title(f'Texture Coefficient Analysis\n{sample_name}',
                fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

    # Add texture interpretation
    is_random = all(0.9 <= tc <= 1.1 for tc in tc_vals)
    dominant = max(tc_values.items(), key=lambda x: x[1])

    if is_random:
        interp_text = "Random texture"
        interp_color = 'green'
    else:
        interp_text = f"Preferred: {dominant[0]} (TC={dominant[1]:.2f})"
        interp_color = 'orange'

    fig.text(0.5, 0.02, interp_text, ha='center', fontsize=12,
            fontweight='bold', color=interp_color)

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_tc_evolution(
    data: List[Dict[str, Any]],
    x_param: str = "concentration",
    output_path: Optional[str] = None,
    dpi: int = 300,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (12, 6),
) -> plt.Figure:
    """Plot TC evolution across samples.
    繪製 TC 隨樣品參數演化圖。
    
    Args:
        data: List of dictionaries with keys:
            - 'name': Sample name
            - 'concentration' or 'time': X-axis value
            - 'tc_values': Dict mapping hkl to TC value
        x_param: X-axis parameter ("concentration" or "time").
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        
    Returns:
        Matplotlib Figure object.
        
    Example:
        >>> data = [
        ...     {'concentration': 0, 'tc_values': {'(111)': 1.2, '(200)': 0.8}},
        ...     {'concentration': 9, 'tc_values': {'(111)': 1.4, '(200)': 0.6}},
        ... ]
        >>> fig = plot_tc_evolution(data, x_param='concentration')

    """
    apply_axcsas_style()

    # Collect all unique hkl values
    all_hkls = set()
    for sample in data:
        all_hkls.update(sample.get('tc_values', {}).keys())

    hkl_list = sorted(list(all_hkls))

    if not hkl_list:
        raise ValueError("No TC data found in input")

    fig, ax = plt.subplots(figsize=figsize)

    # X-axis label
    x_labels = {
        'concentration': 'Leveler Concentration (mL/L)',
        'time': 'Plating Time (hours)',
        'age': 'Sample Age (hours)',
    }
    x_label = x_labels.get(x_param, x_param)

    colors = get_color_palette(len(hkl_list))

    # Plot each hkl direction
    for idx, hkl in enumerate(hkl_list):
        x_values = []
        tc_values = []

        for sample in data:
            x_val = sample.get(x_param, sample.get('concentration', 0))
            tc_dict = sample.get('tc_values', {})
            if hkl in tc_dict:
                x_values.append(x_val)
                tc_values.append(tc_dict[hkl])

        if x_values:
            ax.plot(x_values, tc_values, 'o-', color=colors[idx],
                   linewidth=2, markersize=8, label=hkl)

    # Random orientation reference
    ax.axhline(y=1.0, color='red', linestyle='--', linewidth=2,
              alpha=0.7, label='Random (TC=1)')

    # Random zone (0.9-1.1)
    ax.axhspan(0.9, 1.1, alpha=0.1, color='green', label='Random zone')

    ax.set_xlabel(x_label)
    ax.set_ylabel('Texture Coefficient (TC)')
    ax.set_title('Texture Coefficient Evolution')
    ax.legend(loc='best', fontsize=10)

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig
