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
    dpi: int = 1000,
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
    x_param: str = "time",
    output_path: Optional[str] = None,
    dpi: int = 1000,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (15, 5.5),
) -> plt.Figure:
    """Plot TC evolution across samples with subplots by peak direction.
    繪製 TC 隨樣品參數演化圖（按峰方向分子圖）。
    
    Uses 3 subplots (one for each hkl), grouped by concentration.
    Similar structure to FWHM evolution plots.
    
    Args:
        data: List of dictionaries with keys:
            - 'name': Sample name
            - 'concentration' and 'time': Sample parameters
            - 'tc_values': Dict mapping hkl to TC value
        x_param: X-axis parameter ("time").
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        
    Returns:
        Matplotlib Figure object.
    """
    apply_axcsas_style()

    # Collect all unique hkl values
    all_hkls = set()
    for sample in data:
        all_hkls.update(sample.get('tc_values', {}).keys())

    hkl_list = sorted(list(all_hkls))
    n_hkls = len(hkl_list)

    if n_hkls == 0:
        raise ValueError("No TC data found in input")

    # Create subplots for each hkl (3 side-by-side)
    n_cols = min(n_hkls, 3)
    n_rows = (n_hkls + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    # X-axis label
    x_labels = {
        'time': 'Annealing Time (hours)',
        'concentration': 'Leveler Concentration (mL/1.5L)',
    }
    x_label = x_labels.get(x_param, x_param)

    # Group by concentration (when plotting vs time)
    concentrations = sorted(set(sample.get('concentration', 0) for sample in data))
    colors = get_color_palette(len(concentrations))
    conc_color_map = dict(zip(concentrations, colors))

    # Plot each hkl in its own subplot
    for idx, hkl in enumerate(hkl_list):
        ax = axes[idx]
        
        # Separate high-quality and low-quality data
        high_quality_data = [s for s in data if s.get('high_quality', True)]
        low_quality_data = [s for s in data if not s.get('high_quality', True)]
        
        # First, plot low-quality data as gray points (no lines)
        for sample in low_quality_data:
            x_val = sample.get(x_param, 0)
            tc_dict = sample.get('tc_values', {})
            if hkl in tc_dict:
                ax.scatter(x_val, tc_dict[hkl], 
                          c='gray', s=60, alpha=0.5, zorder=1,
                          marker='o', edgecolors='none')
        
        # Plot high-quality data with colored lines
        for conc in concentrations:
            x_values = []
            tc_values = []
            
            for sample in high_quality_data:
                if sample.get('concentration', 0) != conc:
                    continue
                x_val = sample.get(x_param, 0)
                tc_dict = sample.get('tc_values', {})
                if hkl in tc_dict:
                    x_values.append(x_val)
                    tc_values.append(tc_dict[hkl])
            
            if x_values:
                # Sort by x-value
                sorted_indices = np.argsort(x_values)
                x_sorted = np.array(x_values)[sorted_indices]
                tc_sorted = np.array(tc_values)[sorted_indices]
                
                color = conc_color_map[conc]
                label = f'{conc} mL/1.5L'
                
                ax.plot(x_sorted, tc_sorted, 'o-', color=color,
                       linewidth=2, markersize=8, 
                       markeredgecolor='black', markeredgewidth=0.5,
                       label=label, zorder=2)
        
        # Add legend entry for gray points if any exist
        if low_quality_data:
            ax.scatter([], [], c='gray', s=60, alpha=0.5, 
                      marker='o', label='R² < 0.995 (low quality)')
        
        # Random orientation reference
        ax.axhline(y=1.0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        
        # Random zone (0.9-1.1)
        ax.axhspan(0.9, 1.1, alpha=0.1, color='green', label='Random zone')
        
        ax.set_xlabel(x_label)
        ax.set_ylabel('Texture Coefficient (TC)')
        ax.set_title(f'{hkl} Peak', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_box_aspect(1)

    # Hide unused axes
    for idx in range(len(hkl_list), len(axes)):
        axes[idx].set_visible(False)

    # Set consistent Y-axis limits
    all_tc = []
    for sample in data:
        tc_dict = sample.get('tc_values', {})
        all_tc.extend(tc_dict.values())
    
    if all_tc:
        y_min = min(0.5, min(all_tc) * 0.9)
        y_max = max(1.5, max(all_tc) * 1.1)
        for ax in axes[:len(hkl_list)]:
            ax.set_ylim(y_min, y_max)

    fig.suptitle('Texture Coefficient Evolution by Peak', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig
