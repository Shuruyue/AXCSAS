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
    x_param: str = "time",
    output_path: Optional[str] = None,
    dpi: int = 2400,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (15, 5.5),
    normalize: bool = False,  # NEW: If True, plot as fraction (sum=1)
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
        normalize: If True, convert TC to fraction (f = TC/3, sum=1).
                   If False, plot as standard TC (sum=3).
        
        
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
    
    # Different line styles and markers for different concentrations
    line_styles_list = ['-', '--', '-.', ':']
    markers_list = ['o', 's', '^', 'D']
    
    conc_linestyle_map = {
        conc: line_styles_list[i % len(line_styles_list)]
        for i, conc in enumerate(concentrations)
    }
    conc_marker_map = {
        conc: markers_list[i % len(markers_list)]
        for i, conc in enumerate(concentrations)
    }

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
                # Normalize to fraction if requested
                tc_val = tc_dict[hkl] / 3.0 if normalize else tc_dict[hkl]
                ax.scatter(x_val, tc_val, 
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
                    # Normalize to fraction if requested
                    tc_val = tc_dict[hkl] / 3.0 if normalize else tc_dict[hkl]
                    tc_values.append(tc_val)
            
            if x_values:
                # Sort by x-value
                sorted_indices = np.argsort(x_values)
                x_sorted = np.array(x_values)[sorted_indices]
                tc_sorted = np.array(tc_values)[sorted_indices]
                
                color = conc_color_map[conc]
                linestyle = conc_linestyle_map[conc]
                marker = conc_marker_map[conc]
                label = f'{conc} mL/1.5L'
                
                ax.plot(x_sorted, tc_sorted, 
                       linestyle=linestyle,
                       marker=marker,
                       color=color,
                       linewidth=2, 
                       markersize=7,
                       markeredgecolor='black', 
                       markeredgewidth=0.5,
                       label=label, 
                       zorder=2,
                       alpha=0.9)
        
        # Add legend entry for gray points if any exist
        if low_quality_data:
            ax.scatter([], [], c='gray', s=60, alpha=0.5, 
                      marker='o', label='R² < 0.95 (low quality)')
        
        # Random orientation reference and zone (adjusted for normalization)
        if normalize:
            # For fraction: random = 1/3 ≈ 0.333
            random_ref = 1.0 / 3.0
            random_min = 0.9 / 3.0
            random_max = 1.1 / 3.0
            y_label = 'Texture Fraction (f)'
        else:
            # For standard TC: random = 1.0
            random_ref = 1.0
            random_min = 0.9
            random_max = 1.1
            y_label = 'Texture Coefficient (TC)'
        
        if not normalize:
             ax.axhline(y=random_ref, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='Random')
        # ax.axhspan(random_min, random_max, alpha=0.15, color='gray', label='Random zone') # Removed per user request
        
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(f'{hkl} Peak', fontweight='bold', fontsize=13)
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
        if normalize:
            # Scale raw TC values to fraction for limit calculation
            all_vals = [v / 3.0 for v in all_tc]
            # Fraction limits: 0 to 1.0 (fixed)
            y_min = 0.0
            y_max = 1.0
        else:
            all_vals = all_tc
            # TC limits: auto-range with reasonable defaults
            y_min = min(0.5, min(all_vals) * 0.9)
            y_max = max(1.5, max(all_vals) * 1.1)
            
        for ax in axes[:len(hkl_list)]:
            ax.set_ylim(y_min, y_max)

    fig.suptitle('Texture Coefficient Evolution by Peak', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_texture_fraction_single(
    data: List[Dict[str, Any]],
    x_param: str = "time",
    output_path: Optional[str] = None,
    dpi: int = 1200,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (14, 11),
    metric: str = "fraction",  # "fraction" or "tc"
) -> plt.Figure:
    """Plot texture metrics with 2x2 grid (one panel per concentration).
    
    Args:
        metric: "fraction" (f=TC/N) or "tc" (raw TC values)
    """
    apply_axcsas_style()
    
    # Separate high and low quality data
    high_quality_data = [s for s in data if s.get('high_quality', True)]
    
    # Get all concentrations and hkl values
    concentrations = sorted(set(s.get('concentration', 0) for s in high_quality_data))
    all_hkls = sorted(set(
        hkl for s in high_quality_data for hkl in s.get('tc_values', {}).keys()
    ))
    n_peaks = len(all_hkls) if all_hkls else 3
    
    # Create 2x2 subplot grid
    # Adjust rows/cols based on number of concentrations
    n_conc = len(concentrations)
    n_cols = 2
    n_rows = (n_conc + 1) // 2
    
    # Resize figure if we have more rows
    if n_rows > 2:
        figsize = (figsize[0], figsize[1] * (n_rows/2))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()
    
    # Colors for different peaks
    peak_colors = get_color_palette(len(all_hkls))
    hkl_colors = dict(zip(all_hkls, peak_colors))
    
    # Line styles
    line_styles = ['-', '--', '-.', ':']
    markers = ['o', 's', '^', 'D', 'v']
    
    hkl_styles = {hkl: line_styles[i % len(line_styles)] for i, hkl in enumerate(all_hkls)}
    hkl_markers = {hkl: markers[i % len(markers)] for i, hkl in enumerate(all_hkls)}
    
    # Pre-calculate global max for TC metric
    global_max_tc = 1.6  # Default fallback
    if metric == "tc":
        all_tc_vals = [v for s in high_quality_data for v in s.get('tc_values', {}).values()]
        if all_tc_vals:
             global_max_tc = max(max(all_tc_vals) * 1.15, 1.6) # Ensure at least 1.6 if data is high

    # Plot each concentration in its own panel
    for idx, conc in enumerate(concentrations):
        ax = axes[idx]
        
        # Plot each peak for this concentration
        for hkl in all_hkls:
            x_values = []
            y_values = []
            
            # Collect data for this conc-hkl pair
            for sample in high_quality_data:
                if sample.get('concentration', 0) != conc:
                    continue
                x_val = sample.get(x_param, 0)
                tc_dict = sample.get('tc_values', {})
                if hkl in tc_dict:
                    x_values.append(x_val)
                    tc_val = tc_dict[hkl]
                    
                    if metric == "fraction":
                        y_values.append(tc_val / n_peaks)
                    else:
                        y_values.append(tc_val)
            
            if x_values:
                # Sort by x
                sorted_idx = np.argsort(x_values)
                x_sorted = np.array(x_values)[sorted_idx]
                y_sorted = np.array(y_values)[sorted_idx]
                
                color = hkl_colors.get(hkl, 'gray')
                linestyle = hkl_styles.get(hkl, '-')
                marker = hkl_markers.get(hkl, 'o')
                
                ax.plot(x_sorted, y_sorted, 
                       linestyle=linestyle,
                       marker=marker,
                       color=color,
                       linewidth=2,
                       markersize=8,
                       markeredgecolor='black',
                       markeredgewidth=0.5,
                       label=hkl,
                       alpha=0.85)
        
        # Labels and formatting
        x_label = 'Annealing Time (hours)' if x_param == 'time' else 'Concentration (mL/1.5L)'
        
        if metric == "fraction":
            y_label_text = 'Texture Fraction (f)'
            title_suffix = '(Fraction)'
            # Y limit for fraction (0-1)
            ax.set_ylim(0, 1.0)
            # Reference lines for random fraction
            # random_f = 1.0 / n_peaks
            # ax.axhline(random_f, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='Random')
        else:
            y_label_text = 'Texture Coefficient (TC)'
            title_suffix = '(TC Value)'
            # Y limit for TC (0-max)
            ax.set_ylim(0, global_max_tc)
            # Reference line for TC=1
            ax.axhline(1.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='Random (TC=1)')

        ax.set_xlabel(x_label, fontsize=12)
        ax.set_ylabel(y_label_text, fontsize=12)
        ax.set_title(f'{conc} mL/1.5L', fontsize=13, fontweight='bold')
        
        # Legend
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9, ncol=2)
        ax.grid(True, alpha=0.3)
        
    # Hide unused subplots
    for idx in range(len(concentrations), len(axes)):
        axes[idx].set_visible(False)
    
    title_metric = "Texture Fraction" if metric == "fraction" else "Texture Coefficient"
    fig.suptitle(f'{title_metric} Comparison by Concentration', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)
    
    if show:
        plt.show()
    
    return fig
    



