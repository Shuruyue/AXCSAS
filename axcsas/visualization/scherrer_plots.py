"""Scherrer Size Visualization Module
==================================

Plots for crystallite size analysis results.
晶粒尺寸分析結果繪圖模組。
"""

from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .style import (
    apply_axcsas_style,
    get_color_palette,
    save_figure,
)


def plot_scherrer_evolution_by_peak(
    data: List[Dict[str, Any]],
    x_param: str = "time",
    output_path: Optional[str] = None,
    dpi: int = 1200,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (12, 8),
    instrument_limit: Optional[float] = None, # Not really applicable for size, but kept for interface consistency
) -> plt.Figure:
    """Plot Scherrer size evolution by peak (direction).
    
    Creates 3 subplots (111, 200, 220), each showing 4 concentration lines.
    X-axis: Annealing Time, Y-axis: Crystallite Size (nm).
    """
    apply_axcsas_style()

    # Collect all unique hkl values
    all_hkls = set()
    for sample in data:
        for peak in sample.get('peaks', []):
            all_hkls.add(peak['hkl'])

    hkl_list = sorted(list(all_hkls))
    n_hkls = len(hkl_list)

    if n_hkls == 0:
        raise ValueError("No peak data found in input")

    # Create subplots for each hkl (square-like aspect)
    n_cols = min(n_hkls, 3)
    n_rows = (n_hkls + n_cols - 1) // n_cols

    # Each subplot is 5x5 inches for square aspect
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5.5*n_rows), squeeze=False)
    axes = axes.flatten()

    x_labels = {
        'concentration': 'Leveler Concentration (mL/1.5L)',
        'time': 'Annealing Time (hours)',
    }
    x_label = x_labels.get(x_param, x_param)

    # Determine grouping variable (opposite of x_param)
    # Usually we plot vs Time, grouped by Concentration
    group_key = 'concentration'
    group_values = sorted(set(sample.get('concentration', 0) for sample in data))
    colors_by_group = get_color_palette(len(group_values))
    group_color_map = dict(zip(group_values, colors_by_group))
    label_format = lambda g: f'{g} mL/1.5L'
    
    # Different line styles and markers for different concentrations
    line_styles_list = ['-', '--', '-.', ':']
    markers_list = ['o', 's', '^', 'D']
    
    group_linestyle_map = {
        g: line_styles_list[i % len(line_styles_list)]
        for i, g in enumerate(group_values)
    }
    group_marker_map = {
        g: markers_list[i % len(markers_list)]
        for i, g in enumerate(group_values)
    }
    
    # Plot each hkl
    for idx, hkl in enumerate(hkl_list):
        ax = axes[idx]
        
        # Plot each group separately
        for group_val in group_values:
            x_values = []
            y_values = []
            
            for sample in data:
                if sample.get(group_key, 0) != group_val:
                    continue
                x_val = sample.get(x_param, 0)
                for peak in sample.get('peaks', []):
                    if peak['hkl'] == hkl:
                        val = peak.get('size_nm')
                        if val is not None and not np.isnan(val):
                            x_values.append(x_val)
                            y_values.append(val)
            
            if x_values:
                color = group_color_map[group_val]
                linestyle = group_linestyle_map[group_val]
                marker = group_marker_map[group_val]
                # Sort by x-value
                sorted_indices = np.argsort(x_values)
                x_sorted = np.array(x_values)[sorted_indices]
                y_sorted = np.array(y_values)[sorted_indices]
                
                # Retrieve errors (need to match x_values logic)
                y_errs = []
                for sample in data:
                    if sample.get(group_key, 0) == group_val:
                         for peak in sample.get('peaks', []):
                             if peak['hkl'] == hkl and peak.get('size_nm') is not None and not np.isnan(peak.get('size_nm')):
                                 y_errs.append(peak.get('size_err', 0.0))
                
                # Ensure length matches
                if len(y_errs) == len(x_values):
                    e_sorted = np.array(y_errs)[sorted_indices]
                    
                    label = label_format(group_val)
                    ax.errorbar(x_sorted, y_sorted, yerr=e_sorted, fmt=marker, c=color, markersize=7,
                                ecolor=color, elinewidth=1.5, capsize=4, alpha=0.8,
                                markeredgecolor='black', markeredgewidth=0.5, label=label)
                    ax.plot(x_sorted, y_sorted, c=color, alpha=0.7, linestyle=linestyle, linewidth=2.0)

        ax.set_xlabel(x_label)
        ax.set_ylabel('Crystallite Size (nm)')
        hkl_label = f"({hkl[0]}{hkl[1]}{hkl[2]})"
        ax.set_title(f'{hkl_label} Peak')
        
        # Instrument limit visualization
        # 1. Solid red for >300 nm (fully unreliable - beyond instrument resolution)
        ax.axhspan(300, 1000, color='red', alpha=0.15, label='>300 nm (Unreliable)')
        # 2. Hatched pattern for 110-300 nm (transition zone)
        ax.axhspan(110, 300, facecolor='none', edgecolor='red', hatch='///', alpha=0.3, label='110-300 nm (Caution)')
        # 3. Center line at 205 nm
        ax.axhline(y=205, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
        
        ax.legend(loc='upper right', fontsize=8, ncol=2)
        ax.set_box_aspect(1)  # Make subplot square
        ax.grid(True, alpha=0.3, linestyle=':')

    # Hide unused axes
    for idx in range(len(hkl_list), len(axes)):
        axes[idx].set_visible(False)

    # Set same Y-axis limits for all subplots
    all_y_values = []
    for sample in data:
        for peak in sample.get('peaks', []):
            val = peak.get('size_nm')
            if val is not None and not np.isnan(val):
                all_y_values.append(val)
    
    if all_y_values:
        y_min = 0
        y_max = 400  # Explicitly requested by user
        for ax in axes[:len(hkl_list)]:
            ax.set_ylim(y_min, y_max)

    fig.suptitle('Scherrer Size Evolution by Peak', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_scherrer_by_concentration(
    data: List[Dict[str, Any]],
    output_path: Optional[str] = None,
    dpi: int = 1200,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (11, 11), # Square figure to help subplots be square
    instrument_limit: Optional[float] = None,
) -> plt.Figure:
    """Plot Scherrer size evolution with subplots by concentration.
    
    Creates 4 subplots (one per concentration), each showing 3 peak lines (111, 200, 220).
    X-axis: Annealing Time, Y-axis: Crystallite Size (nm).
    Subplots are square.
    """
    apply_axcsas_style()
    
    # Get unique concentrations and peaks
    concentrations = sorted(set(sample.get('concentration', 0) for sample in data))
    all_hkls = set()
    for sample in data:
        for peak in sample.get('peaks', []):
            all_hkls.add(peak['hkl'])
    hkl_list = sorted(list(all_hkls))
    
    if len(concentrations) == 0:
        raise ValueError("No concentration data found")
    
    # Create subplots: one per concentration (2x2 layout for 4 concentrations)
    n_conc = len(concentrations)
    n_cols = 2
    n_rows = (n_conc + 1) // 2
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()
    
    # Colors for each peak direction
    peak_colors = get_color_palette(len(hkl_list))
    hkl_color_map = dict(zip(hkl_list, peak_colors))
    
    # Collect all size values for consistent Y-axis
    all_sizes = []
    for sample in data:
        for peak in sample.get('peaks', []):
            val = peak.get('size_nm')
            if val is not None and not np.isnan(val):
                all_sizes.append(val)
    
    y_min = 0
    y_max = max(all_sizes) * 1.1 if all_sizes else 100.0
    
    # Plot each concentration
    for idx, conc in enumerate(concentrations):
        ax = axes[idx]
        
        # Plot each peak direction
        for hkl in hkl_list:
            x_values = []
            y_values = []
            
            for sample in data:
                if sample.get('concentration', 0) != conc:
                    continue
                time_val = sample.get('time', 0)
                for peak in sample.get('peaks', []):
                    if peak['hkl'] == hkl:
                        val = peak.get('size_nm')
                        if val is not None and not np.isnan(val):
                            x_values.append(time_val)
                            y_values.append(val)
            
                # Define styles using standard Miller Index notation (111)
                # Matches style in compare_ka_methods.py
                
                # Format: (1, 1, 1) -> "(111)"
                hkl_label = f"({hkl[0]}{hkl[1]}{hkl[2]})"
                
                color = hkl_color_map[hkl]
                
                marker_map = {
                    '(111)': 'o', '(200)': 's', '(220)': '^', 
                    '(311)': 'D', '(222)': 'v'
                }
                linestyle_map = {
                    '(111)': '-', '(200)': '--', '(220)': ':', 
                    '(311)': '-.', '(222)': '-'
                }
                
                marker = marker_map.get(hkl_label, 'o')
                linestyle = linestyle_map.get(hkl_label, '-')

                # Sort by time
                sorted_indices = np.argsort(x_values)
                x_sorted = np.array(x_values)[sorted_indices]
                y_sorted = np.array(y_values)[sorted_indices]
                
                # Retrieve errors
                y_errs = []
                for sample in data:
                    if sample.get('concentration', 0) == conc:
                         for peak in sample.get('peaks', []):
                             if peak['hkl'] == hkl and peak.get('size_nm') is not None and not np.isnan(peak.get('size_nm')):
                                 y_errs.append(peak.get('size_err', 0.0))
                
                if len(y_errs) == len(x_values):
                    e_sorted = np.array(y_errs)[sorted_indices]
                    
                    ax.errorbar(x_sorted, y_sorted, yerr=e_sorted, fmt=marker, c=color, markersize=7,
                                ecolor=color, elinewidth=1.5, capsize=4, alpha=0.8,
                                markeredgecolor='black', markeredgewidth=0.5, label=hkl_label)
                    ax.plot(x_sorted, y_sorted, c=color, alpha=0.7, linestyle=linestyle, linewidth=2.0)
        
        ax.set_xlabel('Annealing Time (hours)')
        ax.set_ylabel('Crystallite Size (nm)')
        ax.set_title(f'{conc} mL/1.5L', fontsize=12, fontweight='bold')
        
        # Instrument limit visualization
        # 1. Solid red for >300 nm (fully unreliable - beyond instrument resolution)
        ax.axhspan(300, 1000, color='red', alpha=0.15, label='>300 nm (Unreliable)')
        # 2. Hatched pattern for 110-300 nm (transition zone)
        ax.axhspan(110, 300, facecolor='none', edgecolor='red', hatch='///', alpha=0.3, label='110-300 nm (Caution)')
        # 3. Center line at 205 nm
        ax.axhline(y=205, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
        
        ax.legend(loc='upper right', fontsize=8)
        ax.set_ylim(0, 400)  # Explicitly requested by user
        ax.grid(True, alpha=0.3, linestyle=':')
        ax.set_box_aspect(1)  # Make subplot square
    
    # Hide unused axes
    for idx in range(len(concentrations), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle('Scherrer Size Evolution by Concentration', fontsize=16, fontweight='bold', y=1.02)
    # Add Physical Limitation Footnote
    LIMITATION_NOTE = (
        "Note: Size > 200 nm approaches instrumental resolution limit (FWHM ≈ 0.04°).\n"
        "註：晶粒尺寸 > 200 nm 接近儀器解析極限，數值僅供趨勢參考。"
    )
    fig.text(0.5, 0.01, LIMITATION_NOTE, ha='center', va='bottom', 
             fontsize=10, style='italic', color='#555555', 
             bbox=dict(facecolor='#f0f0f0', alpha=0.5, edgecolor='none', pad=5))

    plt.tight_layout(rect=[0, 0.05, 1, 0.96]) # Adjust for footer
    
    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)
    
    if show:
        plt.show()
    
    return fig
