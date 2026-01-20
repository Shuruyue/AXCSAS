"""FWHM Visualization Module
=========================

Plots for FWHM (Full Width at Half Maximum) analysis.
FWHM（半高寬）分析繪圖模組。
"""

from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .style import (
    apply_axcsas_style,
    get_color_palette,
    get_peak_color,
    save_figure,
)


def plot_fwhm_evolution(
    data: List[Dict[str, Any]],
    x_param: str = "concentration",
    output_path: Optional[str] = None,
    dpi: int = 2400,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (12, 8),
    instrument_limit: Optional[float] = 0.05,
) -> plt.Figure:
    """Plot FWHM evolution across samples.
    繪製 FWHM 隨樣品參數演化圖。
    
    Args:
        data: List of dictionaries with keys:
            - 'name': Sample name
            - 'concentration' or 'time': X-axis value
            - 'peaks': List of dicts with 'hkl' (str) and 'fwhm' (float)
        x_param: X-axis parameter ("concentration" or "time").
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format (png, svg, pdf).
        show: Whether to display figure.
        figsize: Figure size (width, height).
        instrument_limit: Instrument broadening limit line (degrees).
        
    Returns:
        Matplotlib Figure object.
        
    Example:
        >>> data = [
        ...     {'name': 'Sample1', 'concentration': 0, 
        ...      'peaks': [{'hkl': '(111)', 'fwhm': 0.25}]},
        ...     {'name': 'Sample2', 'concentration': 9, 
        ...      'peaks': [{'hkl': '(111)', 'fwhm': 0.28}]},
        ... ]
        >>> fig = plot_fwhm_evolution(data, x_param='concentration')

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

    # X-axis label
    x_labels = {
        'concentration': 'Leveler Concentration (mL/1.5L)',
        'time': 'Annealing Time (hours)',
        'age': 'Sample Age (hours)',
    }
    x_label = x_labels.get(x_param, x_param)

    # Determine grouping variable (opposite of x_param)
    if x_param == 'concentration':
        # Group by time when plotting against concentration
        group_key = 'time'
        group_values = sorted(set(sample.get('time', 0) for sample in data))
        colors_by_group = get_color_palette(len(group_values))
        group_color_map = dict(zip(group_values, colors_by_group))
        label_format = lambda g: f'{g:.0f}h'
    else:
        # Group by concentration when plotting against time
        group_key = 'concentration'
        group_values = sorted(set(sample.get('concentration', 0) for sample in data))
        colors_by_group = get_color_palette(len(group_values))
        group_color_map = dict(zip(group_values, colors_by_group))
        label_format = lambda g: f'{g} mL/1.5L'
    
    # Plot each hkl
    for idx, hkl in enumerate(hkl_list):
        ax = axes[idx]
        
        # Plot each group separately
        for group_val in group_values:
            x_values = []
            y_values = []
            
            # Re-collect x, y, y_err
            x_vals = []
            y_vals = []
            y_errs = []
            
            for sample in data:
                if sample.get(group_key, 0) != group_val:
                    continue
                x_val = sample.get(x_param, sample.get('time', 0))
                for peak in sample.get('peaks', []):
                    if peak['hkl'] == hkl:
                        x_vals.append(x_val)
                        y_vals.append(peak['fwhm'])
                        y_errs.append(peak.get('fwhm_error', 0.0))
                            
            if x_vals:
                color = group_color_map[group_val]
                # Sort
                indices = np.argsort(x_vals)
                x_s = np.array(x_vals)[indices]
                y_s = np.array(y_vals)[indices]
                e_s = np.array(y_errs)[indices]
                
                label = label_format(group_val)
                
                # Plot with error bars
                ax.errorbar(x_s, y_s, yerr=e_s, fmt='o', c=color, markersize=8, 
                           ecolor=color, elinewidth=1.5, capsize=4, alpha=0.8,
                           markeredgecolor='black', markeredgewidth=0.5, label=label)
                ax.plot(x_s, y_s, c=color, alpha=0.6, linestyle='-', linewidth=2.0)

        # Instrument limit visualization
        # 1. Solid red for <0.03° (fully unreliable - below instrument resolution)
        ax.axhspan(0, 0.03, color='red', alpha=0.15, label='<0.03° (Unreliable)')
        # 2. Hatched pattern for 0.03-0.08° (transition zone)
        ax.axhspan(0.03, 0.08, facecolor='none', edgecolor='red', hatch='///', alpha=0.3, label='0.03-0.08° (Caution)')
        # 3. Center line at 0.055°
        ax.axhline(y=0.055, color='red', linestyle='--', alpha=0.5, linewidth=1.5)



        ax.set_xlabel(x_label)
        ax.set_ylabel('FWHM (°)')
        ax.set_title(f'{hkl} Peak')
        ax.legend(loc='upper right', fontsize=8, ncol=2)
        ax.set_box_aspect(1)  # Make subplot square

    # Hide unused axes
    for idx in range(len(hkl_list), len(axes)):
        axes[idx].set_visible(False)

    # Set same Y-axis limits for all subplots
    all_y_values = []
    for sample in data:
        for peak in sample.get('peaks', []):
            if peak['fwhm'] is not None and not np.isnan(peak['fwhm']):
                all_y_values.append(peak['fwhm'])
    
    if all_y_values:
        y_min = 0
        y_max = max(all_y_values) * 1.1
        for ax in axes[:len(hkl_list)]:
            ax.set_ylim(y_min, y_max)

    fig.suptitle('FWHM Evolution by Peak', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_fwhm_by_peak(
    results: List[Dict[str, Any]],
    output_path: Optional[str] = None,
    dpi: int = 2400,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (10, 6),
    y_limit: Optional[Tuple[float, float]] = None,
    instrument_limit: Optional[float] = 0.05,
) -> plt.Figure:
    """Plot FWHM comparison across different peaks for all samples.
    繪製所有樣品各峰 FWHM 比較圖。
    
    Args:
        results: List of dictionaries with keys:
            - 'name': Sample name
            - 'peaks': List of dicts with 'hkl' (str) and 'fwhm' (float)
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        y_limit: Optional Y-axis limits (min, max).
        instrument_limit: Instrument broadening limit.
        
    Returns:
        Matplotlib Figure object.

    """
    apply_axcsas_style()

    # Collect all unique hkl values
    all_hkls = set()
    for sample in results:
        for peak in sample.get('peaks', []):
            all_hkls.add(peak['hkl'])

    hkl_list = sorted(list(all_hkls))
    n_samples = len(results)
    n_peaks = len(hkl_list)

    if n_peaks == 0:
        raise ValueError("No peak data found in input")

    fig, ax = plt.subplots(figsize=figsize)

    # Bar width and positions
    bar_width = 0.8 / n_samples
    x_positions = np.arange(n_peaks)
    colors = get_color_palette(n_samples)

    for sample_idx, sample in enumerate(results):
        fwhm_values = []
        for hkl in hkl_list:
            fwhm = None
            for peak in sample.get('peaks', []):
                if peak['hkl'] == hkl:
                    fwhm = peak['fwhm']
                    break
            fwhm_values.append(fwhm if fwhm is not None else 0)

        offset = (sample_idx - n_samples / 2 + 0.5) * bar_width
        bars = ax.bar(
            x_positions + offset,
            fwhm_values,
            bar_width,
            label=sample.get('name', f'Sample {sample_idx + 1}'),
            color=colors[sample_idx % len(colors)],
            alpha=0.8,
            edgecolor='black',
            linewidth=0.5
        )

    # Instrument limit range (Typical: 0.03-0.08)
    ax.axhspan(0.03, 0.08, color='red', alpha=0.1, label='Instrument limit (0.03°-0.08°)')
    if instrument_limit is not None:
        ax.axhline(y=instrument_limit, color='red', linestyle='--', alpha=0.3, linewidth=1)

    ax.set_xlabel('Peak (hkl)')
    ax.set_ylabel('FWHM (°)')
    ax.set_title('FWHM Comparison by Peak')
    ax.set_xticks(x_positions)
    ax.set_xticklabels(hkl_list)
    ax.legend(loc='upper right', fontsize=9, ncol=2)

    if y_limit:
        ax.set_ylim(y_limit)

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_fwhm_by_concentration(
    data: List[Dict[str, Any]],
    output_path: Optional[str] = None,
    dpi: int = 2400,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (16, 10),
    instrument_limit: Optional[float] = 0.05,
) -> plt.Figure:
    """Plot FWHM evolution with subplots by concentration.
    
    Creates 4 subplots (one per concentration), each showing 3 peak lines (111, 200, 220).
    X-axis: Annealing Time, Y-axis: FWHM.
    
    Args:
        data: List of sample dictionaries with 'time', 'concentration', 'peaks'.
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        instrument_limit: Instrument limit line.
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
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11, 11), squeeze=False)
    axes = axes.flatten()
    
    # Colors for each peak direction
    peak_colors = get_color_palette(len(hkl_list))
    hkl_color_map = dict(zip(hkl_list, peak_colors))
    
    # Collect all FWHM values for consistent Y-axis
    all_fwhm = []
    for sample in data:
        for peak in sample.get('peaks', []):
            if peak['fwhm'] is not None and not np.isnan(peak['fwhm']):
                all_fwhm.append(peak['fwhm'])
    
    y_min = 0
    y_max = max(all_fwhm) * 1.1 if all_fwhm else 1.0
    
    # Plot each concentration
    for idx, conc in enumerate(concentrations):
        ax = axes[idx]
        
        # Plot each peak direction
        for hkl in hkl_list:
            x_values = []
            y_values = []
            
            x_vals = []
            y_vals = []
            y_errs = []
            
            for sample in data:
                if sample.get('concentration', 0) != conc:
                    continue
                time_val = sample.get('time', 0)
                for peak in sample.get('peaks', []):
                    if peak['hkl'] == hkl:
                        x_vals.append(time_val)
                        y_vals.append(peak['fwhm'])
                        y_errs.append(peak.get('fwhm_error', 0.0))
            
            if x_vals:
                color = hkl_color_map[hkl]
                # Sort by time
                indices = np.argsort(x_vals)
                x_s = np.array(x_vals)[indices]
                y_s = np.array(y_vals)[indices]
                e_s = np.array(y_errs)[indices]
                
                ax.errorbar(x_s, y_s, yerr=e_s, fmt='o', c=color, markersize=8,
                           ecolor=color, elinewidth=1.5, capsize=4, alpha=0.8,
                           markeredgecolor='black', markeredgewidth=0.5, label=hkl)
                ax.plot(x_s, y_s, c=color, alpha=0.6, linestyle='-', linewidth=2.0)
        
        # Instrument limit visualization
        # 1. Solid red for <0.03° (fully unreliable - below instrument resolution)
        ax.axhspan(0, 0.03, color='red', alpha=0.15, label='<0.03° (Unreliable)')
        # 2. Hatched pattern for 0.03-0.08° (transition zone)
        ax.axhspan(0.03, 0.08, facecolor='none', edgecolor='red', hatch='///', alpha=0.3, label='0.03-0.08° (Caution)')
        # 3. Center line at 0.055°
        ax.axhline(y=0.055, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
        
        ax.set_xlabel('Annealing Time (hours)')
        ax.set_ylabel('FWHM (°)')
        ax.set_title(f'{conc} mL/1.5L', fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', fontsize=8)
        ax.set_ylim(y_min, y_max)
        ax.grid(True, alpha=0.3)
        ax.set_box_aspect(1)  # Make subplot square
    
    # Hide unused axes
    for idx in range(len(concentrations), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle('FWHM Evolution by Concentration', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)
    
    if show:
        plt.show()
    
    return fig
