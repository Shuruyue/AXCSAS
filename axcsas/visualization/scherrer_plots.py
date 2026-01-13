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


def plot_scherrer_sizes(
    results: List[Dict[str, Any]],
    output_path: Optional[str] = None,
    dpi: int = 300,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (10, 6),
    show_validity: bool = True,
) -> plt.Figure:
    """Plot crystallite sizes from Scherrer analysis for each peak.
    繪製各峰 Scherrer 晶粒尺寸分析結果。
    
    Args:
        results: List of dictionaries with keys:
            - 'name': Sample name
            - 'peaks': List of dicts with 'hkl', 'size_nm', 'validity'
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        show_validity: Whether to show validity markers.
        
    Returns:
        Matplotlib Figure object.
        
    Example:
        >>> results = [
        ...     {'name': 'Sample1', 'peaks': [
        ...         {'hkl': '(111)', 'size_nm': 45.2, 'validity': 'VALID'},
        ...         {'hkl': '(200)', 'size_nm': 42.8, 'validity': 'VALID'},
        ...     ]},
        ... ]
        >>> fig = plot_scherrer_sizes(results)

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
        size_values = []
        validity_values = []

        for hkl in hkl_list:
            size = None
            validity = 'UNKNOWN'
            for peak in sample.get('peaks', []):
                if peak['hkl'] == hkl:
                    size = peak.get('size_nm', 0)
                    validity = peak.get('validity', 'VALID')
                    break
            size_values.append(size if size is not None else 0)
            validity_values.append(validity)

        offset = (sample_idx - n_samples / 2 + 0.5) * bar_width

        # Base bars
        bars = ax.bar(
            x_positions + offset,
            size_values,
            bar_width,
            label=sample.get('name', f'Sample {sample_idx + 1}'),
            color=colors[sample_idx % len(colors)],
            alpha=0.8,
            edgecolor='black',
            linewidth=0.5
        )

        # Add validity markers if enabled
        if show_validity:
            for bar, validity in zip(bars, validity_values):
                if validity == 'UNRELIABLE':
                    # Add cross-hatch for unreliable
                    bar.set_hatch('//')
                    bar.set_alpha(0.5)
                elif validity == 'WARNING':
                    # Add dot pattern for warning
                    bar.set_hatch('..')

    # Reference lines
    ax.axhline(y=2, color='red', linestyle=':', alpha=0.7,
              label='Precision limit (2 nm)')
    ax.axhline(y=200, color='orange', linestyle=':', alpha=0.7,
              label='Detection limit (200 nm)')

    ax.set_xlabel('Peak (hkl)')
    ax.set_ylabel('Crystallite Size (nm)')
    ax.set_title('Scherrer Crystallite Size Analysis')
    ax.set_xticks(x_positions)
    ax.set_xticklabels(hkl_list)
    ax.legend(loc='upper right', fontsize=9, ncol=2)

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_size_distribution(
    sizes: List[float],
    output_path: Optional[str] = None,
    dpi: int = 300,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (8, 6),
    bins: int = 15,
    sample_name: str = "Sample",
) -> plt.Figure:
    """Plot crystallite size distribution histogram.
    繪製晶粒尺寸分佈直方圖。
    
    Args:
        sizes: List of crystallite sizes in nm.
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        bins: Number of histogram bins.
        sample_name: Name of the sample for title.
        
    Returns:
        Matplotlib Figure object.

    """
    apply_axcsas_style()

    sizes = np.array([s for s in sizes if np.isfinite(s) and s > 0])

    if len(sizes) == 0:
        raise ValueError("No valid size data")

    fig, ax = plt.subplots(figsize=figsize)

    # Histogram
    n, bins_edges, patches = ax.hist(
        sizes, bins=bins,
        color='#0077BB',
        alpha=0.7,
        edgecolor='black',
        linewidth=0.8
    )

    # Statistics
    mean_size = np.mean(sizes)
    std_size = np.std(sizes)
    median_size = np.median(sizes)

    # Add mean and median lines
    ax.axvline(x=mean_size, color='red', linestyle='-', linewidth=2,
              label=f'Mean: {mean_size:.1f} nm')
    ax.axvline(x=median_size, color='green', linestyle='--', linewidth=2,
              label=f'Median: {median_size:.1f} nm')

    # Add standard deviation region
    ax.axvspan(mean_size - std_size, mean_size + std_size,
              alpha=0.2, color='red', label=f'±1σ: {std_size:.1f} nm')

    ax.set_xlabel('Crystallite Size (nm)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Crystallite Size Distribution - {sample_name}')
    ax.legend(loc='upper right')

    # Add text box with statistics
    textstr = f'N = {len(sizes)}\nMean = {mean_size:.1f} ± {std_size:.1f} nm'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=props)

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig
