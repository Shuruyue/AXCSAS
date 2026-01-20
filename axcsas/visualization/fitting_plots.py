"""Peak Fitting Visualization Module
=================================

Diagnostic plots for peak fitting analysis.
峰型擬合診斷繪圖模組。
"""

from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .style import (
    COLORBLIND_SAFE,
    apply_axcsas_style,
    save_figure,
)


def plot_peak_fit(
    two_theta: np.ndarray,
    intensity_obs: np.ndarray,
    intensity_fit: np.ndarray,
    peak_params: Optional[Dict[str, Any]] = None,
    output_path: Optional[str] = None,
    dpi: int = 2400,  # Ultra-high quality (期刊出版標準)
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
    peak_name: str = "",
) -> plt.Figure:
    """Plot peak fitting comparison (observed vs fitted).
    繪製峰型擬合對比圖（觀測值 vs 擬合值）。
    
    Args:
        two_theta: 2θ array (degrees).
        intensity_obs: Observed intensity array.
        intensity_fit: Fitted intensity array.
        peak_params: Optional dict with 'center', 'fwhm', 'eta', 'r_squared'.
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        peak_name: Name/hkl of the peak for title.
        
    Returns:
        Matplotlib Figure object.

    """
    apply_axcsas_style()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize,
                                   gridspec_kw={'height_ratios': [3, 1]},
                                   sharex=True)

    # Main plot: Observed vs Fitted
    ax1.scatter(two_theta, intensity_obs, s=20, c='black', alpha=0.5,
               label='Observed', zorder=3)
    ax1.plot(two_theta, intensity_fit, color=COLORBLIND_SAFE[0],
            linewidth=2.5, label='Fitted (Pseudo-Voigt)', zorder=4)

    # Fill area under fit
    ax1.fill_between(two_theta, 0, intensity_fit, alpha=0.2,
                    color=COLORBLIND_SAFE[0])

    ax1.set_ylabel('Intensity (counts)')
    title = 'Peak Fitting Analysis'
    if peak_name:
        title += f' - {peak_name}'
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right')

    # Add parameters text box
    if peak_params:
        param_text = '\n'.join([
            f"Center: {peak_params.get('center', 0):.3f}°",
            f"FWHM: {peak_params.get('fwhm', 0):.4f}°",
            f"η: {peak_params.get('eta', 0):.3f}",
            f"R²: {peak_params.get('r_squared', 0):.4f}",
        ])
        props = dict(boxstyle='round', facecolor='lightyellow',
                    alpha=0.9, edgecolor='gray')
        ax1.text(0.02, 0.98, param_text, transform=ax1.transAxes,
                fontsize=10, verticalalignment='top',
                family='monospace', bbox=props)

    # Residual plot
    residuals = intensity_obs - intensity_fit
    ax2.scatter(two_theta, residuals, s=15, c=COLORBLIND_SAFE[1], alpha=0.7)
    ax2.axhline(y=0, color='red', linestyle='-', linewidth=1.5)
    ax2.fill_between(two_theta, residuals, 0, alpha=0.3, color=COLORBLIND_SAFE[1])

    ax2.set_xlabel('2θ (degrees)')
    ax2.set_ylabel('Residual')
    ax2.set_title('Residuals (Observed - Fitted)', fontsize=11)

    # Residual statistics
    rmse = np.sqrt(np.mean(residuals**2))
    ax2.text(0.98, 0.95, f'RMSE = {rmse:.1f}', transform=ax2.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_doublet_comparison(
    two_theta: np.ndarray,
    intensity_obs: np.ndarray,
    intensity_single: np.ndarray,
    intensity_doublet: np.ndarray,
    ka1_center: float,
    ka2_center: float,
    output_path: Optional[str] = None,
    dpi: int = 2400,  # Ultra-high quality (期刊出版標準)
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (12, 8),
    peak_name: str = "",
) -> plt.Figure:
    """Compare single peak vs Kα doublet fitting.
    比較單峰與 Kα 雙峰擬合。
    
    Args:
        two_theta: 2θ array (degrees).
        intensity_obs: Observed intensity.
        intensity_single: Single peak fitted intensity.
        intensity_doublet: Doublet fitted intensity.
        ka1_center: Kα₁ peak center.
        ka2_center: Kα₂ peak center.
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        peak_name: Name of the peak for title.
        
    Returns:
        Matplotlib Figure object.

    """
    apply_axcsas_style()

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Top left: Observed data with both fits
    ax1 = axes[0, 0]
    ax1.scatter(two_theta, intensity_obs, s=10, c='black', alpha=0.4,
               label='Observed')
    ax1.plot(two_theta, intensity_single, color=COLORBLIND_SAFE[1],
            linewidth=2, linestyle='--', label='Single peak fit')
    ax1.plot(two_theta, intensity_doublet, color=COLORBLIND_SAFE[0],
            linewidth=2, label='Doublet fit')

    # Mark Kα1 and Kα2 positions
    ax1.axvline(x=ka1_center, color='green', linestyle=':', alpha=0.7)
    ax1.axvline(x=ka2_center, color='purple', linestyle=':', alpha=0.7)
    ax1.text(ka1_center, ax1.get_ylim()[1] * 0.95, 'Kα₁', fontsize=10,
            ha='center', color='green')
    ax1.text(ka2_center, ax1.get_ylim()[1] * 0.9, 'Kα₂', fontsize=10,
            ha='center', color='purple')

    ax1.set_xlabel('2θ (degrees)')
    ax1.set_ylabel('Intensity')
    ax1.set_title(f'Fitting Comparison - {peak_name}')
    ax1.legend(loc='upper right', fontsize=9)

    # Top right: Residual comparison
    ax2 = axes[0, 1]
    res_single = intensity_obs - intensity_single
    res_doublet = intensity_obs - intensity_doublet

    ax2.plot(two_theta, res_single, color=COLORBLIND_SAFE[1],
            alpha=0.7, label=f'Single (RMSE={np.sqrt(np.mean(res_single**2)):.1f})')
    ax2.plot(two_theta, res_doublet, color=COLORBLIND_SAFE[0],
            alpha=0.7, label=f'Doublet (RMSE={np.sqrt(np.mean(res_doublet**2)):.1f})')
    ax2.axhline(y=0, color='red', linestyle='-', linewidth=1)

    ax2.set_xlabel('2θ (degrees)')
    ax2.set_ylabel('Residual')
    ax2.set_title('Residual Comparison')
    ax2.legend(loc='upper right', fontsize=9)

    # Bottom left: R² comparison bar chart
    ax3 = axes[1, 0]
    ss_tot = np.sum((intensity_obs - np.mean(intensity_obs))**2)
    r2_single = 1 - np.sum(res_single**2) / ss_tot if ss_tot > 0 else 0
    r2_doublet = 1 - np.sum(res_doublet**2) / ss_tot if ss_tot > 0 else 0

    methods = ['Single Peak', 'Kα Doublet']
    r2_values = [r2_single, r2_doublet]
    colors = [COLORBLIND_SAFE[1], COLORBLIND_SAFE[0]]

    bars = ax3.bar(methods, r2_values, color=colors, alpha=0.8,
                   edgecolor='black', linewidth=1)

    # Add value labels on bars
    for bar, val in zip(bars, r2_values):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{val:.4f}', ha='center', fontsize=11, fontweight='bold')

    ax3.set_ylabel('R²')
    ax3.set_title('Goodness of Fit Comparison')
    ax3.set_ylim(min(r2_values) * 0.95, 1.01)
    ax3.axhline(y=0.95, color='green', linestyle='--', alpha=0.7,
               label='Target (0.95)')
    ax3.legend(loc='lower right')

    # Bottom right: Summary table
    ax4 = axes[1, 1]
    ax4.axis('off')

    table_data = [
        ['Method', 'R²', 'RMSE', 'Improvement'],
        ['Single Peak', f'{r2_single:.4f}', f'{np.sqrt(np.mean(res_single**2)):.1f}', '-'],
        ['Kα Doublet', f'{r2_doublet:.4f}', f'{np.sqrt(np.mean(res_doublet**2)):.1f}',
         f'{(r2_doublet - r2_single) / r2_single * 100:+.2f}%' if r2_single > 0 else 'N/A'],
    ]

    table = ax4.table(
        cellText=table_data,
        cellLoc='center',
        loc='center',
        colWidths=[0.3, 0.2, 0.2, 0.3]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.8)

    # Style header row
    for i in range(4):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    # Highlight better method
    better_row = 2 if r2_doublet > r2_single else 1
    for i in range(4):
        table[(better_row, i)].set_facecolor('#C6EFCE')

    ax4.set_title('Summary', fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_fit_residuals(
    two_theta: np.ndarray,
    residuals: np.ndarray,
    output_path: Optional[str] = None,
    dpi: int = 2400,  # Ultra-high quality (期刊出版標準)
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (10, 5),
) -> plt.Figure:
    """Plot peak fitting residuals with statistics.
    繪製峰型擬合殘差統計圖。
    
    Args:
        two_theta: 2θ array (degrees).
        residuals: Residual array (observed - fitted).
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        
    Returns:
        Matplotlib Figure object.

    """
    apply_axcsas_style()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Left: Residuals vs 2θ
    ax1.scatter(two_theta, residuals, s=20, c=COLORBLIND_SAFE[0], alpha=0.6)
    ax1.axhline(y=0, color='red', linestyle='-', linewidth=1.5)

    # Rolling mean to show systematic trends
    if len(residuals) > 10:
        window = max(5, len(residuals) // 10)
        rolling_mean = np.convolve(residuals, np.ones(window)/window, mode='valid')
        x_rolling = two_theta[window//2:window//2+len(rolling_mean)]
        ax1.plot(x_rolling, rolling_mean, color='orange', linewidth=2,
                label='Rolling mean')

    ax1.set_xlabel('2θ (degrees)')
    ax1.set_ylabel('Residual')
    ax1.set_title('Residuals vs 2θ')
    ax1.legend(loc='upper right')

    # Right: Residual histogram
    ax2.hist(residuals, bins=20, color=COLORBLIND_SAFE[0], alpha=0.7,
            edgecolor='black', density=True)

    # Overlay normal distribution
    from scipy import stats
    mu, sigma = np.mean(residuals), np.std(residuals)
    x_norm = np.linspace(mu - 4*sigma, mu + 4*sigma, 100)
    ax2.plot(x_norm, stats.norm.pdf(x_norm, mu, sigma),
            color='red', linewidth=2, label=f'Normal (μ={mu:.1f}, σ={sigma:.1f})')

    ax2.axvline(x=0, color='green', linestyle='--', linewidth=1.5)
    ax2.set_xlabel('Residual')
    ax2.set_ylabel('Density')
    ax2.set_title('Residual Distribution')
    ax2.legend(loc='upper right')

    # Add statistics
    rmse = np.sqrt(np.mean(residuals**2))
    skewness = stats.skew(residuals)

    fig.text(0.5, 0.02,
            f'RMSE = {rmse:.2f} | Mean = {mu:.2f} | Skewness = {skewness:.3f}',
            ha='center', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout(rect=[0, 0.05, 1, 1])

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig
