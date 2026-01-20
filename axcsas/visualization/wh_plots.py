"""Williamson-Hall Visualization Module
====================================

Plots for Williamson-Hall size/strain analysis.
Williamson-Hall 尺寸/應變分析繪圖模組。
"""

from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .style import (
    COLORBLIND_SAFE,
    apply_axcsas_style,
    save_figure,
)
from axcsas.core.constants import CU_KA1, SCHERRER_K


def plot_williamson_hall(
    two_theta: np.ndarray,
    fwhm_sample: np.ndarray,
    fit_result: Optional[Dict[str, Any]] = None,
    hkl_labels: Optional[List[str]] = None,
    output_path: Optional[str] = None,
    dpi: int = 2400,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (10, 7),
    sample_name: str = "Sample",
) -> plt.Figure:
    """Plot Williamson-Hall linear fit diagram.
    繪製 Williamson-Hall 線性擬合圖。
    
    W-H Equation: β cos θ = (Kλ/D) + 4ε sin θ
    
    Args:
        two_theta: Array of 2θ peak positions (degrees).
        fwhm_sample: Array of sample FWHM values (degrees).
        fit_result: Optional dict with 'slope', 'intercept', 'r_squared',
                   'size_nm', 'microstrain'.
        hkl_labels: Optional list of hkl labels for each point.
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        sample_name: Name for title.
        
    Returns:
        Matplotlib Figure object.
        
    Example:
        >>> two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        >>> fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        >>> fig = plot_williamson_hall(two_theta, fwhm)

    """
    apply_axcsas_style()

    # Calculate W-H coordinates
    theta_rad = np.radians(np.asarray(two_theta) / 2.0)
    beta_rad = np.radians(np.asarray(fwhm_sample))

    x_data = np.sin(theta_rad)  # sin(θ)
    y_data = beta_rad * np.cos(theta_rad)  # β cos(θ)

    fig, ax = plt.subplots(figsize=figsize)

    # Scatter plot with data points
    scatter = ax.scatter(
        x_data, y_data,
        s=100, c=COLORBLIND_SAFE[0],
        edgecolors='black', linewidths=1,
        zorder=5, label='Experimental data'
    )

    # Add hkl labels if provided
    if hkl_labels:
        for i, label in enumerate(hkl_labels):
            ax.annotate(
                label, (x_data[i], y_data[i]),
                xytext=(5, 5), textcoords='offset points',
                fontsize=10, fontweight='bold'
            )

    # Perform or use provided linear fit
    if fit_result is None:
        # Perform linear regression
        from scipy.stats import linregress
        reg = linregress(x_data, y_data)
        slope = reg.slope
        intercept = reg.intercept
        r_squared = reg.rvalue ** 2
    else:
        slope = fit_result.get('slope', 0)
        intercept = fit_result.get('intercept', 0)
        r_squared = fit_result.get('r_squared', 0)

    # Plot regression line
    x_fit = np.linspace(x_data.min() * 0.9, x_data.max() * 1.1, 100)
    y_fit = slope * x_fit + intercept

    ax.plot(x_fit, y_fit, color=COLORBLIND_SAFE[1], linewidth=2,
           linestyle='--', label='Linear fit', zorder=3)

    # Fill confidence region (optional visualization)
    ax.fill_between(x_fit, y_fit * 0.95, y_fit * 1.05,
                   alpha=0.1, color=COLORBLIND_SAFE[1])

    # Calculate physical quantities
    # Calculate physical quantities
    wavelength = CU_KA1  # Cu Kα1
    K = SCHERRER_K.default

    if intercept > 0:
        size_angstrom = K * wavelength / intercept
        size_nm = size_angstrom / 10.0
    else:
        size_nm = float('inf')

    microstrain = slope / 4.0

    # Axis labels
    ax.set_xlabel(r'sin($\theta$)', fontsize=14)
    ax.set_ylabel(r'$\beta$ cos($\theta$) (rad)', fontsize=14)
    ax.set_title(f'Williamson-Hall Analysis - {sample_name}', fontsize=16)

    # Legend
    ax.legend(loc='upper left', fontsize=11)

    # Results text box
    textstr = '\n'.join([
        'Results:',
        f'  Slope = {slope:.5f}',
        f'  Intercept = {intercept:.5f}',
        f'  R² = {r_squared:.4f}',
        '',
        f'  D = {size_nm:.1f} nm',
        f'  ε = {microstrain:.2e}',
    ])
    props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='gray')
    ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='bottom', horizontalalignment='right',
           bbox=props, family='monospace')

    # Quality indicator
    quality_threshold = 0.95
    if r_squared > quality_threshold:
        quality = "Excellent"
        box_facecolor = 'lightgreen'
    elif r_squared > 0.85:
        quality = 'ACCEPTABLE'
        quality_color = 'orange'
    else:
        quality = 'POOR (Anisotropy?)'
        quality_color = 'red'

    ax.text(0.02, 0.98, f'Quality: {quality}', transform=ax.transAxes,
           fontsize=11, fontweight='bold', color=quality_color,
           verticalalignment='top')

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig


def plot_wh_residuals(
    two_theta: np.ndarray,
    fwhm_sample: np.ndarray,
    fit_result: Optional[Dict[str, Any]] = None,
    hkl_labels: Optional[List[str]] = None,
    output_path: Optional[str] = None,
    dpi: int = 2400,
    format: str = "png",
    show: bool = True,
    figsize: Tuple[float, float] = (10, 5),
) -> plt.Figure:
    """Plot Williamson-Hall fit residuals.
    繪製 Williamson-Hall 擬合殘差圖。
    
    Args:
        two_theta: Array of 2θ peak positions (degrees).
        fwhm_sample: Array of sample FWHM values (degrees).
        fit_result: Optional dict with 'slope', 'intercept'.
        hkl_labels: Optional list of hkl labels.
        output_path: Optional path to save figure.
        dpi: Output resolution.
        format: Output format.
        show: Whether to display figure.
        figsize: Figure size.
        
    Returns:
        Matplotlib Figure object.

    """
    apply_axcsas_style()

    # Calculate W-H coordinates
    theta_rad = np.radians(np.asarray(two_theta) / 2.0)
    beta_rad = np.radians(np.asarray(fwhm_sample))

    x_data = np.sin(theta_rad)
    y_data = beta_rad * np.cos(theta_rad)

    # Get or calculate fit
    if fit_result is None:
        from scipy.stats import linregress
        reg = linregress(x_data, y_data)
        slope = reg.slope
        intercept = reg.intercept
    else:
        slope = fit_result.get('slope', 0)
        intercept = fit_result.get('intercept', 0)

    # Calculate residuals
    y_pred = slope * x_data + intercept
    residuals = y_data - y_pred

    fig, ax = plt.subplots(figsize=figsize)

    # Stem plot for residuals
    markerline, stemlines, baseline = ax.stem(
        x_data, residuals * 1000,  # Convert to mrad for visibility
        basefmt='k-', linefmt='gray', markerfmt='o'
    )
    plt.setp(markerline, markersize=10, color=COLORBLIND_SAFE[0])

    # Zero line
    ax.axhline(y=0, color='red', linestyle='-', linewidth=1.5)

    # Add hkl labels if provided
    if hkl_labels:
        for i, label in enumerate(hkl_labels):
            offset = 5 if residuals[i] >= 0 else -15
            ax.annotate(
                label, (x_data[i], residuals[i] * 1000),
                xytext=(0, offset), textcoords='offset points',
                fontsize=10, fontweight='bold', ha='center'
            )

    ax.set_xlabel(r'sin($\theta$)', fontsize=14)
    ax.set_ylabel('Residual (mrad)', fontsize=14)
    ax.set_title('Williamson-Hall Fit Residuals')

    # Statistics
    rmse = np.sqrt(np.mean(residuals**2)) * 1000
    ax.text(0.98, 0.98, f'RMSE = {rmse:.3f} mrad', transform=ax.transAxes,
           fontsize=11, verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path, dpi=dpi, format=format)

    if show:
        plt.show()

    return fig
