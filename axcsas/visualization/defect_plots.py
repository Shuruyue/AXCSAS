"""
Defect Analysis Plots Module
============================

Visualization for residual stress and stacking fault analysis.
Features:
- Residual Stress Evolution (by Annealing Time)
- Stacking Fault Probability Evolution
- Lattice Constant Monitoring
- Automatic footnote for physical limitations

Style matches: axcsas.visualization.style
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path

from axcsas.visualization.style import apply_axcsas_style, save_figure, get_color_palette

# Physical Limitation Footnotes (Bilingual)
LIMITATION_NOTE_STRESS = (
    "Note: Stress calculated using direction-dependent E & ν (plane stress model). "
    "Values include potential chemical expansion (impurity effect).\n"
    "註：應力採方向相依模數計算(平面應力模形)。數值包含潛在的化學膨脹效應(雜質影響)。"
)

LIMITATION_NOTE_SF = (
    "Note: Stacking fault probability α estimated via Warren method (peak separation). "
    "Coupling with anisotropic stress may affect absolute values.\n"
    "註：層錯機率α採Warren法(峰間距)估算。各向異性應力耦合可能影響絕對數值。"
)


def plot_stress_evolution(
    data: List[Dict[str, Any]],
    output_path: Optional[str] = None,
    dpi: int = 1200,
    show: bool = True
) -> plt.Figure:
    """
    Plot Residual Stress Evolution for (111) and (200) grains.
    
    Args:
        data: List of sample dictionaries containing 'stress_results'
        output_path: Path to save the figure
    """
    apply_axcsas_style()
    
    # Extract data
    high_quality_data = [s for s in data if s.get('high_quality', True)]
    concentrations = sorted(set(s.get('concentration', 0) for s in high_quality_data))
    
    # Create 2x2 grid (similar to texture plots)
    n_conc = len(concentrations)
    n_cols = 2
    n_rows = (n_conc + 1) // 2
    
    figsize = (14, 11)
    if n_rows > 2:
        figsize = (14, 11 * (n_rows/2))
        
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()
    
    # Colors for directions
    # (111) = Blue, (200) = Orange
    # Colors for directions
    # (111) = Blue, (200) = Orange, (220) = Green
    colors = {'(111)': '#0072B2', '(200)': '#D55E00', '(220)': '#009E73'} 
    markers = {'(111)': 'o', '(200)': 's', '(220)': '^'}
    
    global_y_min = 0
    global_y_max = 0
    
    for idx, conc in enumerate(concentrations):
        ax = axes[idx]
        conc_samples = [s for s in high_quality_data if s.get('concentration') == conc]
        
        # Determine unique x-axis (times)
        times = sorted(list(set(s['time'] for s in conc_samples)))
        
        for hkl_key in ['(111)', '(200)', '(220)']:
            x_vals = []
            y_vals = [] # Stress in MPa
            
            for t in times:
                sample = next((s for s in conc_samples if s['time'] == t), None)
                if sample and 'stress_results' in sample:
                    # Find stress result for this HKL
                    # stress_results structure expected: {'(111)': {val: ..., err: ...}, ...}
                    res = sample['stress_results'].get(hkl_key)
                    if res:
                        x_vals.append(t)
                        y_vals.append(res['stress_mpa'])
                        
                        # Track global range
                        global_y_min = min(global_y_min, res['stress_mpa'])
                        global_y_max = max(global_y_max, res['stress_mpa'])
            
            if x_vals:
                ax.plot(x_vals, y_vals, 
                       marker=markers[hkl_key], 
                       color=colors[hkl_key],
                       linestyle='-', 
                       linewidth=2, 
                       markersize=8,
                       label=f"{hkl_key} Stress",
                       alpha=0.9)
                
                # Zero line
                ax.axhline(0, color='gray', linestyle=':', alpha=0.5)

        ax.set_title(f"{conc} mL/1.5L", fontsize=13, fontweight='bold')
        ax.set_xlabel("Annealing Time (hours)", fontsize=12)
        ax.set_ylabel("Residual Stress (MPa)", fontsize=12)
        ax.legend(loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3)
    
    # Padding for Y-axis
    y_range = global_y_max - global_y_min
    if y_range == 0: y_range = 100
    
    # Apply consistent Y-limits
    for ax in axes[:len(concentrations)]:
        ax.set_ylim(global_y_min - y_range*0.1, global_y_max + y_range*0.15)
        
    # Hide unused
    for idx in range(len(concentrations), len(axes)):
        axes[idx].set_visible(False)
        
    fig.suptitle("Residual Stress Evolution (Anisotropic Model)", fontsize=16, fontweight='bold', y=0.98)
    
    # Add Limitation Note at bottom
    # Add Limitation Note at bottom
    # fig.text(0.5, 0.01, LIMITATION_NOTE_STRESS, ha='center', va='bottom', 
    #          fontsize=10, style='italic', color='#555555', 
    #          bbox=dict(facecolor='#f0f0f0', alpha=0.5, edgecolor='none', pad=5))
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.96]) # Adjust for footer
    
    if output_path:
        save_figure(fig, output_path, dpi=dpi)
        
    if show:
        plt.show()
    
    return fig


def plot_stacking_fault_evolution(
    data: List[Dict[str, Any]],
    output_path: Optional[str] = None,
    dpi: int = 1200,
    show: bool = True
) -> plt.Figure:
    """
    Plot Stacking Fault Probability Evolution.
    """
    apply_axcsas_style()
    
    high_quality_data = [s for s in data if s.get('high_quality', True)]
    concentrations = sorted(set(s.get('concentration', 0) for s in high_quality_data))
    
    n_conc = len(concentrations)
    n_cols = 2
    n_rows = (n_conc + 1) // 2
    
    figsize = (14, 11)
    if n_rows > 2:
        figsize = (14, 11 * (n_rows/2))
        
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()
    
    # Color mapping for concentration consistency (though we plot specific concs in panels)
    # Actually, here we just need one line per panel, so standard color is fine.
    line_color = '#CC79A7' # Reddish-purple from colorblind palette
    
    for idx, conc in enumerate(concentrations):
        ax = axes[idx]
        conc_samples = [s for s in high_quality_data if s.get('concentration') == conc]
        times = sorted(list(set(s['time'] for s in conc_samples)))
        
        x_vals = []
        y_vals = [] # Alpha %
        
        for t in times:
            sample = next((s for s in conc_samples if s['time'] == t), None)
            if sample and 'stacking_fault' in sample:
                sf_res = sample['stacking_fault']
                # Expecting {'alpha_percent': ...}
                if sf_res and sf_res.get('alpha_percent') is not None:
                    x_vals.append(t)
                    y_vals.append(sf_res.get('alpha_percent'))
        
        if x_vals:
            ax.plot(x_vals, y_vals, 
                   marker='D', 
                   color=line_color,
                   linestyle='-', 
                   linewidth=2, 
                   markersize=7,
                   label=f"Stacking Fault %",
                   alpha=0.9)
            
        ax.set_title(f"{conc} mL/1.5L", fontsize=13, fontweight='bold')
        ax.set_xlabel("Annealing Time (hours)", fontsize=12)
        ax.set_ylabel("Stacking Fault Probability α (%)", fontsize=12)
        
        # Set realistic limits (0 to 2% is typical for Cu)
        current_max = max(y_vals) if y_vals else 0.5
        y_limit = max(1.0, current_max * 1.2)
        ax.set_ylim(0, y_limit)
        
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
        
    # Hide unused
    for idx in range(len(concentrations), len(axes)):
        axes[idx].set_visible(False)
        
    fig.suptitle("Stacking Fault Probability Evolution (Warren Method)", fontsize=16, fontweight='bold', y=0.98)
    
    # Add Limitation Note
    # Add Limitation Note
    # fig.text(0.5, 0.01, LIMITATION_NOTE_SF, ha='center', va='bottom', 
    #          fontsize=10, style='italic', color='#555555', 
    #          bbox=dict(facecolor='#f0f0f0', alpha=0.5, edgecolor='none', pad=5))
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.96])
    
    if output_path:
        save_figure(fig, output_path, dpi=dpi)
        
    if show:
        plt.show()
        
    return fig
