"""
Generate FWHM plots from XRD analysis results
使用 XRD 分析結果生成 FWHM 圖表
"""

from pathlib import Path
from axcsas.analysis import batch_analyze
from axcsas.visualization import plot_fwhm_evolution, plot_fwhm_by_concentration

from axcsas.core.config_loader import load_config

# Setup paths
project_root = Path(__file__).resolve().parent.parent
config_path = project_root / "config.yaml"
config = load_config(config_path)

# Find all data files
data_dir = project_root / "data/raw/202511"

# Use standardized output path from config
output_root = Path(config["output"]["paths"]["root"])
# Convert string path to Path object if needed
if isinstance(output_root, str):
    output_root = project_root / output_root

# FWHM plots go to outputs/plots/analysis (or similar, based on plan: outputs/plots/analysis)
# Plan said: outputs/plots/analysis
# Setup output subdirectories
plots_root = project_root / config["output"]["paths"]["plots"]
dir_results = plots_root / "results"
dir_results.mkdir(parents=True, exist_ok=True)

# Cleanup old plots
for f in dir_results.glob("*.png"):
    f.unlink()
print(f"Cleaned up old files in {dir_results}")

# Get all .txt files
data_files = sorted(data_dir.glob("*.txt"))
print(f"Found {len(data_files)} data files")

# Run batch analysis
print("Analyzing all samples...")
results = batch_analyze([str(f) for f in data_files])

# Prepare data for plotting
plot_data = []
for result in results:
    sample_data = {
        'name': result.sample_name,
        'concentration': result.leveler_concentration or 0,
        'time': result.sample_age_hours or 0,
        'peaks': []
    }
    
    # Add peak data
    for peak in result.peaks:
        hkl_str = f"({peak.hkl[0]}{peak.hkl[1]}{peak.hkl[2]})"
        sample_data['peaks'].append({
            'hkl': hkl_str,
            'fwhm': peak.fwhm,
            'fwhm_error': 0.01  # Default uncertainty
        })
    
    plot_data.append(sample_data)

print(f"Prepared data for {len(plot_data)} samples")

# 1. FWHM evolution by concentration
print("1. FWHM vs Time (grouped by concentration)...")
fig1 = plot_fwhm_evolution(
    plot_data,
    x_param='time',
    output_path=str(dir_results / "fwhm_evolution.png"),
    show=False
)
print("   Saved: fwhm_evolution.png")

# 2. FWHM by concentration (2x2 grid)
print("2. FWHM evolution grid (one subplot per concentration)...")
fig2 = plot_fwhm_by_concentration(
    plot_data,
    output_path=str(dir_results / "fwhm_compare_grid.png"),
    dpi=2000,
    show=False
)
print("   Saved: fwhm_compare_grid.png")

print(f"\n✓ All plots saved to: {plots_root}")
print(f"  - results/fwhm_evolution.png")
print(f"  - results/fwhm_compare_grid.png")

