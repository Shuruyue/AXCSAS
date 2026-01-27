"""
Generate comprehensive fitting diagnosis plots
生成綜合峰擬合診斷圖
"""

from pathlib import Path
from axcsas.analysis import batch_analyze
from axcsas.visualization.generate_fitting_diagnosis import generate_sample_fitting_plot
import matplotlib.pyplot as plt
import numpy as np

from axcsas.core.config_loader import load_config

# Setup paths
project_root = Path(__file__).resolve().parent.parent
config_path = project_root / "config.yaml"
config = load_config(config_path)

data_dir = project_root / "data/raw/202511"

# Use standardized output path from config
# Output directory
plots_root = project_root / config["output"]["paths"]["plots"]
dir_diagnosis = plots_root / "diagnosis"  # User requested removal of '03_' prefix
dir_diagnosis.mkdir(parents=True, exist_ok=True)

# Individual fits
GENERATE_INDIVIDUAL = True
dir_individual = dir_diagnosis / "individual_fits" # User requested removal of "failed"
if GENERATE_INDIVIDUAL:
    dir_individual.mkdir(parents=True, exist_ok=True)
    # Cleanup old diagnosis plots
    for f in dir_individual.glob("*.png"):
        f.unlink()
    print(f"Cleaned up old diagnosis plots in {dir_individual}")

# Get all data files
data_files = sorted(data_dir.glob("*.txt"))
print(f"Found {len(data_files)} data files")

# Analyze all samples
print("Analyzing samples for fitting diagnosis...")
results = batch_analyze([str(f) for f in data_files])

# Generate individual fitting diagnosis plots
print(f"\nGenerating {len(results)} individual fitting diagnosis plots...")
for i, result in enumerate(results, 1):
    print(f"  {i}/{len(results)}: {result.sample_name}")
    
    # Generate diagnosis plot for this sample
    output_file = dir_individual / f"{result.sample_name}_fitting.png"
    
    try:
        generate_sample_fitting_plot(
            Path(result.filepath), # Pass Path object as expected by the new function
            output_dir=dir_individual
        )
    except Exception as e:
        print(f"    ⚠ Error: {e}")

print(f"\n✓ Individual plots saved to: {dir_individual}")

# Generate summary statistics plot
# Summary plot generation disabled per user request
# print("\nGenerating fitting quality summary...")
# ... (code removed)

print(f"\n{'='*60}")
print("Fitting Diagnosis Complete!")
print(f"{'='*60}")
print(f"Individual plots: {len(results)} files in {dir_individual}")
# print(f"Summary plots: {summary_dir}")
