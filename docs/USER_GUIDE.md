# AXCSAS User Guide

**Getting Started with XRD Crystal Size Analysis**

Version: 0.2.0

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Basic Concepts](#basic-concepts)
4. [Quick Start Tutorial](#quick-start-tutorial)
5. [Common Workflows](#common-workflows)
6. [Understanding Results](#understanding-results)
7. [Customization](#customization)
8. [Troubleshooting](#troubleshooting)
9. [Best Practices](#best-practices)

---

## Introduction

### What is AXCSAS?

AXCSAS (Advanced XRD Crystallite Size Analysis System) is a Python package for analyzing X-ray diffraction (XRD) data to determine:

- **Crystallite size** using the Scherrer equation
- **Microstrain** using Williamson-Hall analysis  
- **Preferred orientation** using Harris texture coefficients
- **Stacking faults** and lattice parameter variations

### Who Should Use AXCSAS?

- Materials scientists studying thin films and nanoparticles
- Researchers analyzing electrodeposited copper
- Anyone performing XRD peak analysis for crystallite size determination

### Key Features

- Automated peak detection and Pseudo-Voigt fitting  
- Instrumental broadening correction (Caglioti equation)  
- Direction-dependent Scherrer K values  
- Comprehensive texture analysis  
- Publication-quality plots (2400 DPI)  
- Batch processing support

---

## Installation

### Prerequisites

- Python 3.8 or newer
- pip package manager

### Step-by-Step Installation

**1. Install from source** (recommended for latest version):

```bash
# Clone the repository
git clone https://github.com/Shuruyue/AXCSAS.git
cd AXCSAS

# Install in editable mode
pip install -e .
```

**2. Verify installation**:

```bash
# Check version
axcsas --version

# Test with sample data (if available)
axcsas analyze tests/data/sample.txt --output test_results/
```

### Dependencies

All dependencies are installed automatically:
- NumPy (numerical arrays)
- SciPy (peak fitting)
- Matplotlib (visualization)
- PyYAML (configuration)

---

## Basic Concepts

### XRD Fundamentals

**What is XRD?**  
X-ray Diffraction measures how X-rays scatter off crystal planes. The resulting pattern shows peaks at specific angles (2θ) that depend on:
- Crystal structure (e.g., FCC copper)
- Lattice spacing (d)
- X-ray wavelength (λ, typically Cu Kα = 1.54 Å)

**Bragg's Law**:
```
n·λ = 2d·sin(θ)
```

**Peak Width and Crystallite Size**:  
Smaller crystallites → broader peaks  
Larger crystallites → narrower peaks

### Key Parameters

| Parameter | Symbol | Typical Range | Meaning |
|-----------|--------|---------------|---------|
| 2θ | Two theta | 20-90° | Diffraction angle |
| FWHM | β | 0.1-1.0° | Peak width |
| Crystallite size | D | 5-200 nm | Coherent domain size |
| Microstrain | ε | 0-0.5% | Lattice distortion |
| Texture coefficient | TC | 0.5-3.0 | Orientation preference |

### The Scherrer Equation

```
D = (K·λ) / (β·cos(θ))
```

Where:
- D = crystallite size (nm)
- K = shape factor (0.89-0.94, depends on hkl)
- λ = X-ray wavelength (1.540562 Å for Cu Kα₁)
- β = FWHM corrected for instrumental broadening (radians)
- θ = Bragg angle (half of 2θ)

---

## Quick Start Tutorial

### Tutorial 1: Analyze a Single File

**Step 1**: Prepare your XRD data file
- Supported formats: `.txt`, `.xy`, `.csv`, `.raw`
- Must have two columns: 2θ (degrees) and Intensity (counts)

**Step 2**: Run analysis

```bash
axcsas analyze my_sample.txt --output results/
```

**What happens**:
1. Data is loaded and validated
2. Peaks are automatically detected
3. Each peak is fitted with Pseudo-Voigt profile
4. Crystallite sizes are calculated
5. Williamson-Hall analysis is performed
6. Texture coefficients are computed
7. Comprehensive report is generated

**Step 3**: Check results

```
results/
├── my_sample_report.txt          # Text report
├── my_sample_fwhm_evolution.png  # FWHM by peak
├── my_sample_size_evolution.png  # Size by peak
├── my_sample_wh_plot.png         # W-H linear fit
└── my_sample_texture_polar.png   # Texture radar plot
```

### Tutorial 2: Python Script

```python
from axcsas.analysis import run_full_analysis

# Run analysis
result = run_full_analysis("data/sample.txt")

# Print summary
print(f"Sample: {result.sample_name}")
print(f"Average crystallite size: {result.average_size_nm:.1f} nm")

# Check each peak
for scherrer in result.scherrer_results:
    hkl = scherrer.hkl
    size = scherrer.size_nm
    print(f"  ({hkl[0]}{hkl[1]}{hkl[2]}): {size:.1f} nm")

# Williamson-Hall results
if result.wh_result:
    print(f"\nW-H size: {result.wh_result.crystallite_size_nm:.1f} nm")
    print(f"Microstrain: {result.wh_result.microstrain*100:.3f}%")

# Texture analysis
if result.texture_result:
    print(f"\nDominant orientation: {result.texture_result.dominant_hkl}")
    print(f"TC value: {result.texture_result.dominant_tc:.2f}")
```

---

## Common Workflows

### Workflow A: Batch Processing Time Series

**Scenario**: You have XRD files from samples measured at 0h, 2h, 4h, 24h (self-annealing study).

```python
from axcsas.analysis import batch_analyze
from pathlib import Path
import pandas as pd

# Collect all files
files = sorted(Path("data/time_series").glob("*.txt"))

# Analyze all
results = batch_analyze(files)

# Extract time-dependent crystallite sizes
data = []
for result in results:
    avg_size = result.average_size_nm
    # Parse time from filename (e.g., "sample_2h.txt" → 2)
    time_hours = float(result.sample_name.split('_')[1].replace('h', ''))
    data.append({'time_hours': time_hours, 'size_nm': avg_size})

df = pd.DataFrame(data)
df.to_csv("size_vs_time.csv", index=False)
print(df)
```

### Workflow B: Texture Analysis

**Scenario**: Analyze preferred orientation in electrodeposited copper.

```python
from axcsas.methods import TextureAnalyzer
from axcsas.preprocessing import load_xrd_data
from axcsas.fitting import find_peaks, LMOptimizer
from axcsas.fitting import assign_hkl

# Load data
two_theta, intensity = load_xrd_data("sample.txt")

# Find and fit peaks
peaks = find_peaks(two_theta, intensity, min_height=100)

optimizer = LMOptimizer()
fitted_peaks = {}

for peak in peaks:
    # Assign (hkl)
    hkl = assign_hkl(peak.two_theta, tolerance=0.5)
    if hkl:
        # Fit  to get accurate integrated area
        fit_result = optimizer.fit_single_peak(
            two_theta, intensity, peak_idx=peak.index
        )
        if fit_result.success:
            fitted_peaks[hkl] = fit_result.area

# Texture analysis
analyzer = TextureAnalyzer(use_area=True)
texture_result = analyzer.analyze(fitted_peaks)

# Interpret results
print("=== Texture Analysis ===")
for tc_detail in texture_result.tc_details:
    hkl_str = f"({tc_detail.hkl[0]}{tc_detail.hkl[1]}{tc_detail.hkl[2]})"
    print(f"{hkl_str}: TC={tc_detail.tc_value:.2f} - {tc_detail.orientation_type.value}")

if texture_result.is_random:
    print("\n[INFO] Random orientation (no preferred texture)")
else:
    print(f"\n[WARNING] Preferred orientation along {texture_result.dominant_hkl}")
```

---

## Understanding Results

### Scherrer Results

**Sample Output**:
```
Peak (111): Size = 28.5 nm [VALID]
  K-factor: 0.94
  FWHM (observed): 0.352°
  FWHM (corrected): 0.345°
  β_obs/β_inst = 6.4 (reliable)
```

**Interpretation**:
- **Size**: Average coherent domain size along (111) direction
- **VALID flag**: β_obs/β_inst > 1.2, correction is reliable
- **UNRELIABLE flag**: Instrumental broadening dominates → size > 200 nm

### Williamson-Hall Results

**Sample Plot**:
```
β·cos(θ) = (K·λ/D) + (4ε·sin(θ))
           ↑ intercept  ↑ slope
```

**Interpretation**:
- **Positive slope**: Microstrain present
- **Zero slope**: Pure size broadening
- **Negative slope**: Likely data issue or stacking faults
- **High R² (>0.95)**: Good linear fit, results reliable
- **Low R² (<0.90)**: Anisotropic broadening or measurement error

### Texture Coefficients (TC)

**Harris TC Formula**:
```
TC(hkl) = [I(hkl)/I₀(hkl)] / [Average of all I/I₀]
```

**Interpretation**:
| TC Range | Meaning |
|----------|---------|
| 0.9-1.1 | Random orientation |
| > 1.5 | Strong (hkl) texture |
| < 0.9 | Suppressed orientation |

**Example**:
```
(111): TC = 1.85 → Strong (111) fiber texture
(200): TC = 0.65 → (200) suppressed
(220): TC = 0.98 → Random
```

---

## Customization

### Custom Configuration File

Create `my_config.yaml`:

```yaml
# Instrumental parameters (from LaB6 calibration)
instrument:
  caglioti:
    U: 0.001
    V: -0.002
    W: 0.003

# Preprocessing options
preprocessing:
  smoothing:
    enable: true
    window_size: 11
    poly_order: 3
  
  background:
    enable: true
    method: chebyshev  # or 'sonneveld_visser'
    poly_degree: 5
  
  kalpha_strip:
    enable: true

# Peak fitting
peak_fitting:
  min_intensity: 100
  peak_window: 2.0  # Search window in degrees

# Quality thresholds
validation:
  max_rwp: 10.0        # Maximum R_wp (%)
  min_r_squared: 0.95  # Minimum R²
```

**Use in Python**:

```python
from axcsas.core import load_config
from axcsas.analysis import AXCSASPipeline, AnalysisConfig

# Load config
config_dict = load_config("my_config.yaml")

# Create analysis config
analysis_config = AnalysisConfig(
    caglioti_u=config_dict['instrument']['caglioti']['U'],
    caglioti_v=config_dict['instrument']['caglioti']['V'],
    caglioti_w=config_dict['instrument']['caglioti']['W'],
)

# Run with custom config
pipeline = AXCSASPipeline(analysis_config)
result = pipeline.analyze("sample.txt")
```

---

## Troubleshooting

### Common Errors

**Error: "Only 1 peaks found"**

**Cause**: Peak intensity too low or peaks outside search range

**Solution**:
1. Check your data quality (signal-to-noise ratio)
2. Adjust `min_intensity` threshold
3. Verify 2θ range covers expected peaks (20-90° for Cu)

---

**Error: "Broadening ratio below threshold"**

**Cause**: β_obs/β_inst < 1.2, crystallites too large (>200 nm)

**Solution**:
1. Check instrumental broadening parameters (U, V, W)
2. Verify peak fitting quality (R² should be >0.95)
3. Consider TEM for large crystallites

---

**Error: "Negative FWHM² after correction"**

**Cause**: Observed FWHM smaller than instrumental

**Solutions**:
1. Re-calibrate instrumental broadening with LaB6
2. Check if peak is correctly identified
3. Verify data preprocessing didn't over-smooth

---

### Tips for Best Results

**Data Quality**:
- Use slow scan speed (0.02-0.05 deg/step)
- Ensure good counting statistics (>1000 counts at peak)
- Minimize background (proper sample preparation)

**Sample Preparation**:
- Flat, smooth surface
- Sufficient thickness (>1 um for Cu)
- Avoid preferred orientation (for powder XRD)

**Analysis Parameters**:
- Always calibrate with LaB6 standard first
- Use integrated area (not peak height) for texture
- Check R-squared values (should be >0.95)

---

## Next Steps

**For Advanced Users**:
- See [API Reference](API_REFERENCE.md) for detailed function documentation
- Check `examples/` folder for complete workflows
- Customize analysis pipeline for specific materials

**For Developers**:
- See `CONTRIBUTING.md` for development guidelines
- Run tests with `pytest tests/`
- Check code style with `ruff check .`

---

**Need Help?**  
- GitHub Issues: https://github.com/Shuruyue/AXCSAS/issues
- Documentation: See `docs/` folder

**Citation**:  
If you use AXCSAS in your research, please cite:
```bibtex
@software{axcsas2026,
  title = {AXCSAS: Advanced XRD Crystallite Size Analysis System},
  author = {[Your Name]},
  year = {2026},
  url = {https://github.com/Shuruyue/AXCSAS}
}
```
