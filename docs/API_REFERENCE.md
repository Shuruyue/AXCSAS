# AXCSAS API Documentation

**Advanced XRD Crystallite Size Analysis System**

Version: 0.2.0

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Core Modules](#core-modules)
4. [Preprocessing](#preprocessing)
5. [Fitting](#fitting)
6. [Analysis Methods](#analysis-methods)
7. [Validation](#validation)
8. [Visualization](#visualization)
9. [CLI Usage](#cli-usage)

---

## Installation

```bash
# From source
git clone https://github.com/Shuruyue/AXCSAS.git
cd AXCSAS
pip install -e .

# Or install from PyPI (when available)
pip install axcsas
```

### Requirements

- Python ≥ 3.8
- NumPy ≥ 1.20
- SciPy ≥ 1.7
- Matplotlib ≥ 3.5
- PyYAML ≥ 6.0

---

## Quick Start

### Basic Analysis

```python
from axcsas.analysis import run_full_analysis

# Analyze a single XRD file
result = run_full_analysis("data/sample.txt")

# Print comprehensive report
print(result.report)

# Access specific results
print(f"Average crystallite size: {result.average_size_nm:.1f} nm")
print(f"Dominant texture: {result.texture_result.dominant_hkl}")
```

### Batch Processing

```python
from axcsas.analysis import batch_analyze
from pathlib import Path

# Process all files in directory
files = list(Path("data/raw").glob("*.txt"))
results = batch_analyze(files)

# Generate summary CSV
from axcsas.analysis.report_generator import generate_csv_summary
generate_csv_summary(results, "output/summary.csv")
```

---

## Core Modules

### `axcsas.core.constants`

Physical constants and empirical thresholds for XRD analysis.

**Key Constants**:

```python
from axcsas.core.constants import (
    CU_KA1,           # 1.540562 Å (Bearden 1967)
    CU_KA2,           # 1.544390 Å
    KA2_KA1_RATIO,    # 0.5 (Burger-Dorgelo rule)
    SCHERRER_K,       # 0.94 (spherical particles)
    MIN_RELIABLE_SIZE, # 2.0 nm
    MAX_RELIABLE_SIZE, # 200.0 nm
)
```

**User-Adjustable Thresholds**:

```python
# Modify for different quality standards
from axcsas.core import constants

# For high-precision publication work
constants.MIN_R_SQUARED = 0.99
constants.MAX_RWP_PERCENT = 5.0

# For routine analysis
constants.MIN_R_SQUARED = 0.90
constants.MAX_RWP_PERCENT = 15.0
```

### `axcsas.core.copper_crystal`

Copper-specific crystallographic parameters with full academic citations.

```python
from axcsas.core.copper_crystal import (
    CU_CRYSTAL,           # Crystal structure
    CU_JCPDS_EXTENDED,    # JCPDS 04-0836 peak data
    get_scherrer_k,       # Direction-dependent K values
    get_poisson_ratio,    # Direction-dependent Poisson ratio
    get_elastic_modulus,  # Direction-dependent E values
)

# Get K value for specific direction
k_111 = get_scherrer_k(1, 1, 1)  # 0.94 for (111)
k_200 = get_scherrer_k(2, 0, 0)  # 0.89 for (200)

# Access JCPDS peak data
peak_111 = CU_JCPDS_EXTENDED[(1, 1, 1)]
print(f"2θ: {peak_111['two_theta']:.3f}°")
print(f"d: {peak_111['d_spacing']:.4f} Å")
print(f"I: {peak_111['intensity']}")
```

### `axcsas.core.config_loader`

Load and merge YAML configuration files.

```python
from axcsas.core.config_loader import load_config

# Load with custom config
config = load_config("my_config.yaml")

# Access nested parameters
smoothing_window = config['preprocessing']['smoothing']['window_size']
```

### `axcsas.core.units`

Unit conversion utilities for XRD calculations.

```python
from axcsas.core.units import (
    degrees_to_radians,
    angstrom_to_nm,
    two_theta_to_d,      # Bragg's law
    d_to_two_theta,
)

# Convert 2θ to d-spacing
d_spacing = two_theta_to_d(43.3, wavelength=1.540562)  # Å
```

---

## Preprocessing

### Data Loading

```python
from axcsas.preprocessing import load_xrd_data

# Auto-detect format (.xy, .csv, .txt, .raw)
two_theta, intensity = load_xrd_data("sample.txt")
```

### Smoothing

```python
from axcsas.preprocessing import SavitzkyGolayFilter

filter = SavitzkyGolayFilter(window_size=11, poly_order=3)
smoothed = filter.apply(intensity)
```

### Background Subtraction

```python
from axcsas.preprocessing import BackgroundSubtractor

# Chebyshev polynomial method
subtractor = BackgroundSubtractor(method="chebyshev", poly_degree=5)
corrected, background = subtractor.subtract(two_theta, intensity)

# Or Sonneveld-Visser second derivative
subtractor = BackgroundSubtractor(method="sonneveld_visser")
corrected, background = subtractor.subtract(two_theta, intensity)
```

### Complete Pipeline

```python
from axcsas.preprocessing import PreprocessingPipeline

pipeline = PreprocessingPipeline(
    enable_smoothing=True,
    enable_background=True,
    enable_kalpha_strip=True
)

result = pipeline.run(two_theta, intensity)
print(result.summary())
```

---

## Fitting

### Peak Detection

```python
from axcsas.fitting import find_peaks

peaks = find_peaks(two_theta, intensity, min_height=100, min_distance=0.5)

for peak in peaks:
    print(f"Peak at {peak.two_theta:.3f}°, FWHM={peak.estimated_fwhm:.3f}°")
```

### Pseudo-Voigt Fitting

```python
from axcsas.fitting import LMOptimizer, PseudoVoigtParams

optimizer = LMOptimizer()

# Auto-fit single peak
result = optimizer.fit_single_peak(theta_region, intensity_region)

if result.success:
    print(f"Center: {result.params.center:.3f}°")
    print(f"FWHM: {result.params.fwhm:.3f}°")
    print(f"η: {result.params.eta:.3f}")
    print(f"R²: {result.r_squared:.4f}")
```

### Kα Doublet Fitting

```python
from axcsas.fitting.ka_doublet import DoubletFitter

fitter = DoubletFitter()
result = fitter.fit(theta_region, intensity_region)

print(f"Kα₁: {result.center_ka1:.3f}°")
print(f"Kα₂: {result.center_ka2:.3f}°")
print(f"FWHM: {result.fwhm:.3f}° (intrinsic)")
```

### HKL Assignment

```python
from axcsas.fitting import assign_hkl, assign_hkl_detailed

# Simple assignment
hkl = assign_hkl(43.3, tolerance=0.5)  # (1, 1, 1)

# Detailed assignment with confidence
assignment = assign_hkl_detailed(43.3)
print(f"{assignment.hkl}: deviation={assignment.deviation:.3f}°")
print(f"Confidence: {assignment.confidence}")
```

---

## Analysis Methods

### Scherrer Analysis

```python
from axcsas.methods import ScherrerCalculator

calc = ScherrerCalculator()

result = calc.calculate(
    two_theta=43.3,
    fwhm_observed=0.35,
    fwhm_instrumental=0.055,
    hkl=(1, 1, 1)
)

print(f"Crystallite size: {result.size_nm:.1f} nm")
print(f"K-factor: {result.k_factor:.3f}")
print(f"Validity: {result.validity_flag.value}")
```

### Williamson-Hall Analysis

```python
from axcsas.methods import WilliamsonHallAnalyzer

analyzer = WilliamsonHallAnalyzer()

# Multiple peaks required
two_theta_arr = [43.3, 50.5, 74.2]
fwhm_arr = [0.35, 0.38, 0.45]
hkl_list = [(1,1,1), (2,0,0), (2,2,0)]

result = analyzer.analyze(two_theta_arr, fwhm_arr, hkl_list)

print(f"Crystallite size: {result.crystallite_size_nm:.1f} nm")
print(f"Microstrain: {result.microstrain*100:.3f}%")
print(f"R²: {result.r_squared:.4f}")
print(f"Quality: {result.quality_level.value}")
```

### Texture Analysis

```python
from axcsas.methods import TextureAnalyzer

analyzer = TextureAnalyzer()

# Integrated intensities for each peak
intensities = {
    (1, 1, 1): 15680,
    (2, 0, 0): 5520,
    (2, 2, 0): 4200,
}

result = analyzer.analyze(intensities)

print(f"Dominant: {result.dominant_hkl}, TC={result.dominant_tc:.2f}")
print(f"Random texture: {result.is_random}")

for detail in result.tc_details:
    print(f"{detail.hkl}: TC={detail.tc_value:.2f} [{detail.orientation_type.value}]")
```

### Caglioti Correction

```python
from axcsas.methods import CagliotiCorrection

correction = CagliotiCorrection(U=0.0, V=0.0, W=0.003)

# Calculate instrumental broadening at 2θ
fwhm_inst = correction.calculate_fwhm_inst(43.3)

# Correct observed FWHM
fwhm_sample, is_reliable, warning = correction.correct_broadening(
    fwhm_observed=0.35,
    two_theta=43.3
)

if is_reliable:
    print(f"Corrected FWHM: {fwhm_sample:.4f}°")
else:
    print(f"WARNING: {warning}")
```

---

## Validation

### Error Analysis

```python
from axcsas.validation import ErrorAnalyzer

analyzer = ErrorAnalyzer()

# Validate crystallite size
size_result = analyzer.validate_size(size_nm=45.2)
if not size_result.is_valid:
    for warning in size_result.warnings:
        print(warning)

# Validate broadening ratio
br_result = analyzer.validate_broadening(
    fwhm_observed=0.35,
    fwhm_instrumental=0.055
)

# Comprehensive validation
result = analyzer.validate_all(
    size_nm=45.2,
    fwhm_observed=0.35,
    fwhm_instrumental=0.055,
    rwp=7.2,
    r_squared=0.9850
)
```

### Goodness of Fit

```python
from axcsas.validation import assess_fit_quality

quality = assess_fit_quality(
    observed=intensity_data,
    calculated=fitted_curve,
    n_params=4,
    rwp_threshold=10.0
)

print(f"Rwp: {quality.rwp:.2f}%")
print(f"R²: {quality.r_squared:.4f}")
print(f"Acceptable: {quality.is_acceptable}")
```

---

## Visualization

### FWHM Plots

```python
from axcsas.visualization import plot_fwhm_evolution

plot_fwhm_evolution(
    results=analysis_results,
    x_param='time',  # or 'concentration'
    output_path='fwhm_evolution.png'
)
```

### Scherrer Plots

```python
from axcsas.visualization import plot_scherrer_evolution_by_peak

plot_scherrer_evolution_by_peak(
    results=analysis_results,
    x_param='time',
    output_path='size_evolution.png'
)
```

### Williamson-Hall Plot

```python
from axcsas.visualization import plot_williamson_hall

plot_williamson_hall(
    wh_result=wh_result,
    output_path='wh_plot.png'
)
```

### Texture Polar Plot

```python
from axcsas.visualization import plot_texture_polar

plot_texture_polar(
    texture_result=texture_result,
    output_path='texture_polar.png',
    sample_name='Sample_01'
)
```

---

## CLI Usage

### Analyze XRD Data

```bash
# Single file analysis
axcsas analyze sample.txt --output results/

# Batch analysis
axcsas analyze data/*.txt --output results/

# Generate comprehensive report
axcsas report data/*.txt --output report/
```

### Instrument Calibration

```bash
# Calibrate using LaB6 standard
axcsas calibrate lab6_standard.txt --output calibration/

# The calibration results (U, V, W) can be used in config.yaml
```

### Configuration

Create `config.yaml`:

```yaml
instrument:
  caglioti:
    U: 0.0
    V: 0.0
    W: 0.003

preprocessing:
  smoothing:
    enable: true
    window_size: 11
    poly_order: 3
  
  background:
    enable: true
    method: chebyshev
    poly_degree: 5

peak_fitting:
  min_intensity: 100
  peak_window: 2.0

validation:
  max_rwp: 10.0
  min_r_squared: 0.95
```

---

## API Reference by Module

### Complete Module Index

| Module | Purpose |
|--------|---------|
| `axcsas.core.constants` | Physical constants and thresholds |
| `axcsas.core.copper_crystal` | Copper crystallography |
| `axcsas.core.config_loader` | Configuration management |
| `axcsas.core.units` | Unit conversions |
| `axcsas.preprocessing` | Data preprocessing pipeline |
| `axcsas.fitting` | Peak detection and fitting |
| `axcsas.methods.scherrer` | Scherrer equation |
| `axcsas.methods.williamson_hall` | W-H strain analysis |
| `axcsas.methods.texture` | Harris texture coefficients |
| `axcsas.methods.caglioti` | Instrumental correction |
| `axcsas.methods.defect_analysis` | Stacking faults & stress |
| `axcsas.validation` | Result validation |
| `axcsas.visualization` | Plotting functions |
| `axcsas.analysis` | Complete pipeline |

---

## References

1. **Scherrer, P.** (1918). Bestimmung der Grösse und der inneren Struktur von Kolloidteilchen mittels Röntgenstrahlen. *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen*, 98-100.

2. **Bearden, J. A.** (1967). X-Ray Wavelengths. *Reviews of Modern Physics*, 39(1), 78-124.

3. **Williamson, G. K., & Hall, W. H.** (1953). X-ray line broadening from filed aluminium and wolfram. *Acta Metallurgica*, 1(1), 22-31.

4. **Harris, G. B.** (1952). Quantitative measurement of preferred orientation in rolled uranium bars. *Philosophical Magazine*, 43(336), 113-123.

5. **Caglioti, G., Paoletti, A., & Ricci, F. P.** (1958). Choice of collimators for a crystal spectrometer for neutron diffraction. *Nuclear Instruments*, 3(4), 223-228.

---

**For more information**: See [GitHub Repository](https://github.com/Shuruyue/AXCSAS)
