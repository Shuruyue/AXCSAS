# AXCSAS Verification & Utility Scripts
## AXCSAS 驗證與公用腳本

This directory contains standalone scripts for validating physical calculations, testing instrument calibration, and generating debug plots.
此目錄包含用於驗證物理計算、測試儀器校準和生成除錯圖表的獨立腳本。

### Verification Scripts (驗證腳本)

#### `verify_physics.py`
**Purpose**: Validates basic crystallographic formulas.
**Uses**:
- Checks Bragg's Law calculations.
- Verifies interplanar spacing ($d_{hkl}$) for Cubic systems.
- Confirms wavelength constants (Cu K$\alpha$).

#### `verify_elastic_moduli.py`
**Purpose**: Verification of Young's Modulus calculations.
**Uses**:
- Calculates directional Young's Modulus ($E_{hkl}$) for Copper.
- Verifies Reuss/Voigt/Hill averaging methods.

### Plot Generation (繪圖工具)

#### `generate_fitting_diagnosis.py`
**Purpose**: Generates visual diagnostics for peak fitting.
**Output**:
- `outputs/fitting_diagnosis/`: Individual peak fits and residual plots.

#### `generate_fwhm_plots.py`
**Purpose**: Plots Full Width at Half Maximum (FWHM) evolution.
**Output**:
- `outputs/fwhm_plots/`: Williamson-Hall style FWHM vs $sin(\theta)$ plots.

### Usage Example

Run scripts from the project root directory:

```bash
# Run physics verification
python scripts/verify_physics.py

# Generate diagnosis plots
python scripts/generate_fitting_diagnosis.py
```
