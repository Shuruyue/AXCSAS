# AXCSAS

## Advanced XRD Crystallite Size Analysis System

Automated crystallite size analysis system using Pseudo-Voigt fitting and instrumental correction.

---

## Project Overview

AXCSAS is an automated XRD analysis system that addresses common issues in manual Scherrer size calculations:

- Subjective baseline selection errors
- Inaccurate instrumental broadening correction
- Incorrect peak profile function selection

### Core Features

- **Pseudo-Voigt Full Spectrum Fitting**: Uses the most physically accurate mathematical description of XRD peak profiles
- **Caglioti Instrumental Correction**: Full-angle instrumental width correction using NIST SRM 660c (LaB6)
- **High Precision**: Applicable to crystallite sizes in the 2-100 nm range

---

## Quick Start

### Installation

```bash
pip install -e .
```

Or for development:

```bash
pip install -e ".[dev]"
```

### CLI Usage

```bash
# Run analysis on a single file
axcsas analyze data/raw/sample.xy

# Generate comprehensive report
axcsas report --input-dir data/raw/ --output-dir outputs/
```

### Validation Scripts

```bash
# Verify 2-theta peak positions (Bragg's Law)
python scripts/verify_physics.py

# Verify directional Young's modulus (Ledbetter & Naimon)
python scripts/verify_elastic_moduli.py
```

---

## Project Structure

AXCSAS/
├── config.yaml              # Global configuration
├── pyproject.toml           # Python project configuration
├── data/                    # Data directory
│   ├── raw/                 # Raw XRD data
│   ├── standards/           # NIST standard data
│   └── processed/           # Preprocessed data
├── axcsas/                  # Core source code
│   ├── core/                # Physical constants and crystallography
│   ├── preprocessing/       # Data preprocessing modules
│   ├── fitting/             # Peak fitting core
│   ├── methods/             # Analysis methods (Scherrer, W-H, Texture)
│   ├── validation/          # Error analysis and validation
│   └── visualization/       # Plotting and reports
├── scripts/                 # Verification & Utility scripts
├── tests/                   # Unit tests
├── docs/                    # User documentation (API, Guides)
├── dev_notes/               # Developer notes & technical details
├── references/              # Literature references (PDFs)
└── outputs/                 # Output directory (Ignored by Git)

---

## Theoretical Background

### Pseudo-Voigt Peak Profile

I(2theta) = I0 * [ eta * L(2theta) + (1-eta) * G(2theta) ] + Background

where L is Lorentzian, G is Gaussian, and eta is the mixing parameter.

### Caglioti Equation

FWHM_inst^2 = U * tan^2(theta) + V * tan(theta) + W

### Scherrer Equation

D = K * lambda / (beta * cos(theta))

---

## References

1. Bearden, J.A. (1967). X-Ray Wavelengths. Rev. Mod. Phys. 39, 78-124.
2. Langford, J.I. & Wilson, A.J.C. (1978). Scherrer after Sixty Years. J. Appl. Cryst. 11, 102-113.
3. NIST Standard Reference Material 660c (LaB6)

---

## License

MIT License
