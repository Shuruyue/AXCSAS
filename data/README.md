# Data Directory

This directory contains XRD data files.

## Structure

- `raw/` - Original XRD data files (.xy, .csv, .txt)
- `standards/` - NIST standard reference materials (SRM 660c, SRM 640e)
- `processed/` - Preprocessed data (after smoothing, background subtraction)

## Supported Formats

- `.xy` - Two-column format (2θ, Intensity)
- `.csv` - CSV with header
- `.txt` - Bruker TXT export format
