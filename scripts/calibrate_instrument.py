#!/usr/bin/env python3
"""
Instrument Calibration Script
Calibrates Caglioti parameters from NIST standard sample data.

Usage:
    python calibrate_instrument.py --standard data/standards/LaB6.xy --output config.yaml
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import yaml
import numpy as np
import matplotlib.pyplot as plt

from src.preprocessing import load_xrd_data, apply_smoothing, subtract_background
from src.fitting import detect_peaks, fit_peaks_lm
from src.physics import CagliotiCorrection


def calibrate_instrument(standard_file: str, output_config: str = None):
    """
    Calibrate instrument using standard sample.
    
    Args:
        standard_file: Path to standard sample XRD data
        output_config: Path to save calibration (default: config.yaml)
    """
    print("=" * 60)
    print("AXCSAS Instrument Calibration")
    print("=" * 60)
    
    # Load data
    print(f"\n[1/4] Loading standard data: {standard_file}")
    two_theta, intensity = load_xrd_data(standard_file)
    print(f"      Loaded {len(two_theta)} data points")
    
    # Preprocess
    print("\n[2/4] Preprocessing...")
    intensity_smooth = apply_smoothing(intensity, window_size=11)
    intensity_bg = subtract_background(two_theta, intensity_smooth)
    
    # Detect and fit peaks
    print("\n[3/4] Fitting peaks...")
    peaks = detect_peaks(two_theta, intensity_bg, min_height=100)
    fit_results = fit_peaks_lm(two_theta, intensity_bg, peaks)
    
    peak_positions = np.array([r.center for r in fit_results])
    peak_fwhms = np.array([r.fwhm for r in fit_results])
    
    print(f"      Found {len(peak_positions)} peaks")
    for i, (pos, fwhm) in enumerate(zip(peak_positions, peak_fwhms)):
        print(f"        Peak {i+1}: 2θ = {pos:.2f}°, FWHM = {fwhm:.4f}°")
    
    # Calibrate Caglioti
    print("\n[4/4] Calibrating Caglioti parameters...")
    caglioti = CagliotiCorrection()
    params = caglioti.calibrate(peak_positions, peak_fwhms)
    
    print(f"\n      Calibration Results:")
    print(f"        U = {params.U:.6f} deg²")
    print(f"        V = {params.V:.6f} deg²")
    print(f"        W = {params.W:.6f} deg²")
    
    # Save to config
    if output_config:
        config_path = Path(output_config)
        
        if config_path.exists():
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
        else:
            config = {}
        
        if 'instrument' not in config:
            config['instrument'] = {}
        if 'caglioti' not in config['instrument']:
            config['instrument']['caglioti'] = {}
        
        config['instrument']['caglioti']['U'] = float(params.U)
        config['instrument']['caglioti']['V'] = float(params.V)
        config['instrument']['caglioti']['W'] = float(params.W)
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        
        print(f"\n✓ Saved calibration to {output_config}")
    
    print("\n" + "=" * 60)
    print("Calibration Complete!")
    print("=" * 60)
    
    return params


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate instrument using NIST standard sample"
    )
    parser.add_argument(
        "--standard", "-s",
        required=True,
        help="Path to standard sample XRD data (e.g., LaB6, Si)"
    )
    parser.add_argument(
        "--output", "-o",
        default="config.yaml",
        help="Output config file (default: config.yaml)"
    )
    
    args = parser.parse_args()
    
    calibrate_instrument(args.standard, args.output)


if __name__ == "__main__":
    main()
