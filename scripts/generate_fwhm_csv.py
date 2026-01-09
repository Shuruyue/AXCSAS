#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FWHM Calculator for Cu XRD Peaks

Calculates Full Width at Half Maximum (FWHM) for Cu XRD peaks:
- Cu (111): ~43.3° 2θ
- Cu (200): ~50.4° 2θ  
- Cu (220): ~74.1° 2θ

Output format:
- One CSV per concentration
- Columns: Time points (0h, 2h, ..., 24h)
- Rows: Peak directions (111, 200, 220)
"""

import re
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.interpolate import interp1d


# Cu peak positions (2θ in degrees, Cu Kα radiation)
CU_PEAKS = {
    '111': 43.3,
    '200': 50.4,
    '220': 74.1
}

# Search range around each peak (±degrees)
PEAK_SEARCH_RANGE = 2.0


def parse_xrd_file(filepath: Path) -> tuple:
    """Parse Bruker XRD TXT file and extract angle-intensity data."""
    angles = []
    intensities = []
    in_data_section = False
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            
            if line == '[Data]':
                in_data_section = True
                continue
            
            if in_data_section:
                if 'Angle' in line:
                    continue
                
                if line and (line[0].isdigit() or line[0] == ' '):
                    parts = line.split(',')
                    if len(parts) >= 2:
                        try:
                            angle = float(parts[0].strip())
                            intensity = int(parts[1].strip())
                            angles.append(angle)
                            intensities.append(intensity)
                        except (ValueError, IndexError):
                            continue
    
    return np.array(angles), np.array(intensities)


def calculate_fwhm(angles: np.ndarray, intensities: np.ndarray, 
                   peak_position: float, search_range: float = 2.0) -> float:
    """
    Calculate FWHM for a peak near the specified position.
    
    Args:
        angles: Array of 2θ angles
        intensities: Array of intensity values
        peak_position: Expected peak position in degrees
        search_range: Range to search around peak position (±degrees)
    
    Returns:
        FWHM in degrees, or NaN if calculation fails
    """
    # Extract region around expected peak
    mask = (angles >= peak_position - search_range) & (angles <= peak_position + search_range)
    region_angles = angles[mask]
    region_intensities = intensities[mask]
    
    if len(region_angles) < 5:
        return np.nan
    
    # Find the actual peak maximum in this region
    peak_idx = np.argmax(region_intensities)
    peak_angle = region_angles[peak_idx]
    peak_intensity = region_intensities[peak_idx]
    
    # Estimate background (minimum in the region)
    background = np.min(region_intensities)
    
    # Half maximum level
    half_max = (peak_intensity + background) / 2
    
    # Find where the curve crosses half maximum
    above_half = region_intensities > half_max
    
    # Find left crossing point
    left_idx = peak_idx
    while left_idx > 0 and above_half[left_idx]:
        left_idx -= 1
    
    # Find right crossing point  
    right_idx = peak_idx
    while right_idx < len(above_half) - 1 and above_half[right_idx]:
        right_idx += 1
    
    # Interpolate to find exact crossing points
    try:
        # Left side interpolation
        if left_idx < peak_idx and region_intensities[left_idx] != region_intensities[left_idx + 1]:
            left_angle = region_angles[left_idx] + (half_max - region_intensities[left_idx]) * \
                         (region_angles[left_idx + 1] - region_angles[left_idx]) / \
                         (region_intensities[left_idx + 1] - region_intensities[left_idx])
        else:
            left_angle = region_angles[left_idx]
        
        # Right side interpolation
        if right_idx > peak_idx and region_intensities[right_idx] != region_intensities[right_idx - 1]:
            right_angle = region_angles[right_idx - 1] + (half_max - region_intensities[right_idx - 1]) * \
                          (region_angles[right_idx] - region_angles[right_idx - 1]) / \
                          (region_intensities[right_idx] - region_intensities[right_idx - 1])
        else:
            right_angle = region_angles[right_idx]
        
        fwhm = right_angle - left_angle
        
        # Sanity check
        if fwhm <= 0 or fwhm > 2 * search_range:
            return np.nan
        
        return fwhm
        
    except Exception:
        return np.nan


def parse_filename(filename: str) -> tuple:
    """Parse filename to extract concentration and time."""
    match = re.match(r'\d+_(\d+\.?\d*ml)_(\d+h)(?:_\d+min)?\.txt', filename, re.IGNORECASE)
    if not match:
        return None, None
    
    concentration = match.group(1)
    time_str = match.group(2)
    time_hours = float(time_str.replace('h', ''))
    
    return concentration, time_hours


def generate_fwhm_csv(files: list, concentration: str, output_dir: Path) -> Path:
    """Generate FWHM CSV for a specific concentration."""
    # Sort files by time
    files_sorted = sorted(files, key=lambda x: x[1])
    
    # Time points as columns
    time_points = [f"{int(t)}h" if t == int(t) else f"{t}h" for _, t in files_sorted]
    
    # Calculate FWHM for each file and each peak
    data = {'Peak': ['111', '200', '220']}
    
    for filepath, time_hours in files_sorted:
        angles, intensities = parse_xrd_file(filepath)
        time_label = f"{int(time_hours)}h" if time_hours == int(time_hours) else f"{time_hours}h"
        
        fwhm_values = []
        for peak_name, peak_pos in CU_PEAKS.items():
            fwhm = calculate_fwhm(angles, intensities, peak_pos, PEAK_SEARCH_RANGE)
            fwhm_values.append(round(fwhm, 4) if not np.isnan(fwhm) else np.nan)
        
        data[time_label] = fwhm_values
    
    df = pd.DataFrame(data)
    output_path = output_dir / f"fwhm_{concentration}.csv"
    df.to_csv(output_path, index=False)
    
    print(f"Generated: {output_path.name}")
    print(f"  - FWHM values (° 2θ): {df.to_dict('list')}")
    
    return output_path


def main():
    # Setup paths
    script_path = Path(__file__).resolve()
    project_dir = script_path.parent.parent
    data_dir = project_dir / "data" / "raw" / "202511"
    output_dir = project_dir / "outputs"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("FWHM Calculator for Cu XRD Peaks")
    print("=" * 60)
    print(f"Data source: {data_dir}")
    print(f"Cu peaks: 111 (43.3°), 200 (50.4°), 220 (74.1°)")
    print()
    
    # Group files by concentration
    files_by_concentration = defaultdict(list)
    
    for txt_file in data_dir.glob("*.txt"):
        concentration, time_hours = parse_filename(txt_file.name)
        if concentration and time_hours is not None:
            files_by_concentration[concentration].append((txt_file, time_hours))
    
    concentrations = ['0ml', '4.5ml', '9ml', '18ml']
    
    print(f"Found {len(files_by_concentration)} concentrations")
    print()
    
    # Generate FWHM CSV for each concentration
    generated_files = []
    for conc in concentrations:
        if conc in files_by_concentration:
            output_path = generate_fwhm_csv(files_by_concentration[conc], conc, output_dir)
            generated_files.append(output_path)
            print()
    
    print("=" * 60)
    print(f"Generated {len(generated_files)} FWHM CSV files")
    print("=" * 60)


if __name__ == "__main__":
    main()
