#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
XRD Data CSV Generator for AXCSAS Project

Generates CSV files for XRD data organized by concentration.
- Columns (X-axis): Time points (0h, 2h, 4h, ..., 24h)
- Rows (Y-axis): 2θ angles (40°~78°)
- Values: PSD intensity signals

Output: 4 CSV files for 0ml, 4.5ml, 9ml, 18ml concentrations
"""

import os
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd


def parse_xrd_file(filepath: Path) -> tuple:
    """
    Parse Bruker XRD TXT file and extract angle-intensity data.
    
    Returns:
        tuple: (angles_list, intensities_list)
    """
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
                # Skip header line
                if 'Angle' in line:
                    continue
                
                # Parse data line: "   40.0015,       645,"
                if line and line[0].isdigit() or (line and line[0] == ' '):
                    parts = line.split(',')
                    if len(parts) >= 2:
                        try:
                            angle = float(parts[0].strip())
                            intensity = int(parts[1].strip())
                            angles.append(angle)
                            intensities.append(intensity)
                        except (ValueError, IndexError):
                            continue
    
    return angles, intensities


def parse_filename(filename: str) -> tuple:
    """
    Parse filename to extract concentration and time.
    
    Example: "20251125_0ml_0h_15min.txt" -> ("0ml", 0.0)
             "20251126_4.5ml_6h.txt" -> ("4.5ml", 6.0)
    """
    # Pattern: DATE_CONC_TIME.txt
    match = re.match(r'\d+_(\d+\.?\d*ml)_(\d+h)(?:_\d+min)?\.txt', filename, re.IGNORECASE)
    if not match:
        return None, None
    
    concentration = match.group(1)
    time_str = match.group(2)
    time_hours = float(time_str.replace('h', ''))
    
    return concentration, time_hours


def generate_csv_for_concentration(
    files: list,
    concentration: str,
    output_dir: Path
) -> Path:
    """
    Generate CSV file for a specific concentration.
    
    Args:
        files: List of (filepath, time_hours) tuples
        concentration: Concentration string (e.g., "0ml")
        output_dir: Output directory path
    
    Returns:
        Path to generated CSV file
    """
    # Sort files by time
    files_sorted = sorted(files, key=lambda x: x[1])
    
    # Parse first file to get angles (should be same for all)
    angles, _ = parse_xrd_file(files_sorted[0][0])
    
    # Create DataFrame with angles as index
    df = pd.DataFrame({'2Theta': angles})
    
    # Add intensity columns for each time point
    for filepath, time_hours in files_sorted:
        _, intensities = parse_xrd_file(filepath)
        
        # Format column name
        time_label = f"{int(time_hours)}h" if time_hours == int(time_hours) else f"{time_hours}h"
        df[time_label] = intensities
    
    # Save to CSV
    output_path = output_dir / f"xrd_{concentration}.csv"
    df.to_csv(output_path, index=False)
    
    print(f"Generated: {output_path.name}")
    print(f"  - Rows (angles): {len(df)}")
    print(f"  - Columns: {list(df.columns)}")
    
    return output_path


def main():
    # Setup paths - use __file__ for reliable path resolution
    script_path = Path(__file__).resolve()
    scripts_dir = script_path.parent
    project_dir = scripts_dir.parent
    data_dir = project_dir / "data" / "raw" / "202511"
    output_dir = project_dir / "outputs"
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("XRD Data CSV Generator")
    print("=" * 60)
    print(f"Script location: {script_path}")
    print(f"Data source: {data_dir}")
    print(f"Data exists: {data_dir.exists()}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Group files by concentration
    files_by_concentration = defaultdict(list)
    
    for txt_file in data_dir.glob("*.txt"):
        concentration, time_hours = parse_filename(txt_file.name)
        if concentration and time_hours is not None:
            files_by_concentration[concentration].append((txt_file, time_hours))
    
    # Sort concentrations for consistent ordering
    concentrations = ['0ml', '4.5ml', '9ml', '18ml']
    
    print(f"Found {len(files_by_concentration)} concentrations:")
    for conc in concentrations:
        if conc in files_by_concentration:
            print(f"  - {conc}: {len(files_by_concentration[conc])} files")
    print()
    
    # Generate CSV for each concentration
    generated_files = []
    for conc in concentrations:
        if conc in files_by_concentration:
            output_path = generate_csv_for_concentration(
                files_by_concentration[conc],
                conc,
                output_dir
            )
            generated_files.append(output_path)
            print()
    
    print("=" * 60)
    print(f"Successfully generated {len(generated_files)} CSV files:")
    for f in generated_files:
        print(f"  - {f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
