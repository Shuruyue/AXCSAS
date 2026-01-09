#!/usr/bin/env python3
"""
Test AXCSAS Pipeline with Real Data
====================================

Integration test using actual XRD data files.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from integration.pipeline import (
    AXCSASPipeline,
    AnalysisConfig,
    run_full_analysis,
)


def main():
    """Test pipeline with real data."""
    
    # Find data files
    data_dir = Path(__file__).parent / 'data' / 'raw' / '202511'
    
    if not data_dir.exists():
        print(f"Data directory not found: {data_dir}")
        return
    
    # Get first file of each concentration
    test_files = [
        data_dir / '20251125_0ml_0h_15min.txt',
        data_dir / '20251125_4.5ml_0h_7min.txt',
        data_dir / '20251125_9ml_0h_5min.txt',
        data_dir / '20251125_18ml_0h_10min.txt',
    ]
    
    # Run pipeline
    pipeline = AXCSASPipeline()
    
    print("=" * 70)
    print("AXCSAS Integration Test - Real Data Analysis")
    print("=" * 70)
    
    for filepath in test_files:
        if not filepath.exists():
            print(f"File not found: {filepath.name}")
            continue
        
        print(f"\nAnalyzing: {filepath.name}")
        print("-" * 50)
        
        # Sample age = 0 for as-deposited
        result = pipeline.analyze(str(filepath), sample_age_hours=0.5)
        
        # Show summary
        print(f"Peaks found: {len(result.peaks)}")
        
        if result.average_size_nm:
            print(f"Scherrer Size: {result.average_size_nm:.1f} nm")
        
        if result.wh_result:
            print(f"W-H Size: {result.wh_result.crystallite_size_nm:.1f} nm")
            print(f"W-H Strain: {result.wh_result.microstrain:.2e}")
            print(f"W-H R²: {result.wh_result.r_squared:.3f}")
        
        if result.texture_result:
            dom = result.texture_result.dominant_hkl
            if dom:
                print(f"Dominant: ({dom[0]}{dom[1]}{dom[2]}) TC={result.texture_result.dominant_tc:.2f}")
            print(f"Random: {result.texture_result.is_random}")
        
        if result.stacking_fault:
            print(f"Peak Sep: {result.stacking_fault.peak_separation_deg:.3f}°")
            print(f"SF α: {result.stacking_fault.alpha_percent:.2f}%")
        
        if result.lattice_result:
            print(f"Lattice a: {result.lattice_result.lattice_constant:.4f} Å")
    
    print("\n" + "=" * 70)
    print("Analysis Complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
