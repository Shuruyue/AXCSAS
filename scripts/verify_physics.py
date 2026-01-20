#!/usr/bin/env python3
"""
AXCSAS Physics Verification Script
AXCSAS 物理驗證腳本
==================================

Calculates theoretical 2-theta positions for Copper FCC crystal
based on fundamental constants to ensure absolute precision.
根據基本常數計算銅 FCC 晶體的理論 2-theta 位置，確保絕對精確。

Constants / 常數 (No rounding / 不四捨五入):
- Wavelength (Cu Kα1): 1.540562 Å (Bearden 1967, Rev. Mod. Phys. 39, 78)
- Lattice Constant (a₀): 3.6150 Å (JCPDS 04-0836)
"""

import math

def calculate_theoretical_values():
    print("=" * 60)
    print("AXCSAS Physics Verification: Copper FCC Peak Positions")
    print("=" * 60)

    # Fundamental Constants (Exact values from literature)
    # Bearden (1967) Table V
    LAMBDA_KA1 = 1.540562  # Å
    # JCPDS 04-0836
    LATTICE_A = 3.6150     # Å
    
    print(f"Constants:")
    print(f"  λ (Cu Kα1) = {LAMBDA_KA1:.6f} Å (Bearden 1967)")
    print(f"  a₀ (Cu)    = {LATTICE_A:.4f} Å   (JCPDS 04-0836)")
    print("-" * 60)
    print(f"{'hkl':<10} {'d (Å)':<15} {'sin(θ)':<15} {'2θ (deg)':<15}")
    print("-" * 60)

    # hkl planes to verify
    planes = [
        (1, 1, 1),
        (2, 0, 0),
        (2, 2, 0),
        (3, 1, 1),
        (2, 2, 2)
    ]

    results = {}

    for h, k, l in planes:
        # 1. Calculate d-spacing
        # d = a / sqrt(h² + k² + l²)
        hkl_sq = h**2 + k**2 + l**2
        d_spacing = LATTICE_A / math.sqrt(hkl_sq)

        # 2. Calculate theta using Bragg's Law
        # 2d sin(θ) = λ  =>  sin(θ) = λ / 2d
        sin_theta = LAMBDA_KA1 / (2 * d_spacing)
        
        # Check physical limit
        if sin_theta > 1.0:
            print(f"({h}{k}{l}) matches unreachable angle (sin(θ) > 1)")
            continue

        theta_rad = math.asin(sin_theta)
        two_theta_deg = 2 * math.degrees(theta_rad)

        results[(h, k, l)] = two_theta_deg

        print(f"({h}{k}{l}):    {d_spacing:<15.5f} {sin_theta:<15.5f} {two_theta_deg:<15.5f}")

    print("-" * 60)
    print("Recommended Values for copper_crystal.py (to 3 decimal places):")
    for hkl, angle in results.items():
        print(f"    {hkl}: {angle:.3f},")
    
    # Return exit code for CI integration
    return 0  # All calculations successful

if __name__ == "__main__":
    import sys
    sys.exit(calculate_theoretical_values())
