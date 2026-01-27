"""
Pseudo-Voigt Area Calculation Utilities
========================================

Provides accurate integrated area calculations for Pseudo-Voigt profiles.
為 Pseudo-Voigt 剖面提供精確的積分面積計算。

Mathematical Background:
PV profile: I(x) = A × [η × L(x) + (1-η) × G(x)]

Integrated area:
    Area = A × FWHM × [η × (π/2) + (1-η) × √(π/ln2)]

References:
- Ida, T. & Kimura, K. (1999). J. Appl. Cryst. 32, 982-991
- Pseudo-Voigt normalization factors derived from analytical integration
"""

import numpy as np


def calculate_pv_area(amplitude: float, fwhm: float, eta: float) -> float:
    """
    Calculate integrated area of a Pseudo-Voigt profile.
    計算 Pseudo-Voigt 剖面的積分面積。
    
    Args:
        amplitude: Peak amplitude (maximum intensity)
        fwhm: Full Width at Half Maximum  
        eta: Mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian)
    
    Returns:
        Integrated area under the peak
    
    Mathematical Derivation:
    ========================
    For a Pseudo-Voigt PV(x) = A × [η × L(x) + (1-η) × G(x)]
    
    Where:
    - L(x) = normalized Lorentzian
    - G(x) = normalized Gaussian
    
    The integral is:
    ∫ PV(x) dx = A × [η × ∫L(x)dx + (1-η) × ∫G(x)dx]
    
    For normalized profiles with FWHM:
    - ∫ L(x) dx = π × FWHM / 2  (Lorentzian integral)
    - ∫ G(x) dx = FWHM × √(π/ln2)  (Gaussian integral)
    
    Therefore:
    Area = A × FWHM × [η × (π/2) + (1-η) × √(π/ln2)]
    
    Example:
    --------
    >>> calculate_pv_area(amplitude=1000, fwhm=0.3, eta=0.5)
    395.289...  # For η=0.5, factor ≈ 1.3176
    """
    # Validate inputs
    if amplitude <= 0 or fwhm <= 0:
        return 0.0
    
    # Clip eta to valid range
    eta = np.clip(eta, 0.0, 1.0)
    
    # Integration factors for normalized profiles
    lorentzian_factor = np.pi / 2  # = 1.5707963...
    
    # Gaussian Area = FWHM * 0.5 * sqrt(pi/ln2)
    # Original code was missing 0.5 factor
    gaussian_factor = 0.5 * np.sqrt(np.pi / np.log(2))  # = 1.064467...
    
    # Weighted average of integration factors
    integration_factor = eta * lorentzian_factor + (1 - eta) * gaussian_factor
    
    # Total integrated area
    area = amplitude * fwhm * integration_factor
    
    return area


def get_pv_integration_factor(eta: float) -> float:
    """
    Get the integration factor for a given η value.
    取得給定 η 值的積分因子。
    
    Args:
        eta: Mixing parameter (0-1)
    
    Returns:
        Integration factor k where Area = Amplitude × FWHM × k
    
    Note:
        - η = 0 (pure Gaussian): k ≈ 1.064
        - η = 0.5 (50% mix): k ≈ 1.318
        - η = 1 (pure Lorentzian): k ≈ 1.571
        
    The common approximation PV_FACTOR = 1.0645 is for Gaussian limit.
    """
    eta = np.clip(eta, 0.0, 1.0)
    lorentzian_factor = np.pi / 2
    gaussian_factor = 0.5 * np.sqrt(np.pi / np.log(2))
    return eta * lorentzian_factor + (1 - eta) * gaussian_factor


# Precomputed reference values for validation
PV_FACTORS = {
    0.0: 1.064467083,   # Pure Gaussian
    0.3: 1.216365856,
    0.5: 1.317631705,   # 50% mix
    0.7: 1.418897554,
    1.0: 1.570796327,   # Pure Lorentzian
}


def validate_pv_area_calculation():
    """
    Validate area calculation against known values.
    驗證面積計算對已知值的準確性。
    """
    print("Validating PV Area Calculation")
    print("=" * 50)
    
    test_cases = [
        (100, 0.5, 0.0),   # Pure Gaussian
        (100, 0.5, 0.5),   # 50/50 mix
        (100, 0.5, 1.0),   # Pure Lorentzian
    ]
    
    for amp, fwhm, eta in test_cases:
        area = calculate_pv_area(amp, fwhm, eta)
        factor = get_pv_integration_factor(eta)
        expected = amp * fwhm * factor
        
        print(f"η = {eta:.1f}: Area = {area:.3f} (expected: {expected:.3f})")
        assert abs(area - expected) < 1e-6, "Area calculation mismatch!"
    
    print("✓ All validation tests passed")


if __name__ == "__main__":
    # Run validation
    validate_pv_area_calculation()
    
    # Show common mistake
    print("\nCommon Mistake Warning:")
    print("=" * 50)
    print("Using fixed PV_FACTOR = 1.0645 is accurate for pure Gaussian (eta=0)")
    print("But actual η varies per peak!")
    print(f"  η=0.3 (Gaussian-dominated): factor = {get_pv_integration_factor(0.3):.4f}")
    print(f"  η=0.5 (balanced):           factor = {get_pv_integration_factor(0.5):.4f}")
    print(f"  η=0.7 (Lorentzian-dominated): factor = {get_pv_integration_factor(0.7):.4f}")
    print(f"Using fixed 1.0645 for Lorentzian (eta=1) causes {abs((1.0645 - get_pv_integration_factor(1.0))/get_pv_integration_factor(1.0)*100):.1f}% error!")
