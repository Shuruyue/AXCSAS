#!/usr/bin/env python3
"""
Advanced tests for copper crystal elastic moduli calculations.
測試銅晶體彈性模量進階計算。
"""

import pytest
from axcsas.core.copper_crystal import (
    calculate_youngs_modulus_from_stiffness,
    calculate_d_spacing,
    CU_ELASTIC
)


class TestElasticModuliCalculation:
    """Test direction-dependent elastic moduli calculations."""
    
    def test_youngs_modulus_100_direction(self):
        """Test E_100 calculation matches verified value."""
        E_100 = calculate_youngs_modulus_from_stiffness(1, 0, 0)
        
        # Verified from scripts/verify_elastic_moduli.py
        assert abs(E_100 - 66.7) < 0.5, f"E_100 should be ~66.7 GPa, got {E_100:.1f}"
        
        # Should match dataclass value
        assert abs(E_100 - CU_ELASTIC.E_100) < 0.5
    
    def test_youngs_modulus_111_direction(self):
        """Test E_111 calculation matches verified value."""
        E_111 = calculate_youngs_modulus_from_stiffness(1, 1, 1)
        
        # Verified from scripts/verify_elastic_moduli.py
        assert abs(E_111 - 191.1) < 0.5, f"E_111 should be ~191.1 GPa, got {E_111:.1f}"
        
        # Should match dataclass value
        assert abs(E_111 - CU_ELASTIC.E_111) < 0.5
    
    def test_youngs_modulus_110_direction(self):
        """Test E_110 calculation matches verified value."""
        E_110 = calculate_youngs_modulus_from_stiffness(1, 1, 0)
        
        # Verified from scripts/verify_elastic_moduli.py
        assert abs(E_110 - 130.3) < 0.5, f"E_110 should be ~130.3 GPa, got {E_110:.1f}"
        
        # Should match dataclass value
        assert abs(E_110 - CU_ELASTIC.E_110) < 0.5
    
    def test_elastic_anisotropy(self):
        """Test that copper shows expected elastic anisotropy."""
        E_100 = calculate_youngs_modulus_from_stiffness(1, 0, 0)
        E_111 = calculate_youngs_modulus_from_stiffness(1, 1, 1)
        
        # [111] should be hardest, [100] softest
        assert E_111 > E_100, "[111] should be stiffer than [100]"
        
        # Ratio should be ~2.86 (191.1 / 66.7)
        ratio = E_111 / E_100
        assert 2.5 < ratio < 3.2, f"Anisotropy ratio should be ~2.86, got {ratio:.2f}"
    
    def test_custom_elastic_constants(self):
        """Test calculation with custom elastic constants."""
        # Use hypothetical values
        E = calculate_youngs_modulus_from_stiffness(
            1, 0, 0,
            C11=200.0,
            C12=100.0,
            C44=80.0
        )
        
        # Should return a positive value
        assert E > 0, "Modulus should be positive"
        
        # Should be physically reasonable (0-300 GPa for metals)
        assert E < 300, "Modulus should be < 300 GPa"
    
    def test_equivalent_directions(self):
        """Test that equivalent crystallographic directions give same modulus."""
        # All <100> directions should be equivalent
        E_100 = calculate_youngs_modulus_from_stiffness(1, 0, 0)
        E_010 = calculate_youngs_modulus_from_stiffness(0, 1, 0)
        E_001 = calculate_youngs_modulus_from_stiffness(0, 0, 1)
        
        assert abs(E_100 - E_010) < 0.1, "<100> and <010> should be equivalent"
        assert abs(E_100 - E_001) < 0.1, "<100> and <001> should be equivalent"
        
        # All <111> directions should be equivalent
        E_111 = calculate_youngs_modulus_from_stiffness(1, 1, 1)
        E_111_neg = calculate_youngs_modulus_from_stiffness(-1, 1, 1)
        
        assert abs(E_111 - E_111_neg) < 0.1, "<111> and <-111> should be equivalent"


class TestDSpacingCalculation:
    """Test d-spacing calculation function."""
    
    def test_d_spacing_111(self):
        """Test d_111 calculation."""
        d_111 = calculate_d_spacing(1, 1, 1)
        
        # d_111 = a / √3 = 3.6150 / 1.732 = 2.087 Å
        expected = 3.6150 / (3 ** 0.5)
        assert abs(d_111 - expected) < 0.001
    
    def test_d_spacing_200(self):
        """Test d_200 calculation."""
        d_200 = calculate_d_spacing(2, 0, 0)
        
        # d_200 = a / 2 = 3.6150 / 2 = 1.8075 Å
        expected = 3.6150 / 2
        assert abs(d_200 - expected) < 0.001
    
    def test_custom_lattice_constant(self):
        """Test d-spacing with custom lattice constant."""
        a_custom = 4.0  # Å
        d = calculate_d_spacing(1, 1, 1, a=a_custom)
        
        expected = a_custom / (3 ** 0.5)
        assert abs(d - expected) < 0.001


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
