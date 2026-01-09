"""
Unit Tests for Copper Crystal Module
=====================================

Tests physical constants, JCPDS data, Scherrer K values,
and lattice constant validation functions.

Run with: pytest tests/test_copper_crystal.py -v
"""

import pytest
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from core.copper_crystal import (
    CopperCrystal,
    CU_CRYSTAL,
    CU_JCPDS_EXTENDED,
    ScherrerCubicK,
    SCHERRER_CUBIC_K,
    CopperElasticModuli,
    CU_ELASTIC,
    get_k_for_hkl,
    get_youngs_modulus,
    is_fcc_allowed,
    validate_lattice_constant,
    explain_lattice_deviation,
    calculate_d_spacing,
    ELECTROPLATED_A_STANDARD,
)

from core.additives import (
    AdditiveType,
    Additive,
    SUPPRESSOR,
    ACCELERATOR,
    LEVELER_JGB,
    get_additive,
)


class TestCopperCrystalConstants:
    """Tests for FCC copper structure constants."""
    
    def test_lattice_constant_standard_value(self):
        """TC-001: Verify standard lattice constant is 3.6150 Å."""
        assert CU_CRYSTAL.lattice_constant == 3.6150
        
    def test_space_group(self):
        """Verify FCC space group."""
        assert CU_CRYSTAL.space_group == "Fm-3m"
        assert CU_CRYSTAL.space_group_number == 225
        
    def test_density(self):
        """Verify copper density."""
        assert CU_CRYSTAL.density == 8.92
        
    def test_packing_factor(self):
        """Verify FCC atomic packing factor."""
        assert CU_CRYSTAL.packing_factor == 0.74
        
    def test_crystal_is_immutable(self):
        """Verify crystal constants are frozen."""
        with pytest.raises(Exception):  # FrozenInstanceError
            CU_CRYSTAL.lattice_constant = 3.62


class TestJCPDSData:
    """Tests for JCPDS 04-0836 standard data."""
    
    def test_jcpds_111_two_theta(self):
        """TC-002: Verify (111) 2θ position."""
        assert CU_JCPDS_EXTENDED[(1, 1, 1)]["two_theta"] == 43.297
        
    def test_jcpds_completeness(self):
        """TC-003: Verify all 5 standard peaks are present."""
        expected_peaks = [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (2, 2, 2)]
        assert len(CU_JCPDS_EXTENDED) == 5
        for peak in expected_peaks:
            assert peak in CU_JCPDS_EXTENDED
            
    def test_jcpds_has_multiplicity(self):
        """Verify multiplicity field is present for all peaks."""
        for hkl, data in CU_JCPDS_EXTENDED.items():
            assert "multiplicity" in data
            assert data["multiplicity"] > 0
            
    def test_jcpds_111_is_strongest(self):
        """Verify (111) has highest intensity."""
        assert CU_JCPDS_EXTENDED[(1, 1, 1)]["intensity"] == 100
        
    def test_jcpds_multiplicities_correct(self):
        """Verify multiplicity values are physically correct."""
        assert CU_JCPDS_EXTENDED[(1, 1, 1)]["multiplicity"] == 8   # {111}
        assert CU_JCPDS_EXTENDED[(2, 0, 0)]["multiplicity"] == 6   # {200}
        assert CU_JCPDS_EXTENDED[(2, 2, 0)]["multiplicity"] == 12  # {220}


class TestFCCExtinctionRule:
    """Tests for FCC systematic absence rules."""
    
    def test_fcc_allowed_all_odd(self):
        """TC-004a: (111) should be allowed (all odd)."""
        assert is_fcc_allowed(1, 1, 1) is True
        
    def test_fcc_allowed_all_even(self):
        """(200) should be allowed (all even)."""
        assert is_fcc_allowed(2, 0, 0) is True
        assert is_fcc_allowed(2, 2, 0) is True
        
    def test_fcc_forbidden_mixed(self):
        """TC-004b: (100) should be forbidden (mixed)."""
        assert is_fcc_allowed(1, 0, 0) is False
        assert is_fcc_allowed(1, 1, 0) is False
        assert is_fcc_allowed(2, 1, 0) is False


class TestScherrerCubicK:
    """Tests for direction-dependent Scherrer K values."""
    
    def test_k_111_cubic_habit(self):
        """TC-005: K(111) should be 1.155 for cubic habit."""
        k = get_k_for_hkl(1, 1, 1, use_cubic_habit=True)
        assert abs(k - 1.155) < 0.001
        
    def test_k_200_cubic_habit(self):
        """TC-006: K(200) should be 1.000 for cubic habit."""
        k = get_k_for_hkl(2, 0, 0, use_cubic_habit=True)
        assert abs(k - 1.000) < 0.001
        
    def test_k_220_cubic_habit(self):
        """K(220) should be 1.061 for cubic habit (文件 04 §2.2)."""
        k = get_k_for_hkl(2, 2, 0, use_cubic_habit=True)
        assert abs(k - 1.061) < 0.001
        
    def test_k_spherical_fallback(self):
        """Spherical assumption should return 0.89."""
        k = get_k_for_hkl(1, 1, 1, use_cubic_habit=False)
        assert k == 0.89
        
    def test_k_constants_dataclass(self):
        """Verify ScherrerCubicK dataclass values."""
        assert SCHERRER_CUBIC_K.K_111 == 1.155
        assert SCHERRER_CUBIC_K.K_200 == 1.000
        assert SCHERRER_CUBIC_K.K_SPHERICAL == 0.89


class TestElasticAnisotropy:
    """Tests for direction-dependent elastic moduli."""
    
    def test_youngs_modulus_111(self):
        """TC-007: E<111> should be 191 GPa."""
        E = get_youngs_modulus(1, 1, 1)
        assert E == 191.0
        
    def test_youngs_modulus_200(self):
        """E<100> should be 66 GPa (softest)."""
        E = get_youngs_modulus(2, 0, 0)
        assert E == 66.0
        
    def test_youngs_modulus_220(self):
        """E<110> should be 130 GPa."""
        E = get_youngs_modulus(2, 2, 0)
        assert E == 130.0
        
    def test_anisotropy_ratio(self):
        """Anisotropy ratio E_111/E_100 should be ~3."""
        ratio = CU_ELASTIC.E_111 / CU_ELASTIC.E_100
        assert 2.8 < ratio < 3.0


class TestLatticeValidation:
    """Tests for lattice constant validation and explanation."""
    
    def test_standard_lattice_is_valid(self):
        """Standard value should pass validation."""
        result = validate_lattice_constant(3.6150)
        assert result.is_normal is True
        assert result.warning_level == "none"
        
    def test_slight_expansion_detected(self):
        """Slight expansion should be flagged as minor."""
        result = validate_lattice_constant(3.6180)
        assert result.warning_level == "minor"
        
    def test_significant_expansion_detected(self):
        """Large expansion should be flagged as significant."""
        result = validate_lattice_constant(3.6200)
        assert result.is_normal is False
        assert result.warning_level == "significant"
        
    def test_contraction_flagged(self):
        """Contraction should be flagged."""
        result = validate_lattice_constant(3.6100)
        assert result.warning_level == "significant"
        
    def test_explanation_includes_causes(self):
        """Explanation should include physical causes."""
        explanation = explain_lattice_deviation(3.6180)
        assert "Sulfur" in explanation or "impurity" in explanation.lower()
        
    def test_explanation_with_sample_age(self):
        """Explanation should include self-annealing context."""
        explanation = explain_lattice_deviation(3.6180, sample_age_hours=0.5)
        assert "as-deposited" in explanation


class TestAdditives:
    """Tests for additive database."""
    
    def test_additive_types_exist(self):
        """TC-008: Verify 3 additive types exist."""
        assert len(AdditiveType) == 3
        
    def test_suppressor_properties(self):
        """Verify suppressor additive properties."""
        assert SUPPRESSOR.type == AdditiveType.SUPPRESSOR
        assert "PEG" in SUPPRESSOR.common_chemicals
        
    def test_get_additive_by_name(self):
        """Test additive lookup."""
        additive = get_additive("jgb")
        assert additive is not None
        assert additive.type == AdditiveType.LEVELER


class TestCalculations:
    """Tests for crystallographic calculations."""
    
    def test_d_spacing_111(self):
        """Verify d-spacing calculation for (111)."""
        d = calculate_d_spacing(1, 1, 1, a=3.6150)
        assert abs(d - 2.088) < 0.01
        
    def test_d_spacing_200(self):
        """Verify d-spacing calculation for (200)."""
        d = calculate_d_spacing(2, 0, 0, a=3.6150)
        assert abs(d - 1.8075) < 0.01


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
