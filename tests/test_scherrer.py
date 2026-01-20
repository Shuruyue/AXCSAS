"""
Unit Tests for Enhanced Scherrer Calculator
============================================

Tests K value lookup, unit conversion, validity flags, and calculation accuracy.

Run with: pytest tests/test_scherrer.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path


from axcsas.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
    calculate_scherrer,
    generate_scherrer_report,
    FWHM_RATIO_THRESHOLD,
)
from axcsas.core.copper_crystal import (
    get_k_for_hkl,
    SCHERRER_CUBIC_K,
)


class TestKValueLookup:
    """Tests for dynamic K value selection (K_w for FWHM)."""
    
    def test_k_111_value(self):
        """K(111) should be 0.855 (K_w from L&W 1978 Table 2)."""
        k = get_k_for_hkl(1, 1, 1)
        assert abs(k - 0.855) < 0.001
    
    def test_k_200_value(self):
        """K(200) should be 0.886 (K_w from L&W 1978 Table 2)."""
        k = get_k_for_hkl(2, 0, 0)
        assert abs(k - 0.886) < 0.001
    
    def test_k_220_value(self):
        """K(220) should be 0.834 (K_w from L&W 1978 Table 2, via 110)."""
        k = get_k_for_hkl(2, 2, 0)
        assert abs(k - 0.834) < 0.001
    
    def test_k_311_value(self):
        """K(311) should be 0.908 (K_w from L&W 1978 Table 2)."""
        k = get_k_for_hkl(3, 1, 1)
        assert abs(k - 0.908) < 0.001
    
    def test_k_spherical_fallback(self):
        """Non-cubic habit should use K=0.829 (L&W 1978 spherical)."""
        k = get_k_for_hkl(1, 1, 1, use_cubic_habit=False)
        assert abs(k - 0.829) < 0.001
    
    def test_scherrer_cubic_k_class(self):
        """ScherrerCubicK class should have K_w values."""
        assert abs(SCHERRER_CUBIC_K.K_111 - 0.855) < 0.001
        assert abs(SCHERRER_CUBIC_K.K_200 - 0.886) < 0.001
        assert abs(SCHERRER_CUBIC_K.K_220 - 0.834) < 0.001
        assert abs(SCHERRER_CUBIC_K.K_311 - 0.908) < 0.001


class TestUnitConversion:
    """Tests for critical unit conversion."""
    
    def test_degrees_to_radians(self):
        """Verify degree to radian conversion."""
        # 0.25° should be ~0.00436 rad
        deg = 0.25
        rad = deg * np.pi / 180
        assert abs(rad - 0.00436) < 0.0001
    
    def test_angstrom_to_nm(self):
        """Verify Å to nm conversion."""
        angstrom = 490
        nm = angstrom / 10
        assert nm == 49.0


class TestDocumentExample:
    """Test case from document 04 §4 (updated for K_w)."""
    
    def test_document_example_calculation(self):
        """
        Verify calculation using K_w (FWHM) values.
        
        Input:
            2θ = 43.32°
            FWHM_obs = 0.25°
            FWHM_inst = 0.08°
            K = 0.855 (K_w for cubic (111) from L&W 1978)
            
        Calculation:
            β_sample = √(0.25² - 0.08²) = 0.237° = 0.00414 rad
            D = Kλ / (β cosθ) = 0.855 × 1.540562 / (0.00414 × cos(21.66°))
            D = 1.317 / 0.00385 = 342 Å = 34.2 nm
        """
        result = calculate_scherrer(
            two_theta=43.32,
            fwhm_observed=0.25,
            fwhm_instrumental=0.08,
            use_cubic_habit=True
        )
        
        # Check sample broadening (quadratic subtraction)
        # β_sample = √(0.0625 - 0.0064) = √0.0561 = 0.2369
        assert abs(result.fwhm_sample - 0.237) < 0.01
        
        # Check crystallite size (K_w = 0.855, expect ~34 nm)
        assert 30 < result.size_nm < 40
        
        # Check K value used (K_w for (111))
        assert abs(result.k_factor - 0.855) < 0.01
        
        # Check validity flag
        assert result.validity_flag == ValidityFlag.VALID
    
    def test_spherical_gives_smaller_size(self):
        """Using spherical K=0.829 should give slightly smaller size."""
        result = calculate_scherrer(
            two_theta=43.32,
            fwhm_observed=0.25,
            fwhm_instrumental=0.08,
            use_cubic_habit=False
        )
        
        # K_spherical = 0.829, expect ~33 nm
        assert 28 < result.size_nm < 38
        assert abs(result.k_factor - 0.829) < 0.01


class TestValidityFlags:
    """Tests for validity flag system."""
    
    def test_valid_flag_normal(self):
        """Normal calculation should have VALID flag."""
        result = calculate_scherrer(
            two_theta=43.32,
            fwhm_observed=0.25,
            fwhm_instrumental=0.08
        )
        
        assert result.validity_flag == ValidityFlag.VALID
        assert result.is_reliable
    
    def test_unreliable_flag_narrow_peak(self):
        """Narrow peak (ratio < 1.2) should be UNRELIABLE."""
        # FWHM_obs = 0.09, FWHM_inst = 0.08 → ratio = 1.125 < 1.2
        result = calculate_scherrer(
            two_theta=43.32,
            fwhm_observed=0.09,
            fwhm_instrumental=0.08
        )
        
        assert result.validity_flag == ValidityFlag.UNRELIABLE
        assert not result.is_reliable
    
    def test_warning_flag_large_size(self):
        """Very large size (>200 nm) should trigger WARNING."""
        # Very narrow peak → large size
        result = calculate_scherrer(
            two_theta=43.32,
            fwhm_observed=0.02,
            fwhm_instrumental=0.001
        )
        
        # Either WARNING or UNRELIABLE
        assert result.validity_flag in [ValidityFlag.WARNING, ValidityFlag.UNRELIABLE]


class TestBatchCalculation:
    """Tests for batch processing."""
    
    def test_batch_multiple_peaks(self):
        """Calculate sizes for multiple peaks."""
        calc = ScherrerCalculator()
        
        peaks = [
            (43.32, 0.25),   # (111)
            (50.45, 0.28),   # (200)
            (74.16, 0.32),   # (220)
        ]
        
        results = calc.batch_calculate(peaks, fwhm_instrumental=0.08)
        
        assert len(results) == 3
        
        # Each result should have correct hkl
        assert results[0].hkl == (1, 1, 1)
        assert results[1].hkl == (2, 0, 0)
        assert results[2].hkl == (2, 2, 0)
    
    def test_average_size_calculation(self):
        """Test average size from multiple peaks."""
        calc = ScherrerCalculator()
        
        peaks = [
            (43.32, 0.25),
            (50.45, 0.28),
            (74.16, 0.32),
        ]
        
        results = calc.batch_calculate(peaks, fwhm_instrumental=0.08)
        avg, std = calc.average_size(results)
        
        assert avg > 0
        assert std >= 0


class TestReportGeneration:
    """Tests for report generation."""
    
    def test_report_contains_headers(self):
        """Report should contain expected headers."""
        calc = ScherrerCalculator()
        results = [
            calc.calculate(43.32, 0.25, 0.08),
            calc.calculate(50.45, 0.28, 0.08),
        ]
        
        report = generate_scherrer_report(results, "Test Sample")
        
        assert "Scherrer Crystallite Size Analysis" in report
        assert "Test Sample" in report
        assert "(111)" in report
        assert "(200)" in report


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
