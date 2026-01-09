"""
Unit Tests for Enhanced Scherrer Calculator
============================================

Tests K value lookup, unit conversion, validity flags, and calculation accuracy.

Run with: pytest tests/test_scherrer_enhanced.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from physics.scherrer_enhanced import (
    ScherrerCalculatorEnhanced,
    ScherrerResultEnhanced,
    ValidityFlag,
    calculate_scherrer_enhanced,
    generate_scherrer_report,
    FWHM_RATIO_THRESHOLD,
)
from core.copper_crystal import (
    get_k_for_hkl,
    SCHERRER_CUBIC_K,
)


class TestKValueLookup:
    """Tests for dynamic K value selection."""
    
    def test_k_111_value(self):
        """K(111) should be 1.155."""
        k = get_k_for_hkl(1, 1, 1)
        assert abs(k - 1.155) < 0.001
    
    def test_k_200_value(self):
        """K(200) should be 1.000."""
        k = get_k_for_hkl(2, 0, 0)
        assert abs(k - 1.000) < 0.001
    
    def test_k_220_value_updated(self):
        """K(220) should be 1.061 (not 0.707)."""
        k = get_k_for_hkl(2, 2, 0)
        assert abs(k - 1.061) < 0.001  # Updated value
    
    def test_k_311_value_updated(self):
        """K(311) should be 1.116 (not 0.89)."""
        k = get_k_for_hkl(3, 1, 1)
        assert abs(k - 1.116) < 0.001  # Updated value
    
    def test_k_spherical_fallback(self):
        """Non-cubic habit should use K=0.89."""
        k = get_k_for_hkl(1, 1, 1, use_cubic_habit=False)
        assert abs(k - 0.89) < 0.001
    
    def test_scherrer_cubic_k_class_updated(self):
        """ScherrerCubicK class should have updated values."""
        assert abs(SCHERRER_CUBIC_K.K_220 - 1.061) < 0.001
        assert abs(SCHERRER_CUBIC_K.K_311 - 1.116) < 0.001


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
    """Test case from document 04 §4."""
    
    def test_document_example_calculation(self):
        """
        Verify calculation matches document 04 §4 example.
        
        Input:
            2θ = 43.32°
            FWHM_obs = 0.25°
            FWHM_inst = 0.08°
            
        Expected:
            β_sample = 0.224°
            D = 49.0 nm
        """
        result = calculate_scherrer_enhanced(
            two_theta=43.32,
            fwhm_observed=0.25,
            fwhm_instrumental=0.08,
            use_cubic_habit=True
        )
        
        # Check sample broadening
        assert abs(result.fwhm_sample - 0.224) < 0.01
        
        # Check crystallite size (allow 5% tolerance)
        assert abs(result.size_nm - 49.0) < 2.5
        
        # Check K value used
        assert abs(result.k_factor - 1.155) < 0.01
        
        # Check validity flag
        assert result.validity_flag == ValidityFlag.VALID
    
    def test_spherical_gives_smaller_size(self):
        """Using spherical K=0.89 should give ~37.8 nm (30% less)."""
        result = calculate_scherrer_enhanced(
            two_theta=43.32,
            fwhm_observed=0.25,
            fwhm_instrumental=0.08,
            use_cubic_habit=False
        )
        
        # Should be ~37.8 nm per document
        assert abs(result.size_nm - 37.8) < 3.0


class TestValidityFlags:
    """Tests for validity flag system."""
    
    def test_valid_flag_normal(self):
        """Normal calculation should have VALID flag."""
        result = calculate_scherrer_enhanced(
            two_theta=43.32,
            fwhm_observed=0.25,
            fwhm_instrumental=0.08
        )
        
        assert result.validity_flag == ValidityFlag.VALID
        assert result.is_reliable
    
    def test_unreliable_flag_narrow_peak(self):
        """Narrow peak (ratio < 1.2) should be UNRELIABLE."""
        # FWHM_obs = 0.09, FWHM_inst = 0.08 → ratio = 1.125 < 1.2
        result = calculate_scherrer_enhanced(
            two_theta=43.32,
            fwhm_observed=0.09,
            fwhm_instrumental=0.08
        )
        
        assert result.validity_flag == ValidityFlag.UNRELIABLE
        assert not result.is_reliable
    
    def test_warning_flag_large_size(self):
        """Very large size (>200 nm) should trigger WARNING."""
        # Very narrow peak → large size
        result = calculate_scherrer_enhanced(
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
        calc = ScherrerCalculatorEnhanced()
        
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
        calc = ScherrerCalculatorEnhanced()
        
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
        calc = ScherrerCalculatorEnhanced()
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
