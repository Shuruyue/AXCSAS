"""
Unit Tests for Peak Fitting Module
===================================

Tests quality metrics, hkl assignment, and fitting validation.

Run with: pytest tests/test_peak_fitting.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path


from axcsas.fitting import (
    PseudoVoigt,
    PseudoVoigtParams,
    calculate_r_wp,
    calculate_rss,
    calculate_r_squared,
    validate_fit_parameters,
    generate_quality_report,
    QualityLevel,
    assign_hkl,
    assign_hkl_detailed,
    format_hkl,
    JCPDS_COPPER_PEAKS,
)


class TestRwpCalculation:
    """Tests for R_wp weighted profile R-factor."""
    
    def test_perfect_fit_zero_rwp(self):
        """Perfect fit should have R_wp = 0."""
        observed = np.array([100, 200, 150, 180, 120])
        calculated = observed.copy()  # Perfect fit
        
        r_wp = calculate_r_wp(observed, calculated)
        
        assert r_wp == 0.0
    
    def test_poor_fit_high_rwp(self):
        """Poor fit should have high R_wp."""
        observed = np.array([100, 200, 150, 180, 120])
        calculated = np.array([50, 100, 75, 90, 60])  # 50% off
        
        r_wp = calculate_r_wp(observed, calculated)
        
        assert r_wp > 10  # Low confidence threshold
    
    def test_rwp_with_custom_weights(self):
        """Test R_wp with custom weights."""
        observed = np.array([100, 200, 150])
        calculated = np.array([95, 190, 145])
        weights = np.array([1.0, 2.0, 1.0])  # Double weight on middle
        
        r_wp = calculate_r_wp(observed, calculated, weights)
        
        assert r_wp > 0


class TestRSquaredCalculation:
    """Tests for R² coefficient of determination."""
    
    def test_perfect_fit_r_squared_one(self):
        """Perfect fit should have R² = 1."""
        observed = np.array([100, 200, 150, 180, 120])
        calculated = observed.copy()
        
        r_sq = calculate_r_squared(observed, calculated)
        
        assert abs(r_sq - 1.0) < 1e-10
    
    def test_poor_fit_low_r_squared(self):
        """Poor fit should have low R²."""
        observed = np.array([100, 200, 150, 180, 120])
        calculated = np.ones(5) * 150  # Constant (mean)
        
        r_sq = calculate_r_squared(observed, calculated)
        
        assert abs(r_sq) < 0.01  # Near zero


class TestParameterValidation:
    """Tests for fit parameter validation."""
    
    def test_valid_parameters_pass(self):
        """Valid parameters should pass validation."""
        is_valid, warnings = validate_fit_parameters(
            center=43.3,
            amplitude=1000,
            fwhm=0.25,
            eta=0.5
        )
        
        assert is_valid
        assert len(warnings) == 0
    
    def test_eta_out_of_bounds_fails(self):
        """η outside [0,1] should fail."""
        is_valid, warnings = validate_fit_parameters(
            center=43.3,
            amplitude=1000,
            fwhm=0.25,
            eta=1.5  # Out of bounds
        )
        
        assert not is_valid
        assert any("η" in w for w in warnings)
    
    def test_negative_fwhm_fails(self):
        """Negative FWHM should fail."""
        is_valid, warnings = validate_fit_parameters(
            center=43.3,
            amplitude=1000,
            fwhm=-0.1,  # Negative
            eta=0.5
        )
        
        assert not is_valid
        assert any("FWHM" in w for w in warnings)
    
    def test_zero_amplitude_fails(self):
        """Zero amplitude should fail."""
        is_valid, warnings = validate_fit_parameters(
            center=43.3,
            amplitude=0,  # Zero
            fwhm=0.25,
            eta=0.5
        )
        
        assert not is_valid


class TestQualityReport:
    """Tests for comprehensive quality report generation."""
    
    def test_excellent_quality(self):
        """Good fit should report excellent quality."""
        observed = np.array([100, 500, 1000, 500, 100])
        # Very close fit
        calculated = np.array([98, 495, 995, 502, 101])
        
        report = generate_quality_report(
            observed, calculated,
            center=43.3, amplitude=1000, fwhm=0.25, eta=0.5
        )
        
        assert report.is_valid
        assert report.r_wp < 10  # At least "good"
    
    def test_failed_quality_invalid_params(self):
        """Invalid parameters should report failed quality."""
        observed = np.array([100, 500, 1000, 500, 100])
        calculated = np.array([100, 500, 1000, 500, 100])
        
        report = generate_quality_report(
            observed, calculated,
            center=43.3, amplitude=1000, fwhm=-0.1, eta=0.5  # Invalid FWHM
        )
        
        assert not report.is_valid
        assert report.quality_level == QualityLevel.FAILED


class TestHklAssignment:
    """Tests for hkl Miller indices assignment."""
    
    def test_assign_111_peak(self):
        """43.3° should be assigned to (111)."""
        hkl = assign_hkl(43.3)
        
        assert hkl == (1, 1, 1)
    
    def test_assign_200_peak(self):
        """50.4° should be assigned to (200)."""
        hkl = assign_hkl(50.4)
        
        assert hkl == (2, 0, 0)
    
    def test_assign_220_peak(self):
        """74.1° should be assigned to (220)."""
        hkl = assign_hkl(74.1)
        
        assert hkl == (2, 2, 0)
    
    def test_assign_311_peak(self):
        """89.9° should be assigned to (311)."""
        hkl = assign_hkl(89.9)
        
        assert hkl == (3, 1, 1)
    
    def test_unassigned_peak(self):
        """30° should not be assigned to any Cu peak."""
        hkl = assign_hkl(30.0)
        
        assert hkl is None
    
    def test_detailed_assignment_confidence(self):
        """Detailed assignment should include confidence level."""
        assignment = assign_hkl_detailed(43.30)  # Very close to 43.297
        
        assert assignment.hkl == (1, 1, 1)
        assert assignment.confidence == "high"
    
    def test_format_hkl(self):
        """Format hkl tuple as string."""
        assert format_hkl((1, 1, 1)) == "(111)"
        assert format_hkl((2, 2, 0)) == "(220)"
        assert format_hkl(None) == "(?)"


class TestJCPDSData:
    """Tests for JCPDS reference data."""
    
    def test_jcpds_111_position(self):
        """JCPDS (111) should be at 43.297°."""
        assert abs(JCPDS_COPPER_PEAKS[(1, 1, 1)] - 43.297) < 0.001
    
    def test_jcpds_200_position(self):
        """JCPDS (200) should be at 50.433°."""
        assert abs(JCPDS_COPPER_PEAKS[(2, 0, 0)] - 50.433) < 0.001
    
    def test_all_peaks_present(self):
        """All main Cu peaks should be in JCPDS data."""
        expected_hkls = [(1,1,1), (2,0,0), (2,2,0), (3,1,1)]
        for hkl in expected_hkls:
            assert hkl in JCPDS_COPPER_PEAKS


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
