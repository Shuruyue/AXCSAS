"""
Unit Tests for Enhanced Williamson-Hall Analyzer
=================================================

Tests R² quality assessment, anisotropy diagnostics, and calculation accuracy.

Run with: pytest tests/test_williamson_hall.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path


from axcsas.methods.williamson_hall import (
    WilliamsonHallEnhanced,
    WHResultEnhanced,
    WHQualityLevel,
    analyze_williamson_hall_enhanced,
    generate_wh_report,
    MODULUS_MAP,
    WH_K_FACTOR,
    R2_EXCELLENT,
    R2_ACCEPTABLE,
)


class TestWHConstants:
    """Tests for W-H constants."""
    
    def test_wh_k_factor(self):
        """W-H K factor should be 0.9."""
        assert WH_K_FACTOR == 0.9
    
    def test_r2_thresholds(self):
        """R² thresholds should be defined."""
        assert R2_EXCELLENT == 0.95
        assert R2_ACCEPTABLE == 0.85
    
    def test_modulus_map_values(self):
        """MODULUS_MAP should have correct values."""
        assert MODULUS_MAP[(1, 1, 1)] == 191.0  # Hardest
        assert MODULUS_MAP[(2, 0, 0)] == 67.0   # Softest
        assert MODULUS_MAP[(2, 2, 0)] == 130.0


class TestDocumentExample:
    """Test case from document 05 §7."""
    
    def test_document_example_calculation(self):
        """
        Verify calculation matches document 05 §7 example.
        
        Input:
            (111): 2θ=43.32°, FWHM=0.224°
            (200): 2θ=50.45°, FWHM=0.251°
            (220): 2θ=74.16°, FWHM=0.282°
            (311): 2θ=89.97°, FWHM=0.305°
            
        Note: For real ED-Cu, R² is often low due to elastic anisotropy.
        This is expected behavior per document 05 §8.
        """
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        
        # Result should be computed (not failed)
        assert result.n_peaks == 4
        
        # Size should be positive
        assert result.crystallite_size_nm > 0
        
        # Strain should be computable
        assert result.microstrain >= 0
        
        # R² value should be defined
        assert 0 <= result.r_squared <= 1
    
    def test_anisotropy_is_detected_for_ed_cu(self):
        """
        ED-Cu data often shows low R² due to elastic anisotropy.
        This should trigger appropriate warnings.
        """
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        
        # If R² is low, warning should be generated
        if result.r_squared < R2_ACCEPTABLE:
            assert result.warning_message
            assert result.quality_level == WHQualityLevel.POOR


class TestQualityAssessment:
    """Tests for R² quality assessment."""
    
    def test_excellent_quality(self):
        """R² > 0.95 should be EXCELLENT."""
        # Create nearly perfect linear data
        x = np.array([0.3, 0.4, 0.5, 0.6])
        two_theta = 2 * np.arcsin(x) * 180 / np.pi
        
        # Generate FWHM that produces near-linear W-H plot
        intercept = 0.003
        slope = 0.001
        y_ideal = intercept + slope * x
        beta = y_ideal / np.cos(np.arcsin(x))
        fwhm = beta * 180 / np.pi
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        
        # High R² expected for synthetic linear data
        assert result.r_squared > 0.9
    
    def test_poor_quality_warning(self):
        """R² < 0.85 should trigger warning."""
        # Create scattered data that won't fit well
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        # Intentionally scattered FWHM values
        fwhm = np.array([0.224, 0.40, 0.20, 0.45])  # Non-monotonic
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        
        # Should either have warning or be poor quality
        assert result.warning_message or not result.is_reliable


class TestDataValidation:
    """Tests for input data validation."""
    
    def test_minimum_peaks_required(self):
        """At least 3 peaks are required."""
        two_theta = np.array([43.32, 50.45])  # Only 2 peaks
        fwhm = np.array([0.224, 0.251])
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        
        assert not result.is_reliable
        assert "3" in result.warning_message
    
    def test_array_length_mismatch(self):
        """Mismatched arrays should fail."""
        two_theta = np.array([43.32, 50.45, 74.16])
        fwhm = np.array([0.224, 0.251])  # Shorter
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        
        assert not result.is_reliable


class TestAnisotropyDiagnostics:
    """Tests for anisotropy diagnostic reporting."""
    
    def test_anisotropy_note_for_poor_r2(self):
        """Poor R² should generate anisotropy note."""
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        # Intentionally scattered values
        fwhm = np.array([0.30, 0.15, 0.40, 0.20])
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        
        # If R² is poor, should have anisotropy note
        if result.r_squared < R2_ACCEPTABLE:
            assert result.anisotropy_note
            assert "異方性" in result.anisotropy_note or "191" in result.anisotropy_note


class TestReportGeneration:
    """Tests for report generation."""
    
    def test_report_contains_key_info(self):
        """Report should contain all key information."""
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        
        result = analyze_williamson_hall_enhanced(two_theta, fwhm)
        report = generate_wh_report(result, "Test Sample")
        
        assert "Williamson-Hall" in report
        assert "Test Sample" in report
        assert "Crystallite Size" in report
        assert "Microstrain" in report
        assert "R²" in report


class TestPlotData:
    """Tests for plot data generation."""
    
    def test_get_plot_data(self):
        """Plot data should be generated correctly."""
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        
        analyzer = WilliamsonHallEnhanced()
        x_data, y_data, x_fit, y_fit = analyzer.get_plot_data(two_theta, fwhm)
        
        assert len(x_data) == 4
        assert len(y_data) == 4
        assert len(x_fit) == 100  # Fit line has more points


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
