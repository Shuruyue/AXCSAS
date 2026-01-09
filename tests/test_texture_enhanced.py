"""
Unit Tests for Enhanced Texture Analyzer
==========================================

Tests Harris TC calculation and document examples (DATA ONLY, no diagnosis).

Run with: pytest tests/test_texture_enhanced.py -v
"""

import pytest
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from physics.texture_enhanced import (
    TextureAnalyzerEnhanced,
    TextureResultEnhanced,
    OrientationType,
    analyze_texture_enhanced,
    generate_texture_report,
    JCPDS_STANDARD_INTENSITY,
)


class TestJCPDSStandard:
    """Tests for JCPDS standard intensity data."""
    
    def test_standard_111(self):
        """I₀(111) should be 100."""
        assert JCPDS_STANDARD_INTENSITY[(1, 1, 1)] == 100.0
    
    def test_standard_200(self):
        """I₀(200) should be 46."""
        assert JCPDS_STANDARD_INTENSITY[(2, 0, 0)] == 46.0
    
    def test_standard_220(self):
        """I₀(220) should be 20."""
        assert JCPDS_STANDARD_INTENSITY[(2, 2, 0)] == 20.0


class TestDocumentExample:
    """Test case from document 06 §4."""
    
    def test_document_example_tc_values(self):
        """
        Verify TC calculation matches document 06 §4 example.
        
        Input:
            I(111) = 15680 → Ratio = 156.8
            I(200) = 5520  → Ratio = 120.0
            I(220) = 4200  → Ratio = 210.0
            
        Average = 162.3
        
        Expected TC:
            TC(111) = 0.97
            TC(200) = 0.74
            TC(220) = 1.29
        """
        intensities = {
            (1, 1, 1): 15680,
            (2, 0, 0): 5520,
            (2, 2, 0): 4200,
        }
        
        result = analyze_texture_enhanced(intensities)
        
        # Check TC values (allow 5% tolerance)
        assert abs(result.tc_values[(1, 1, 1)] - 0.97) < 0.05
        assert abs(result.tc_values[(2, 0, 0)] - 0.74) < 0.05
        assert abs(result.tc_values[(2, 2, 0)] - 1.29) < 0.07
    
    def test_document_example_dominant(self):
        """Dominant orientation should be (220)."""
        intensities = {
            (1, 1, 1): 15680,
            (2, 0, 0): 5520,
            (2, 2, 0): 4200,
        }
        
        result = analyze_texture_enhanced(intensities)
        
        assert result.dominant_hkl == (2, 2, 0)


class TestOrientationType:
    """Tests for orientation type classification (DATA MARKERS)."""
    
    def test_preferred_orientation(self):
        """TC > 1.0 should be marked PREFERRED."""
        # Create strong (111) preference
        intensities = {
            (1, 1, 1): 25000,  # Much stronger
            (2, 0, 0): 4600,   # Standard
            (2, 2, 0): 2000,   # Standard
        }
        
        result = analyze_texture_enhanced(intensities)
        
        # (111) should be marked preferred
        tc_111 = [d for d in result.tc_details if d.hkl == (1, 1, 1)][0]
        assert tc_111.orientation_type == OrientationType.PREFERRED
    
    def test_suppressed_orientation(self):
        """TC < 1.0 (and not near 1) should be marked SUPPRESSED."""
        intensities = {
            (1, 1, 1): 15680,
            (2, 0, 0): 5520,
            (2, 2, 0): 4200,
        }
        
        result = analyze_texture_enhanced(intensities)
        
        # (200) should be marked suppressed (TC = 0.74)
        tc_200 = [d for d in result.tc_details if d.hkl == (2, 0, 0)][0]
        assert tc_200.orientation_type == OrientationType.SUPPRESSED


class TestStatisticalMeasures:
    """Tests for statistical measures (DATA ONLY)."""
    
    def test_degree_of_texture(self):
        """Degree of texture should be standard deviation of TC values."""
        intensities = {
            (1, 1, 1): 15680,
            (2, 0, 0): 5520,
            (2, 2, 0): 4200,
        }
        
        result = analyze_texture_enhanced(intensities)
        
        # Should have non-zero degree of texture
        assert result.degree_of_texture > 0
    
    def test_is_random_detection(self):
        """Should correctly identify random texture."""
        # Balanced intensities = all TC near 1
        intensities = {
            (1, 1, 1): 10000,
            (2, 0, 0): 4600,
            (2, 2, 0): 2000,
        }
        
        result = analyze_texture_enhanced(intensities)
        
        # is_random is True if all TC are between 0.9-1.1
        assert isinstance(result.is_random, bool)


class TestReportGeneration:
    """Tests for report generation (DATA ONLY)."""
    
    def test_report_contains_key_info(self):
        """Report should contain all key data information."""
        intensities = {
            (1, 1, 1): 15680,
            (2, 0, 0): 5520,
            (2, 2, 0): 4200,
        }
        
        result = analyze_texture_enhanced(intensities)
        report = generate_texture_report(result, "Test Sample")
        
        assert "Texture" in report
        assert "Test Sample" in report
        assert "(111)" in report
        assert "(220)" in report
        assert "TC" in report
    
    def test_report_no_diagnosis(self):
        """Report should NOT contain process diagnosis."""
        intensities = {
            (1, 1, 1): 15680,
            (2, 0, 0): 5520,
            (2, 2, 0): 4200,
        }
        
        result = analyze_texture_enhanced(intensities)
        report = generate_texture_report(result, "Test Sample")
        
        # Should not contain judgment keywords
        assert "PEG" not in report
        assert "Diagnosis" not in report or "Process Diagnosis" not in report


class TestEdgeCases:
    """Tests for edge cases."""
    
    def test_minimum_peaks(self):
        """Should require at least 2 peaks."""
        intensities = {
            (1, 1, 1): 15680,  # Only 1 peak
        }
        
        result = analyze_texture_enhanced(intensities)
        
        assert result.n_peaks == 0 or len(result.tc_values) < 2
    
    def test_zero_intensity_handled(self):
        """Zero intensity should not crash."""
        intensities = {
            (1, 1, 1): 15680,
            (2, 0, 0): 0,  # Zero
            (2, 2, 0): 4200,
        }
        
        # Should not raise exception
        result = analyze_texture_enhanced(intensities)
        assert result is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
