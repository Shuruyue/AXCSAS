"""
Unit Tests for Physics Module
Tests Scherrer, Williamson-Hall, and Texture calculations.

Run with: pytest tests/test_physics.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from physics import (
    ScherrerCalculator,
    WilliamsonHallAnalyzer,
    TextureAnalyzer,
    CagliotiCorrection,
    calculate_crystallite_size,
    calculate_texture_coefficient
)


class TestScherrerCalculator:
    """Tests for Scherrer equation calculations."""
    
    def test_basic_calculation(self):
        """Test basic crystallite size calculation."""
        calc = ScherrerCalculator()
        result = calc.calculate(two_theta=43.3, fwhm=0.25)
        
        # Size should be positive and reasonable
        assert result.size_nm > 0
        assert 10 < result.size_nm < 100  # Reasonable range for this FWHM
        
    def test_known_values(self):
        """Test against known reference values."""
        calc = ScherrerCalculator(wavelength=1.54056, k_factor=0.89)
        
        # For 2θ=43.3°, FWHM=0.5°, expected ~17 nm
        result = calc.calculate(two_theta=43.3, fwhm=0.5)
        assert 15 < result.size_nm < 20
        
    def test_reliability_flags(self):
        """Test reliability detection for extreme sizes."""
        calc = ScherrerCalculator()
        
        # Very small FWHM -> large size -> should be flagged
        result = calc.calculate(two_theta=43.3, fwhm=0.02)
        assert result.size_nm > 200
        assert not result.is_reliable
        
        # Very large FWHM -> small size -> should be flagged
        result = calc.calculate(two_theta=43.3, fwhm=5.0)
        assert result.size_nm < 2
        assert not result.is_reliable
        
    def test_fwhm_units(self):
        """Test both degree and radian inputs."""
        calc = ScherrerCalculator()
        
        fwhm_deg = 0.25
        fwhm_rad = np.radians(fwhm_deg)
        
        result_deg = calc.calculate(43.3, fwhm_deg, fwhm_is_radians=False)
        result_rad = calc.calculate(43.3, fwhm_rad, fwhm_is_radians=True)
        
        assert abs(result_deg.size_nm - result_rad.size_nm) < 0.1
        
    def test_batch_calculation(self):
        """Test batch processing of multiple peaks."""
        calc = ScherrerCalculator()
        
        peaks = [(43.3, 0.25), (50.4, 0.28), (74.1, 0.35)]
        results = calc.batch_calculate(peaks)
        
        assert len(results) == 3
        assert all(r.size_nm > 0 for r in results)
        
    def test_invalid_inputs(self):
        """Test error handling for invalid inputs."""
        calc = ScherrerCalculator()
        
        with pytest.raises(ValueError):
            calc.calculate(43.3, -0.25)  # Negative FWHM
            
        with pytest.raises(ValueError):
            calc.calculate(200, 0.25)  # Invalid 2θ


class TestWilliamsonHallAnalyzer:
    """Tests for Williamson-Hall analysis."""
    
    def test_basic_analysis(self):
        """Test basic W-H analysis with synthetic data."""
        wh = WilliamsonHallAnalyzer()
        
        # Sample data with known linear relationship
        two_theta = np.array([43.3, 50.4, 74.1, 89.9])
        fwhm = np.array([0.25, 0.27, 0.35, 0.42])
        
        result = wh.analyze(two_theta, fwhm)
        
        assert result.crystallite_size_nm > 0
        assert result.microstrain >= 0
        assert 0 <= result.r_squared <= 1
        
    def test_zero_strain_case(self):
        """Test case with no strain broadening."""
        wh = WilliamsonHallAnalyzer(wavelength=1.54056, k_factor=0.89)
        
        # Pure size broadening (constant β)
        two_theta = np.array([43.3, 50.4, 74.1, 89.9])
        
        # Calculate β that gives constant βcosθ (no strain)
        theta = np.radians(two_theta / 2)
        target_intercept = 0.01  # Constant
        fwhm = target_intercept / np.cos(theta)
        fwhm_deg = np.degrees(fwhm)
        
        result = wh.analyze(two_theta, fwhm_deg)
        
        # Strain should be near zero
        assert abs(result.microstrain) < 0.01
        
    def test_reliability_with_few_peaks(self):
        """Test warning for insufficient peaks."""
        wh = WilliamsonHallAnalyzer()
        
        # Only 2 peaks - should flag as unreliable
        two_theta = np.array([43.3, 50.4])
        fwhm = np.array([0.25, 0.27])
        
        result = wh.analyze(two_theta, fwhm)
        
        assert not result.is_reliable
        assert result.warning is not None
        
    def test_plot_data(self):
        """Test plot data generation."""
        wh = WilliamsonHallAnalyzer()
        
        two_theta = np.array([43.3, 50.4, 74.1, 89.9])
        fwhm = np.array([0.25, 0.27, 0.35, 0.42])
        
        x_data, y_data, x_fit, y_fit = wh.get_plot_data(two_theta, fwhm)
        
        assert len(x_data) == len(two_theta)
        assert len(x_fit) == 2  # Line defined by 2 points


class TestTextureAnalyzer:
    """Tests for texture coefficient analysis."""
    
    def test_random_texture(self):
        """Test detection of random texture."""
        analyzer = TextureAnalyzer()
        
        # Intensities matching standard ratios = random texture
        intensities = {
            (1, 1, 1): 100,
            (2, 0, 0): 46,
            (2, 2, 0): 20,
        }
        
        result = analyzer.analyze(intensities)
        
        # TC should be approximately 1 for all peaks
        for tc in result.tc_values.values():
            assert 0.9 < tc < 1.1
            
        assert result.is_random
        
    def test_preferred_orientation(self):
        """Test detection of (111) preferred orientation."""
        analyzer = TextureAnalyzer()
        
        # Enhanced (111) intensity
        intensities = {
            (1, 1, 1): 200,  # 2× standard
            (2, 0, 0): 23,   # 0.5× standard
            (2, 2, 0): 10,   # 0.5× standard
        }
        
        result = analyzer.analyze(intensities)
        
        assert result.tc_values[(1, 1, 1)] > 1.5
        assert result.preferred_orientation == (1, 1, 1)
        assert not result.is_random
        
    def test_tc_sum_equals_n(self):
        """Test that sum of TC equals number of reflections."""
        analyzer = TextureAnalyzer()
        
        intensities = {
            (1, 1, 1): 150,
            (2, 0, 0): 30,
            (2, 2, 0): 25,
        }
        
        result = analyzer.analyze(intensities)
        
        tc_sum = sum(result.tc_values.values())
        n = len(result.tc_values)
        
        assert abs(tc_sum - n) < 0.01
        
    def test_hkl_matching(self):
        """Test automatic hkl assignment from angles."""
        analyzer = TextureAnalyzer()
        
        # Test Cu (111) at 43.3°
        hkl = analyzer.get_hkl_for_angle(43.3, tolerance=0.5)
        assert hkl == (1, 1, 1)
        
        # Test Cu (200) at 50.4°
        hkl = analyzer.get_hkl_for_angle(50.4, tolerance=0.5)
        assert hkl == (2, 0, 0)


class TestCagliotiCorrection:
    """Tests for Caglioti instrumental correction."""
    
    def test_calibration(self):
        """Test Caglioti parameter calibration."""
        caglioti = CagliotiCorrection()
        
        # Synthetic standard data
        two_theta = np.array([30, 50, 70, 90, 110])
        fwhm = np.array([0.08, 0.10, 0.13, 0.17, 0.22])
        
        params = caglioti.calibrate(two_theta, fwhm)
        
        assert caglioti.is_calibrated
        assert params.U is not None
        assert params.V is not None
        assert params.W is not None
        
    def test_fwhm_calculation(self):
        """Test instrumental FWHM calculation."""
        caglioti = CagliotiCorrection(U=0.01, V=-0.01, W=0.005)
        
        fwhm = caglioti.calculate_fwhm_inst(45.0)
        
        assert fwhm > 0
        assert fwhm < 0.5  # Reasonable instrumental broadening
        
    def test_broadening_correction(self):
        """Test broadening correction logic."""
        caglioti = CagliotiCorrection(U=0.005, V=-0.005, W=0.003)
        
        fwhm_obs = 0.3
        two_theta = 45.0
        
        fwhm_corr, is_reliable, warning = caglioti.correct_broadening(fwhm_obs, two_theta)
        
        assert fwhm_corr > 0
        assert fwhm_corr < fwhm_obs  # Corrected should be smaller
        assert isinstance(warning, str)  # Warning should be a string


class TestConvenienceFunctions:
    """Tests for convenience functions."""
    
    def test_calculate_crystallite_size(self):
        """Test quick size calculation function."""
        size = calculate_crystallite_size(43.3, 0.25)
        
        assert 20 < size < 50
        
    def test_calculate_texture_coefficient(self):
        """Test quick TC calculation function."""
        intensities = {
            (1, 1, 1): 100,
            (2, 0, 0): 46,
        }
        
        tc = calculate_texture_coefficient(intensities)
        
        assert (1, 1, 1) in tc
        assert (2, 0, 0) in tc


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
