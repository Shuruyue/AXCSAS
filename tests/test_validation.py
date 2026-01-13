"""
Unit Tests for Validation Module
Tests goodness of fit and error analysis functions.

Run with: pytest tests/test_validation.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path


from axcsas.validation import (
    calculate_rwp,
    calculate_r_squared,
    calculate_chi_squared,
    FitQualityResult,
    ErrorAnalyzer,
    validate_size_range,
    check_broadening_ratio
)
from axcsas.validation.goodness_of_fit import assess_fit_quality, calculate_rp


class TestGoodnessOfFit:
    """Tests for fit quality metrics."""
    
    def test_rwp_perfect_fit(self):
        """Test Rwp for perfect fit."""
        observed = np.array([100, 200, 300, 400, 500])
        calculated = observed.copy()
        
        rwp = calculate_rwp(observed, calculated)
        
        assert rwp == 0.0
        
    def test_rwp_with_residuals(self):
        """Test Rwp with known residuals."""
        observed = np.array([100, 200, 300, 400])
        calculated = np.array([110, 190, 310, 390])  # ±10 error
        
        rwp = calculate_rwp(observed, calculated)
        
        assert 0 < rwp < 20  # Should be small but non-zero
        
    def test_r_squared_perfect_fit(self):
        """Test R² for perfect fit."""
        observed = np.array([100, 200, 300, 400])
        calculated = observed.copy()
        
        r2 = calculate_r_squared(observed, calculated)
        
        assert r2 == 1.0
        
    def test_r_squared_known_value(self):
        """Test R² calculation."""
        observed = np.array([1, 2, 3, 4, 5])
        calculated = np.array([1.1, 1.9, 3.2, 3.8, 5.1])
        
        r2 = calculate_r_squared(observed, calculated)
        
        assert 0.95 < r2 < 1.0  # Good but not perfect
        
    def test_chi_squared(self):
        """Test chi-squared calculation."""
        observed = np.array([100, 200, 300])
        calculated = np.array([105, 195, 305])
        
        chi_sq, reduced_chi_sq = calculate_chi_squared(
            observed, calculated, n_params=1
        )
        
        assert chi_sq > 0
        assert reduced_chi_sq > 0
        
    def test_assess_fit_quality(self):
        """Test comprehensive fit assessment."""
        observed = np.array([100, 200, 300, 400, 500])
        calculated = np.array([102, 198, 303, 398, 502])
        
        result = assess_fit_quality(observed, calculated)
        
        assert isinstance(result, FitQualityResult)
        assert result.rwp < 10  # Should be acceptable
        assert result.r_squared > 0.99
        assert result.is_acceptable


class TestErrorAnalyzer:
    """Tests for error analysis and validation."""
    
    def test_validate_size_normal_range(self):
        """Test size validation for normal values."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_size(50.0)  # 50 nm
        
        assert result.is_valid
        assert len(result.warnings) == 0
        
    def test_validate_size_too_small(self):
        """Test size validation for small values."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_size(1.0)  # 1 nm - too small
        
        assert not result.is_valid
        assert len(result.warnings) > 0
        
    def test_validate_size_too_large(self):
        """Test size validation for large values."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_size(300.0)  # 300 nm - too large
        
        assert not result.is_valid
        assert len(result.warnings) > 0
        
    def test_validate_broadening_sufficient(self):
        """Test broadening ratio validation for good case."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_broadening(0.3, 0.1)  # ratio = 3.0
        
        assert result.is_valid
        
    def test_validate_broadening_insufficient(self):
        """Test broadening ratio validation for bad case."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_broadening(0.11, 0.1)  # ratio = 1.1
        
        assert not result.is_valid
        
    def test_validate_broadening_impossible(self):
        """Test when observed < instrumental."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_broadening(0.08, 0.1)
        
        assert not result.is_valid
        assert any(w.level.value == 'critical' for w in result.warnings)
        
    def test_validate_fit_quality_good(self):
        """Test fit quality validation for good values."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_fit_quality(rwp=5.0, r_squared=0.98)
        
        assert result.is_valid
        
    def test_validate_fit_quality_bad(self):
        """Test fit quality validation for poor values."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_fit_quality(rwp=15.0, r_squared=0.80)
        
        assert not result.is_valid
        
    def test_validate_all(self):
        """Test comprehensive validation."""
        analyzer = ErrorAnalyzer()
        
        result = analyzer.validate_all(
            size_nm=50.0,
            fwhm_observed=0.3,
            fwhm_instrumental=0.1,
            rwp=5.0,
            r_squared=0.98
        )
        
        assert result.is_valid


class TestConvenienceFunctions:
    """Tests for convenience validation functions."""
    
    def test_validate_size_range(self):
        """Test quick size range check."""
        assert validate_size_range(50.0)
        assert not validate_size_range(1.0)
        assert not validate_size_range(300.0)
        
    def test_check_broadening_ratio(self):
        """Test quick broadening ratio check."""
        assert check_broadening_ratio(0.3, 0.1)
        assert not check_broadening_ratio(0.11, 0.1)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
