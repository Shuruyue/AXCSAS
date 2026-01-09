"""
Unit Tests for Preprocessing Pipeline
======================================

Tests validation, smoothing, background, Kα2 stripping, and pipeline integration.

Run with: pytest tests/test_preprocessing_pipeline.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from preprocessing import (
    XRDDataset,
    DataValidationResult,
    validate_xrd_data,
    check_negative_values,
    SavitzkyGolayFilter,
    BackgroundSubtractor,
    KalphaStripper,
    PreprocessingPipeline,
    PreprocessingResult,
    should_apply_kalpha_stripping,
)


class TestXRDDataset:
    """Tests for XRDDataset container."""
    
    def test_dataset_creation(self):
        """Test basic dataset creation."""
        two_theta = np.linspace(20, 80, 1000)
        intensity = np.random.rand(1000) * 100
        
        ds = XRDDataset(two_theta, intensity)
        
        assert ds.n_points == 1000
        assert ds.theta_range == (20.0, 80.0)
        
    def test_step_size_calculation(self):
        """Test step size property."""
        two_theta = np.arange(20, 80, 0.02)
        intensity = np.ones_like(two_theta)
        
        ds = XRDDataset(two_theta, intensity)
        
        assert abs(ds.step_size - 0.02) < 0.001


class TestValidation:
    """Tests for data validation functions."""
    
    def test_valid_data_passes(self):
        """Test that good data passes validation."""
        two_theta = np.linspace(20, 80, 1000)
        intensity = np.abs(np.random.randn(1000)) * 100 + 10
        
        result = validate_xrd_data(two_theta, intensity)
        
        assert result.is_valid
        
    def test_negative_intensity_detected(self):
        """Test detection of negative intensity values."""
        two_theta = np.linspace(20, 80, 1000)
        intensity = np.random.randn(1000) * 100  # Has negative values
        
        result = validate_xrd_data(two_theta, intensity)
        
        assert any(w.code == "NEGATIVE_INTENSITY" for w in result.warnings)
        
    def test_low_angle_warning(self):
        """Test warning for data below 10°."""
        two_theta = np.linspace(5, 80, 1000)
        intensity = np.abs(np.random.randn(1000)) * 100
        
        result = validate_xrd_data(two_theta, intensity)
        
        assert any(w.code == "LOW_ANGLE" for w in result.warnings)
        
    def test_insufficient_points(self):
        """Test warning for too few data points."""
        two_theta = np.linspace(20, 80, 50)  # Only 50 points
        intensity = np.ones(50)
        
        result = validate_xrd_data(two_theta, intensity)
        
        assert any(w.code == "INSUFFICIENT_POINTS" for w in result.warnings)
        
    def test_non_monotonic_fails(self):
        """Test that non-monotonic 2θ fails validation."""
        two_theta = np.array([20, 30, 25, 40, 50])  # Not monotonic
        intensity = np.ones(5)
        
        result = validate_xrd_data(two_theta, intensity)
        
        assert not result.is_valid
        assert any(w.code == "NON_MONOTONIC" for w in result.warnings)


class TestNegativeValueHandling:
    """Tests for negative value detection and correction."""
    
    def test_negative_detection(self):
        """Test detection of negative values."""
        intensity = np.array([100, -10, 200, -5, 150])
        
        corrected, indices = check_negative_values(intensity, auto_correct=False)
        
        assert indices == [1, 3]
        assert np.array_equal(corrected, intensity)  # No correction
        
    def test_auto_correction(self):
        """Test automatic correction of negative values."""
        intensity = np.array([100, -10, 200, -5, 150])
        
        corrected, indices = check_negative_values(intensity, auto_correct=True)
        
        assert indices == [1, 3]
        assert corrected[1] == 0
        assert corrected[3] == 0


class TestSmoothingBoundaryCheck:
    """Tests for smoothing window boundary validation."""
    
    def test_window_larger_than_data_fails(self):
        """Test that window > data length raises error."""
        intensity = np.array([1, 2, 3, 4, 5])  # Only 5 points
        sg = SavitzkyGolayFilter(window_size=11)  # Window = 11
        
        with pytest.raises(ValueError, match="must be less than"):
            sg.apply(intensity)


class TestKalphaStrippingCondition:
    """Tests for Kα2 stripping decision logic."""
    
    def test_high_angle_triggers_stripping(self):
        """Test that 2θ > 40° triggers stripping."""
        two_theta = np.linspace(20, 80, 1000)
        
        assert should_apply_kalpha_stripping(two_theta) is True
        
    def test_low_angle_skips_stripping(self):
        """Test that 2θ < 40° skips stripping."""
        two_theta = np.linspace(20, 35, 1000)
        
        assert should_apply_kalpha_stripping(two_theta) is False


class TestPreprocessingPipeline:
    """Tests for complete pipeline."""
    
    def test_pipeline_runs_without_error(self):
        """Test that pipeline completes on valid data."""
        two_theta = np.linspace(20, 80, 1000)
        # Simulate a simple XRD pattern
        intensity = 100 + 50 * np.exp(-((two_theta - 43)**2) / 2)
        intensity += np.random.randn(1000) * 5  # Add noise
        
        pipeline = PreprocessingPipeline()
        result = pipeline.run(two_theta, intensity)
        
        assert isinstance(result, PreprocessingResult)
        assert len(result.two_theta) == 1000
        assert len(result.intensity) == 1000
        
    def test_pipeline_records_steps(self):
        """Test that processing steps are recorded."""
        two_theta = np.linspace(20, 80, 1000)
        intensity = np.abs(np.random.randn(1000)) * 100
        
        pipeline = PreprocessingPipeline()
        result = pipeline.run(two_theta, intensity)
        
        # Should have at least: validation, smoothing, background
        assert len(result.steps) >= 3
        
        step_names = [s.name for s in result.steps]
        assert "Data Validation" in step_names
        assert "Savitzky-Golay Smoothing" in step_names
        
    def test_pipeline_summary_generation(self):
        """Test that summary string is generated."""
        two_theta = np.linspace(20, 80, 1000)
        intensity = np.abs(np.random.randn(1000)) * 100
        
        pipeline = PreprocessingPipeline()
        result = pipeline.run(two_theta, intensity)
        summary = result.summary()
        
        assert "Preprocessing Summary" in summary
        assert "Steps performed" in summary


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
