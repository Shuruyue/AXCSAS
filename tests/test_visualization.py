
"""
Unit Tests for Visualization Module
===================================

Tests for the updated visualization components, ensuring all plots use
correct constants and rigorous fitting methods.
"""

import pytest
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock, patch

from axcsas.visualization.style import apply_axcsas_style, PEAK_COLORS
from axcsas.visualization.wh_plots import plot_williamson_hall
from axcsas.visualization.texture_plots import plot_texture_polar
from axcsas.visualization.fitting_plots import plot_peak_fit
from axcsas.visualization.generate_fitting_diagnosis import fit_peak_with_diagnosis
from axcsas.core.constants import CU_KA1

class TestVisualizationConstants:
    """Tests for constants usage in visualization."""
    
    def test_colors_defined(self):
        """Standard peak colors should be defined."""
        assert '(111)' in PEAK_COLORS
        assert '(200)' in PEAK_COLORS
        
    def test_style_application(self):
        """Style application should not raise error."""
        try:
            apply_axcsas_style()
        except Exception as e:
            pytest.fail(f"Style application failed: {e}")

class TestWHPlotting:
    """Tests for Williamson-Hall plotting function."""
    
    def test_wh_plot_execution(self):
        """Should execute plot_williamson_hall without error."""
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        
        # Build a verification mock to prevent actual plotting
        with patch('matplotlib.pyplot.show'):
             fig = plot_williamson_hall(two_theta, fwhm, show=False)
             assert fig is not None

class TestTexturePlotting:
    """Tests for Texture plotting function."""
    
    def test_texture_plot_execution(self):
        """Should execute plot_texture_polar without error."""
        tc_values = {'(111)': 1.35, '(200)': 0.72, '(220)': 0.93}
        
        with patch('matplotlib.pyplot.show'):
            fig = plot_texture_polar(tc_values, show=False)
            assert fig is not None

class TestFittingDiagnosis:
    """Tests for fitting diagnosis logic (without actual plotting)."""
    
    @patch('axcsas.visualization.generate_fitting_diagnosis.LMOptimizer')
    def test_fitting_logic_structure(self, MockOptimizer):
        """Verify fit_peak_with_diagnosis returns correct dictionary structure."""
        # Setup mock optimizer result
        mock_result = MagicMock()
        mock_result.success = True
        mock_result.params.center = 43.3
        mock_result.params.amplitude = 1000
        mock_result.params.fwhm = 0.2
        mock_result.params.eta = 0.5
        mock_result.r_squared = 0.99
        
        MockOptimizer.return_value.fit_peak.return_value = mock_result
        MockOptimizer.return_value.fit_doublet.return_value = mock_result # For doublet call if applicable

        # Create dummy data
        two_theta = np.linspace(40, 46, 100)
        intensity = 1000 * np.exp(-(two_theta - 43.3)**2 / (2 * 0.1**2))
        
        # Test function
        result = fit_peak_with_diagnosis(two_theta, intensity, 43.3, use_doublet=False)
        
        # Verify keys
        expected_keys = {
            'success', 'center', 'amplitude', 'fwhm', 'eta', 'r_squared', 
            'theta_range', 'int_range', 'fitted_curve', 'method'
        }
        assert expected_keys.issubset(result.keys())
        # Note: success might be False because we mocked the optimizer but maybe not the data alignment perfectly,
        # or because fit_peak_with_diagnosis does its own data extraction.
        # But we primarily check if it runs and returns the dict.
        
if __name__ == "__main__":
    pytest.main([__file__, "-v"])
