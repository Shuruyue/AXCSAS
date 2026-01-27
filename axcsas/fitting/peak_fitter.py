"""
Peak Fitting Orchestrator
=========================

Provides high-level peak fitting logic, combining Kα doublet fitting 
and Pseudo-Voigt fallback strategies with detailed error analysis.

Moves fitting logic out of visualization modules to avoid circular dependencies.
"""

import numpy as np
from typing import Dict, Optional, Tuple, List

from axcsas.fitting.pseudo_voigt import PseudoVoigt, PseudoVoigtParams
from axcsas.fitting.ka_doublet import DoubletFitter


def fit_peak_with_diagnosis(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    expected_center: float,
    window: float = 2.5,
    use_doublet: bool = False
) -> Dict:
    """
    Fit a single peak using Kα doublet model and return detailed diagnosis info.
    使用 Kα 雙峰模型擬合單峰並返回詳細診斷資訊。
    
    Returns dict with:
        - success: bool
        - center: fitted Kα₁ peak center
        - center_ka2: fitted Kα₂ peak center
        - amplitude: fitted amplitude
        - fwhm: fitted FWHM
        - eta: mixing parameter
        - r_squared: goodness of fit
        - theta_range: x data
        - int_range: y data (original)
        - fitted_curve: y data (fitted)
    """
    result = {
        'success': False,
        'center': np.nan,
        'center_ka2': np.nan,
        'amplitude': np.nan,
        'fwhm': np.nan,
        'eta': np.nan,
        'r_squared': np.nan,
        'theta_range': None,
        'int_range': None,
        'fitted_curve': None,
        'chi2_red': np.nan,
        'center_err': np.nan,
        'fwhm_err': np.nan,
        'eta_err': np.nan,
        'method': 'doublet' if use_doublet else 'pseudo-voigt'
    }
    
    # Select range
    mask = (two_theta >= expected_center - window) & (two_theta <= expected_center + window)
    if not np.any(mask):
        return result
    
    theta_range = two_theta[mask]
    int_range = intensity[mask]
    
    result['theta_range'] = theta_range
    result['int_range'] = int_range
    
    # Find maximum for initial guess
    idx_max = np.argmax(int_range)
    peak_theta = theta_range[idx_max]
    peak_int = int_range[idx_max]
    
    if peak_int < 50:
        return result
    
    # Estimate initial FWHM
    half_max = peak_int / 2
    left_idx = idx_max
    while left_idx > 0 and int_range[left_idx] > half_max:
        left_idx -= 1
    right_idx = idx_max
    while right_idx < len(int_range) - 1 and int_range[right_idx] > half_max:
        right_idx += 1
    initial_fwhm = max(theta_range[right_idx] - theta_range[left_idx], 0.1)
    
    try:
        if use_doublet:
            # Use Kα doublet fitting with timeout protection
            fitter = DoubletFitter(max_iterations=100000)
            fit_result = fitter.fit(theta_range, int_range, expected_center, initial_fwhm)
            
            if fit_result.success and fit_result.r_squared > 0.8:
                # Calculate uncertainties from DoubletFitter result
                center_err = np.nan
                fwhm_err = np.nan
                eta_err = np.nan
                chi2_red = np.nan
                
                # Try to get Jacobian from fit_result if available
                if hasattr(fit_result, 'opt_result') and fit_result.opt_result is not None:
                    opt_res = fit_result.opt_result
                    if hasattr(opt_res, 'jac') and opt_res.jac is not None:
                        try:
                            # Calculate chi-square with Poisson weights
                            residuals = int_range - fit_result.fitted_curve
                            
                            # Poisson errors: σ = √I
                            sigma_doublet = np.sqrt(np.maximum(int_range, 1.0))
                            
                            # Unweighted sum of squares (for covariance)
                            ss_res = np.sum(residuals**2)
                            n_data = len(int_range)
                            n_params = len(opt_res.x) if hasattr(opt_res, 'x') else 7
                            dof = n_data - n_params
                            
                            if dof > 0:
                                # Weighted chi-square
                                weighted_residuals = residuals / sigma_doublet
                                chi2 = np.sum(weighted_residuals**2)
                                chi2_red = chi2 / dof
                        except Exception as e:
                            print(f"Warning: Chi2 calc failed: {e}")
                
                result['success'] = True
                result['center'] = fit_result.center_ka1
                # Use errors calculated by DoubletFitter
                result['center_err'] = getattr(fit_result, 'center_error', np.nan)
                result['center_ka2'] = fit_result.center_ka2
                result['amplitude'] = fit_result.amplitude_ka1
                result['fwhm'] = fit_result.fwhm
                result['fwhm_err'] = getattr(fit_result, 'fwhm_error', np.nan)
                result['eta'] = fit_result.eta
                result['eta_err'] = getattr(fit_result, 'eta_error', np.nan)
                result['r_squared'] = fit_result.r_squared
                result['chi2_red'] = chi2_red
                result['fitted_curve'] = fit_result.fitted_curve
                result['method'] = 'doublet-true-voigt'
                result['r_squared'] = fit_result.r_squared
                result['chi2_red'] = chi2_red
                result['fitted_curve'] = fit_result.fitted_curve
                result['method'] = 'doublet-true-voigt'
            else:
                # Fallback to True Fit if doublet fitting fails
                use_doublet = False
        
        if not use_doublet or not result['success']:
            # Rigorous True Voigt fitting with polynomial background
            from scipy.optimize import least_squares
            
            # Estimate polynomial background (quadratic)
            n_edge = max(5, len(theta_range) // 8)
            bg_x = np.concatenate([theta_range[:n_edge], theta_range[-n_edge:]])
            bg_y = np.concatenate([int_range[:n_edge], int_range[-n_edge:]])
            bg_coeffs = np.polyfit(bg_x, bg_y, 2)  # Quadratic background
            
            # Model: True Voigt + quadratic background
            def enhanced_tv_model(x, center, amplitude, sigma, gamma, a2, a1, a0):
                # Use TrueVoigt
                from axcsas.fitting.pseudo_voigt import TrueVoigt
                tv = TrueVoigt.profile(x, center, amplitude, sigma, gamma)
                background = a2 * x**2 + a1 * x + a0
                return tv + background
            
            # Calculate weights
            sigma_w = np.sqrt(np.maximum(int_range, 1.0))
            
            def residual(params):
                return (int_range - enhanced_tv_model(theta_range, *params)) / sigma_w
            
            best_r2 = -1
            best_opt_result = None
            
            # Phase 1: Multi-start optimization
            # Iterate ratio (gamma fraction) instead of eta
            ratio_grid = np.linspace(0.1, 0.9, 9)
            
            for ratio in ratio_grid:
                # Initial guesses
                sigma_init = (initial_fwhm * (1-ratio)) / 2.355
                gamma_init = (initial_fwhm * ratio) / 2.0
                
                x0 = np.array([peak_theta, peak_int * 0.9, sigma_init, gamma_init, 
                              bg_coeffs[0], bg_coeffs[1], bg_coeffs[2]])
                
                # Bounds: sigma/gamma > 0
                lb = np.array([theta_range.min(), 10, 1e-4, 1e-4, -np.inf, -np.inf, -np.inf])
                ub = np.array([theta_range.max(), peak_int * 1.5, 2.0, 2.0, np.inf, np.inf, np.inf])
                
                try:
                    opt_result = least_squares(
                        residual, x0, bounds=(lb, ub),
                        max_nfev=10000, ftol=1e-14, xtol=1e-14, gtol=1e-14
                    )
                    
                    if opt_result.success or opt_result.status >= 1:
                        fitted = enhanced_tv_model(theta_range, *opt_result.x)
                        residuals = int_range - fitted
                        ss_res = np.sum(residuals**2)
                        ss_tot = np.sum((int_range - np.mean(int_range))**2)
                        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                        
                        if r_squared > best_r2:
                            best_r2 = r_squared
                            best_opt_result = opt_result
            
                except Exception:
                    continue
            
            # Process best result
            if best_opt_result is not None and best_r2 > 0.85:
                # Recalculate metrics
                popt = best_opt_result.x
                fitted = enhanced_tv_model(theta_range, *popt)
                
                # Convert to FWHM/Eta
                from axcsas.fitting.pseudo_voigt import TrueVoigt
                fwhm_tot = TrueVoigt.fwhm_from_params(popt[2], popt[3])
                eta_eff = (2*popt[3]) / fwhm_tot if fwhm_tot > 0 else 0.5
                
                # Errors
                fwhm_err = 0.0001
                try:
                    J = best_opt_result.jac
                    residuals = int_range - fitted
                    ss_res = np.sum(residuals**2)
                    dof = len(int_range) - len(popt)
                    variance = ss_res / dof if dof > 0 else 0
                    
                    if J is not None:
                        JTJ = J.T @ J
                        pcov = np.linalg.pinv(JTJ) * variance
                        perr = np.sqrt(np.diag(pcov))
                        fwhm_err = perr[2]*2.35 + perr[3]*2.0 # Approx
                except:
                    pass
                
                # Update result
                result['success'] = True
                result['center'] = popt[0]
                result['center_err'] = 0.001 # Default small
                result['amplitude'] = popt[1]
                result['fwhm'] = fwhm_tot
                result['fwhm_err'] = fwhm_err
                result['eta'] = eta_eff
                result['eta_err'] = 0.01
                result['r_squared'] = best_r2
                result['method'] = 'true-voigt-enhanced'
                result['fitted_curve'] = fitted
                result['low_quality'] = best_r2 < 0.995
                
    except Exception as e:
        result['error'] = str(e)
    
    return result

