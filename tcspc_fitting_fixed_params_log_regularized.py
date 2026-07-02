import os
from datetime import datetime
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares, nnls
from sdtfile import SdtFile
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import filedialog, messagebox

# Preview of data then select exponentials, see if there's a one size fits all exponential based on a numerical input
# Giving one exponential or fit until chi square is below 2
# Give guesses for exponential values for fitting 
# Background: ignore the baselines 
# Shows the life time and amplitude ratios in the plot, and chi squared as a part of the plot 
# fancy error functions, A1 A2, amplitutde and life time errors 
# Amplitude ratios, what is the intensity ratio between A1, A2, and A3 

# https://github.com/PicoQuant/snAPI
# https://rna-fretools.github.io/Lifefit/

# Work until one bin in the histogram fully saturates 

# SDT loading


def load_sdt_decay(path, block=0):
    """
    Load decay and time axis from an SDT file.

    sdtfile often returns time in seconds, while FluoFit reports use ns.
    This function converts seconds -> nanoseconds when needed.
    """
    with SdtFile(path) as sdt:
        data = np.asarray(sdt.data[block], dtype=float)
        time = np.asarray(sdt.times[block], dtype=float)

    # Collapse non-time axes into one decay.
    if data.ndim == 1:
        decay = data
    else:
        decay = data.reshape(-1, data.shape[-1]).sum(axis=0)

    # Make sure time is 1D.
    time = np.squeeze(time)

    # If time is in seconds, convert to ns.
    # A 50 ns window would appear as 5e-8 seconds.
    if np.nanmax(time) < 1e-3:
        time_ns = time * 1e9
    else:
        time_ns = time

    return time_ns, decay



# Model functions


def shift_curve(y, shift_bins):
    """
    Fractionally shift a curve using linear interpolation.

    Positive shift_bins shifts the curve to the right.
    Negative shift_bins shifts the curve to the left.
    """
    y = np.asarray(y, dtype=float)
    x = np.arange(len(y), dtype=float)
    return np.interp(x - shift_bins, x, y, left=0.0, right=0.0)


def unpack_multi_exp_params(params, n_exponentials):
    """
    Backward-compatible parser for the GUI's full parameter list.

    Parameter order:
        A1, tau1_ns, A2, tau2_ns, ..., An, taun_ns,
        shift_irf_ns, background_decay, background_irf

    In the LifeFit-style fit below, only tau values and shift_irf_ns are used as
    nonlinear variables. Amplitudes and the constant decay background are solved
    by scipy.optimize.nnls at every optimizer step.
    """
    params = np.asarray(params, dtype=float)
    amplitudes = params[0 : 2 * n_exponentials : 2]
    taus_ns = params[1 : 2 * n_exponentials : 2]
    shift_irf_ns = params[2 * n_exponentials]
    background_decay = params[2 * n_exponentials + 1]
    background_irf = params[2 * n_exponentials + 2]
    return amplitudes, taus_ns, shift_irf_ns, background_decay, background_irf


def full_params_to_lifefit_params(params, n_exponentials):
    """
    Convert GUI-style initial parameters to LifeFit-style nonlinear parameters.

    Returned order:
        tau1_ns, tau2_ns, ..., taun_ns, shift_irf_ns
    """
    _, taus_ns, shift_irf_ns, _, _ = unpack_multi_exp_params(params, n_exponentials)
    return np.r_[taus_ns, shift_irf_ns].astype(float)


def unpack_lifefit_params(params, n_exponentials):
    """
    LifeFit-style nonlinear parameter order:
        tau1_ns, tau2_ns, ..., taun_ns, shift_irf_ns
    """
    params = np.asarray(params, dtype=float)
    taus_ns = params[:n_exponentials]
    shift_irf_ns = float(params[n_exponentials])
    return taus_ns, shift_irf_ns


def lifefit_bounds_from_full_bounds(n_exponentials):
    """
    Bounds for the LifeFit-style nonlinear variables.
    Amplitudes and background do not need bounds because NNLS enforces x >= 0.
    """
    lower = [0.01] * n_exponentials + [-5.0]
    upper = [100.0] * n_exponentials + [5.0]
    return np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)




def full_fixed_mask_to_lifefit_fixed_mask(fixed_params, n_exponentials):
    """
    Convert a full GUI-style fixed-parameter mask into the LifeFit nonlinear mask.

    Full GUI parameter order:
        A1, tau1, A2, tau2, ..., An, taun, shift_irf, background_decay, background_irf

    LifeFit nonlinear order:
        tau1, tau2, ..., taun, shift_irf

    True means that parameter is fixed at the user's initial value.
    """
    if fixed_params is None:
        return np.zeros(n_exponentials + 1, dtype=bool)

    fixed_params = np.asarray(fixed_params, dtype=bool)
    expected_param_count = 2 * n_exponentials + 3
    if fixed_params.size != expected_param_count:
        raise ValueError(
            f"Expected {expected_param_count} fixed-parameter flags for "
            f"{n_exponentials} exponentials, but received {fixed_params.size}."
        )

    tau_fixed = fixed_params[1 : 2 * n_exponentials : 2]
    shift_fixed = fixed_params[2 * n_exponentials]
    return np.r_[tau_fixed, shift_fixed].astype(bool)


def expand_free_lifefit_params(free_params, fixed_values, fixed_mask):
    """
    Rebuild the full LifeFit nonlinear vector from free optimizer variables.
    """
    fixed_values = np.asarray(fixed_values, dtype=float)
    fixed_mask = np.asarray(fixed_mask, dtype=bool)
    full = fixed_values.copy()
    full[~fixed_mask] = np.asarray(free_params, dtype=float)
    return full


def lifefit_full_to_free_params(full_params, fixed_mask):
    """
    Extract only the non-fixed nonlinear parameters for scipy.optimize.least_squares.
    """
    return np.asarray(full_params, dtype=float)[~np.asarray(fixed_mask, dtype=bool)]


def lifefit_free_bounds(lower, upper, fixed_mask):
    """
    Extract bounds for the non-fixed nonlinear parameters.
    """
    fixed_mask = np.asarray(fixed_mask, dtype=bool)
    return np.asarray(lower, dtype=float)[~fixed_mask], np.asarray(upper, dtype=float)[~fixed_mask]

def single_exponential(time_ns, tau_ns):
    """LifeFit-style single exponential: exp(-t / tau)."""
    time_ns = np.asarray(time_ns, dtype=float)
    tau_ns = max(float(tau_ns), np.finfo(float).tiny)
    return np.exp(-time_ns / tau_ns)


def fft_convolution(irf, sgl_exp):
    """
    FFT convolution used by LifeFit's reconvolution workflow.

    The result has the same length as the input arrays. This circular FFT form
    follows the package's documented ifft(fft(exp) * fft(irf)) calculation.
    """
    irf = np.asarray(irf, dtype=float)
    sgl_exp = np.asarray(sgl_exp, dtype=float)
    return np.real(np.fft.ifft(np.fft.fft(sgl_exp) * np.fft.fft(irf)))


def build_lifefit_design_matrix(time_ns, irf, taus_ns, shift_irf_ns, dt):
    """
    Build LifeFit's NNLS design matrix.

    Columns 0..n-1 are IRF-convolved single-exponential decays.
    Final column is ones for the constant decay background.
    """
    time_ns = np.asarray(time_ns, dtype=float)
    irf = np.asarray(irf, dtype=float)
    taus_ns = np.asarray(taus_ns, dtype=float)

    shifted_irf = shift_curve(irf, shift_irf_ns / dt)
    t = time_ns - time_ns[0]

    columns = []
    for tau_ns in taus_ns:
        sgl_exp = single_exponential(t, tau_ns)
        columns.append(fft_convolution(shifted_irf, sgl_exp))

    columns.append(np.ones_like(time_ns, dtype=float))
    return np.column_stack(columns), shifted_irf


def solve_lifefit_nnls(
    time_ns,
    decay,
    irf,
    taus_ns,
    shift_irf_ns,
    dt,
    fit_mask=None,
    sigma=None,
    fixed_amplitude_mask=None,
    fixed_amplitude_values=None,
    fixed_background_decay=False,
    fixed_background_decay_value=0.0,
):
    """
    Solve LifeFit's linear subproblem for fixed lifetimes and IRF shift.

    Amplitudes and the constant decay background are normally solved by NNLS.
    If the user checks "Fix" for an amplitude or background decay, that value is
    held at the initial estimate and the NNLS fit is done on the remaining linear
    coefficients only.

    Returns:
        model, amplitudes, background_decay, shifted_irf, nnls_residual_norm,
        variable_linear_count
    """
    decay = np.asarray(decay, dtype=float)
    taus_ns = np.asarray(taus_ns, dtype=float)
    n_exponentials = len(taus_ns)

    A_full, shifted_irf = build_lifefit_design_matrix(
        time_ns=time_ns,
        irf=irf,
        taus_ns=taus_ns,
        shift_irf_ns=shift_irf_ns,
        dt=dt,
    )

    if fit_mask is None:
        fit_mask = np.ones_like(decay, dtype=bool)

    if fixed_amplitude_mask is None:
        fixed_amplitude_mask = np.zeros(n_exponentials, dtype=bool)
    else:
        fixed_amplitude_mask = np.asarray(fixed_amplitude_mask, dtype=bool)

    if fixed_amplitude_values is None:
        fixed_amplitude_values = np.zeros(n_exponentials, dtype=float)
    else:
        fixed_amplitude_values = np.asarray(fixed_amplitude_values, dtype=float)

    if fixed_amplitude_mask.size != n_exponentials:
        raise ValueError("fixed_amplitude_mask must have one entry per exponential.")
    if fixed_amplitude_values.size != n_exponentials:
        raise ValueError("fixed_amplitude_values must have one entry per exponential.")

    fixed_coeffs = np.zeros(n_exponentials + 1, dtype=float)
    linear_fixed_mask = np.zeros(n_exponentials + 1, dtype=bool)

    linear_fixed_mask[:n_exponentials] = fixed_amplitude_mask
    fixed_coeffs[:n_exponentials] = np.where(fixed_amplitude_mask, fixed_amplitude_values, 0.0)

    linear_fixed_mask[-1] = bool(fixed_background_decay)
    fixed_coeffs[-1] = float(fixed_background_decay_value) if fixed_background_decay else 0.0

    fixed_contribution_full = A_full @ fixed_coeffs
    y_remaining_full = decay - fixed_contribution_full

    A_fit_all = A_full[fit_mask]
    y_fit = y_remaining_full[fit_mask]
    free_linear_mask = ~linear_fixed_mask
    A_fit = A_fit_all[:, free_linear_mask]

    if sigma is not None:
        sigma_fit = np.asarray(sigma, dtype=float)[fit_mask]
        sigma_fit = np.maximum(sigma_fit, np.finfo(float).eps)
        A_fit = A_fit / sigma_fit[:, None]
        y_fit = y_fit / sigma_fit

    coeffs = fixed_coeffs.copy()
    if A_fit.shape[1] > 0:
        free_coeffs, nnls_residual_norm = nnls(A_fit, y_fit)
        coeffs[free_linear_mask] = free_coeffs
    else:
        nnls_residual_norm = float(np.linalg.norm(y_fit))

    model = A_full @ coeffs
    amplitudes = coeffs[:-1]
    background_decay = float(coeffs[-1])
    variable_linear_count = int(np.sum(free_linear_mask))

    return model, amplitudes, background_decay, shifted_irf, nnls_residual_norm, variable_linear_count

def reconvolved_multi_exp(params, time_ns, irf, dt, n_exponentials):
    """
    Backward-compatible reconvolution model using supplied amplitudes/background.

    The main fit path no longer optimizes these amplitudes nonlinearly; it uses
    solve_lifefit_nnls instead. This helper remains available for plotting or
    comparing a manually supplied full parameter vector.
    """
    amplitudes, taus_ns, shift_irf_ns, background_decay, background_irf = (
        unpack_multi_exp_params(params, n_exponentials)
    )
    A_full, _ = build_lifefit_design_matrix(
        time_ns=time_ns,
        irf=irf + background_irf,
        taus_ns=taus_ns,
        shift_irf_ns=shift_irf_ns,
        dt=dt,
    )
    coeffs = np.r_[amplitudes, background_decay]
    return A_full @ coeffs


def lifefit_residuals(
    params,
    time_ns,
    decay,
    irf,
    dt,
    fit_mask,
    n_exponentials,
    sigma=None,
    fixed_amplitude_mask=None,
    fixed_amplitude_values=None,
    fixed_background_decay=False,
    fixed_background_decay_value=0.0,
):
    taus_ns, shift_irf_ns = unpack_lifefit_params(params, n_exponentials)
    model, _, _, _, _, _ = solve_lifefit_nnls(
        time_ns=time_ns,
        decay=decay,
        irf=irf,
        taus_ns=taus_ns,
        shift_irf_ns=shift_irf_ns,
        dt=dt,
        fit_mask=fit_mask,
        sigma=sigma,
        fixed_amplitude_mask=fixed_amplitude_mask,
        fixed_amplitude_values=fixed_amplitude_values,
        fixed_background_decay=fixed_background_decay,
        fixed_background_decay_value=fixed_background_decay_value,
    )

    residual = decay[fit_mask] - model[fit_mask]
    if sigma is not None:
        sigma_fit = np.maximum(np.asarray(sigma, dtype=float)[fit_mask], np.finfo(float).eps)
        residual = residual / sigma_fit
    return residual


def lifefit_free_residuals(
    free_params,
    fixed_nonlinear_values,
    fixed_nonlinear_mask,
    time_ns,
    decay,
    irf,
    dt,
    fit_mask,
    n_exponentials,
    sigma=None,
    fixed_amplitude_mask=None,
    fixed_amplitude_values=None,
    fixed_background_decay=False,
    fixed_background_decay_value=0.0,
):
    full_params = expand_free_lifefit_params(
        free_params=free_params,
        fixed_values=fixed_nonlinear_values,
        fixed_mask=fixed_nonlinear_mask,
    )
    return lifefit_residuals(
        params=full_params,
        time_ns=time_ns,
        decay=decay,
        irf=irf,
        dt=dt,
        fit_mask=fit_mask,
        n_exponentials=n_exponentials,
        sigma=sigma,
        fixed_amplitude_mask=fixed_amplitude_mask,
        fixed_amplitude_values=fixed_amplitude_values,
        fixed_background_decay=fixed_background_decay,
        fixed_background_decay_value=fixed_background_decay_value,
    )


def poisson_residuals(params, time_ns, decay, irf, dt, fit_mask, n_exponentials):
    """
    Compatibility wrapper: LifeFit-style residuals with Poisson sigma.
    """
    sigma = np.sqrt(np.maximum(decay, 1.0))
    return lifefit_residuals(
        params=full_params_to_lifefit_params(params, n_exponentials),
        time_ns=time_ns,
        decay=decay,
        irf=irf,
        dt=dt,
        fit_mask=fit_mask,
        n_exponentials=n_exponentials,
        sigma=sigma,
    )


# Fitting workflow


def make_initial_params(decay, n_exponentials):
    """
    Create reasonable generic initial guesses for any number of exponentials.
    """
    max_decay = float(np.max(decay))
    amplitude_guess = max_decay / max(n_exponentials, 1)

    # Spread lifetime guesses from fast to slow. These are generic and can be
    # adjusted for a specific instrument/sample if desired.
    if n_exponentials == 1:
        tau_guesses = np.array([5.5], dtype=float)
    else:
        tau_guesses = np.geomspace(0.5, 12.0, n_exponentials)

    params = []
    for tau_guess in tau_guesses:
        params.extend([amplitude_guess, float(tau_guess)])

    params.extend([
        -0.04,  # IRF shift, ns
        0.0,    # decay background
        0.0,    # IRF background
    ])

    return np.asarray(params, dtype=float)


def make_bounds(n_exponentials):
    lower_bounds = []
    upper_bounds = []

    for _ in range(n_exponentials):
        lower_bounds.extend([0.0, 0.01])
        upper_bounds.extend([np.inf, 100.0])

    lower_bounds.extend([-5.0, 0.0, 0.0])
    upper_bounds.extend([5.0, np.inf, np.inf])

    return lower_bounds, upper_bounds


def normalize_irf_max_to_decay_by_offset(irf, decay):
    """
    Scale-normalize the IRF so its maximum equals the decay maximum.

    This uses proportional scaling, not flat offset subtraction:
        normalized_irf = irf * (max(decay) / max(irf))

    Returns:
        normalized_irf, scale_factor
    """
    irf = np.asarray(irf, dtype=float)
    decay = np.asarray(decay, dtype=float)

    irf_max = float(np.max(irf))
    decay_max = float(np.max(decay))

    if irf_max <= 0:
        raise ValueError("Cannot normalize IRF because its maximum is zero or negative.")

    scale_factor = decay_max / irf_max
    normalized_irf = irf * scale_factor

    return normalized_irf, scale_factor




def regularize_decay_min_one(decay):
    """
    Optional decay regularization to set min to 1 or 10^0 to match tcspc outputs

    """
    decay = np.asarray(decay, dtype=float)
    decay_min = float(np.min(decay))
    print(decay_min)
    
    # Safety check: Prevent division by zero if the minimum is exactly 0
    if decay_min == 0:
        raise ValueError("The minimum decay value is 0. Cannot scale by a ratio to 1.")
        
    # Apply the ratio to all values
    regularized_decay = decay / decay_min
    
    info = {
        "enabled": True,
        "original_min": decay_min,
        "original_max": float(np.max(decay)),
        "regularized_min": float(np.min(regularized_decay)),
        "regularized_max": float(np.max(regularized_decay)),
        "method": "divided by original minimum (ratio scaling)",
    }
    return regularized_decay, info

def fit_multi_exp_reconv(
    decay_sdt_path,
    irf_sdt_path,
    n_exponentials=1,
    decay_block=0,
    irf_block=0,
    fit_start_ns=5,
    fit_end_ns=48,
    initial_params=None,
    fixed_params=None,
    use_poisson_weights=False,
):
    n_exponentials = int(n_exponentials)
    if n_exponentials < 1:
        raise ValueError("Number of exponentials must be at least 1.")

    time_ns, decay = load_sdt_decay(decay_sdt_path, decay_block)
    irf_time_ns, irf = load_sdt_decay(irf_sdt_path, irf_block)
    print("Decay time first 10:", time_ns[:10])
    print("Decay time last:", time_ns[-1])
    print("Decay dt:", np.median(np.diff(time_ns)))

    print("IRF time first 10:", irf_time_ns[:10])
    print("IRF time last:", irf_time_ns[-1])
    print("IRF dt:", np.median(np.diff(irf_time_ns)))

    # Interpolate IRF onto decay time axis if needed.
    if len(irf_time_ns) != len(time_ns) or not np.allclose(irf_time_ns, time_ns):
        irf = np.interp(time_ns, irf_time_ns, irf, left=0.0, right=0.0)

    decay_raw = decay.copy()
    decay_regularization_info = {"enabled": False}

    # --------------------------------
    # Comment out the next line to fit the original raw decay instead.
    #------------------------------

    #decay, decay_regularization_info = regularize_decay_min_one(decay)SS


    # Scale-normalize the loaded/interpolated IRF so its maximum equals the
    # current decay maximum using proportional scaling.
    irf, irf_scale_factor = normalize_irf_max_to_decay_by_offset(irf, decay)
    irf = np.maximum(irf, 0.0)
    print("IRF scale factor to match decay max:", irf_scale_factor)
    print("Decay max after loading:", np.max(decay))
    print("IRF max after offset normalization:", np.max(irf))

    dt = float(np.median(np.diff(time_ns)))

    if fit_start_ns is None:
        # Better dynamic masking
        peak_idx = np.argmax(decay)
        # Find the index where the rising edge crosses 5% of the peak
        threshold = 0.05 * decay[peak_idx]
        start_idx = np.where(decay[:peak_idx] > threshold)[0][0] 
        fit_start_ns = float(time_ns[start_idx])
        #fit_start_ns = float(time_ns[np.argmax(decay)])
    if fit_end_ns is None:
        fit_end_ns = float(time_ns[-1])

    fit_mask = (time_ns >= fit_start_ns) & (time_ns <= fit_end_ns)

    if initial_params is None:
        p0 = make_initial_params(decay, n_exponentials)
    else:
        p0 = np.asarray(initial_params, dtype=float)
        expected_param_count = 2 * n_exponentials + 3
        if p0.size != expected_param_count:
            raise ValueError(
                f"Expected {expected_param_count} initial parameters for "
                f"{n_exponentials} exponentials, but received {p0.size}."
            )

    if fixed_params is None:
        fixed_params = np.zeros_like(p0, dtype=bool)
    else:
        fixed_params = np.asarray(fixed_params, dtype=bool)
        if fixed_params.size != p0.size:
            raise ValueError(
                f"Expected {p0.size} fixed-parameter flags, but received {fixed_params.size}."
            )

    # Validate the full GUI-style initial vector, then convert it to LifeFit's
    # nonlinear vector: tau1..taun and IRF shift. Amplitudes/background are not
    # nonlinear fit variables; they are re-solved by NNLS unless checked as fixed.
    lower_bounds, upper_bounds = make_bounds(n_exponentials)
    lower_bounds_arr = np.asarray(lower_bounds, dtype=float)
    upper_bounds_arr = np.asarray(upper_bounds, dtype=float)
    if np.any(p0 < lower_bounds_arr) or np.any(p0 > upper_bounds_arr):
        raise ValueError(
            "One or more initial estimates are outside the allowed bounds. "
            "Use positive amplitudes, lifetimes between 0.01 and 100 ns, "
            "IRF shift between -5 and 5 ns, and non-negative backgrounds."
        )

    nonlinear_p0 = full_params_to_lifefit_params(p0, n_exponentials)
    nonlinear_lower, nonlinear_upper = lifefit_bounds_from_full_bounds(n_exponentials)
    fixed_nonlinear_mask = full_fixed_mask_to_lifefit_fixed_mask(fixed_params, n_exponentials)

    initial_amplitudes, _, _, initial_bg_decay, initial_bg_irf = unpack_multi_exp_params(
        p0, n_exponentials
    )
    fixed_amplitude_mask = np.asarray(fixed_params[0 : 2 * n_exponentials : 2], dtype=bool)
    fixed_background_decay = bool(fixed_params[2 * n_exponentials + 1])
    fixed_background_irf = bool(fixed_params[2 * n_exponentials + 2])
    fixed_background_irf_value = float(initial_bg_irf) if fixed_background_irf else 0.0

    # Background IRF is not independently identifiable in the NNLS subproblem here.
    # When checked as fixed, it is added to the IRF before reconvolution. When it
    # is not checked, the existing behavior is preserved: background_irf = 0.
    irf_for_fit = irf + fixed_background_irf_value

    free_nonlinear_p0 = lifefit_full_to_free_params(nonlinear_p0, fixed_nonlinear_mask)
    free_lower, free_upper = lifefit_free_bounds(
        nonlinear_lower, nonlinear_upper, fixed_nonlinear_mask
    )

    sigma = np.sqrt(np.maximum(decay, 1.0)) if use_poisson_weights else None

    if free_nonlinear_p0.size > 0:
        result = least_squares(
            lifefit_free_residuals,
            free_nonlinear_p0,
            bounds=(free_lower, free_upper),
            args=(
                nonlinear_p0,
                fixed_nonlinear_mask,
                time_ns,
                decay,
                irf_for_fit,
                dt,
                fit_mask,
                n_exponentials,
                sigma,
                fixed_amplitude_mask,
                initial_amplitudes,
                fixed_background_decay,
                initial_bg_decay,
            ),
            loss="linear",
            max_nfev=20000,
            ftol=1e-8,
            xtol=1e-8,
            gtol=1e-8,
        )
        nonlinear_params = expand_free_lifefit_params(
            free_params=result.x,
            fixed_values=nonlinear_p0,
            fixed_mask=fixed_nonlinear_mask,
        )
    else:
        nonlinear_params = nonlinear_p0.copy()
        result = SimpleNamespace(
            x=np.asarray([], dtype=float),
            success=True,
            message="All nonlinear parameters were fixed; skipped least_squares.",
        )
    taus_ns, shift_irf_ns = unpack_lifefit_params(nonlinear_params, n_exponentials)

    model, amplitudes, bg_decay, shifted_irf, nnls_residual_norm, variable_linear_count = solve_lifefit_nnls(
        time_ns=time_ns,
        decay=decay,
        irf=irf_for_fit,
        taus_ns=taus_ns,
        shift_irf_ns=shift_irf_ns,
        dt=dt,
        fit_mask=fit_mask,
        sigma=sigma,
        fixed_amplitude_mask=fixed_amplitude_mask,
        fixed_amplitude_values=initial_amplitudes,
        fixed_background_decay=fixed_background_decay,
        fixed_background_decay_value=initial_bg_decay,
    )
    bg_irf = fixed_background_irf_value

    full_mask = np.ones_like(time_ns, dtype=bool)
    residuals = lifefit_residuals(
        nonlinear_params,
        time_ns,
        decay,
        irf_for_fit,
        dt,
        full_mask,
        n_exponentials,
        sigma,
        fixed_amplitude_mask,
        initial_amplitudes,
        fixed_background_decay,
        initial_bg_decay,
    )

    fit_residuals = lifefit_residuals(
        nonlinear_params,
        time_ns,
        decay,
        irf_for_fit,
        dt,
        fit_mask,
        n_exponentials,
        sigma,
        fixed_amplitude_mask,
        initial_amplitudes,
        fixed_background_decay,
        initial_bg_decay,
    )
    n_fit_points = int(fit_mask.sum())
    # Degrees of freedom count only free fitted variables, not fixed user values.
    n_params = int(np.sum(~fixed_nonlinear_mask)) + int(variable_linear_count)

    reduced_chi2 = np.sum(fit_residuals ** 2) / max(n_fit_points - n_params, 1)

    amplitude_sum = float(np.sum(amplitudes))
    amplitude_fractions = (
        amplitudes / amplitude_sum if amplitude_sum > 0 else np.zeros_like(amplitudes)
    )

    intensity_weights_raw = amplitudes * taus_ns
    intensity_sum = float(np.sum(intensity_weights_raw))
    intensity_fractions = (
        intensity_weights_raw / intensity_sum
        if intensity_sum > 0
        else np.zeros_like(intensity_weights_raw)
    )

    amplitude_weighted_lifetime = (
        float(np.sum(amplitudes * taus_ns) / amplitude_sum) if amplitude_sum > 0 else np.nan
    )
    intensity_weighted_lifetime = (
        float(np.sum(amplitudes * taus_ns ** 2) / np.sum(amplitudes * taus_ns))
        if np.sum(amplitudes * taus_ns) > 0
        else np.nan
    )

    return {
        "time_ns": time_ns,
        "decay": decay,
        "decay_raw": decay_raw,
        "decay_regularization_info": decay_regularization_info,
        "irf": irf_for_fit,
        "shifted_irf": shifted_irf,
        #"irf_offset_subtracted": irf_offset_subtracted,
        "fit": model,
        "residuals": residuals,
        "fit_mask": fit_mask,
        "n_exponentials": n_exponentials,
        "initial_params": p0.copy(),
        "nonlinear_params": nonlinear_params.copy(),
        "fixed_params": fixed_params.copy(),
        "fixed_nonlinear_mask": fixed_nonlinear_mask.copy(),
        "nnls_residual_norm": nnls_residual_norm,
        "fit_method": "LifeFit-style NNLS reconvolution",
        "use_poisson_weights": use_poisson_weights,
        "amplitudes_counts": amplitudes,
        "taus_ns": taus_ns,
        "amplitude_fractions": amplitude_fractions,
        "intensity_fractions": intensity_fractions,
        "amplitude_weighted_lifetime_ns": amplitude_weighted_lifetime,
        "intensity_weighted_lifetime_ns": intensity_weighted_lifetime,
        "shift_irf_ns": shift_irf_ns,
        "background_decay_counts": bg_decay,
        "background_irf_counts": bg_irf,
        "reduced_chi2": reduced_chi2,
        "fitted_points": n_fit_points,
        "success": result.success,
        "fit_message": result.message,
    }


# DAT writer

def write_fluofit_style_dat(
    fit_result,
    output_dat_path,
    decay_sdt_path,
    irf_sdt_path,
    title=None,
):
    time_ns = np.asarray(fit_result["time_ns"])
    decay = np.asarray(fit_result["decay"])
    irf = np.asarray(fit_result["irf"])
    model = np.asarray(fit_result["fit"])
    residuals = np.asarray(fit_result["residuals"])

    n_exponentials = int(fit_result["n_exponentials"])
    amplitudes = np.asarray(fit_result["amplitudes_counts"])
    taus_ns = np.asarray(fit_result["taus_ns"])
    amplitude_fractions = np.asarray(fit_result["amplitude_fractions"])
    intensity_fractions = np.asarray(fit_result["intensity_fractions"])
    shift_irf = fit_result["shift_irf_ns"]
    bg_decay = fit_result["background_decay_counts"]
    bg_irf = fit_result["background_irf_counts"]
    reduced_chi2 = fit_result["reduced_chi2"]
    fitted_points = fit_result["fitted_points"]

    decay_name = os.path.basename(decay_sdt_path)
    irf_name = os.path.basename(irf_sdt_path)

    if title is None:
        title = f"Python LifeFit-style {n_exponentials}-exponential reconvolution fit"

    now = datetime.now().strftime("%m/%d/%Y %I:%M:%S %p")

    lines = []

    lines.append("PicoQuant FluoFit-style Python Export")
    lines.append(f"Saved : {now}")
    lines.append("")
    lines.append(title)
    lines.append("")
    lines.append(f"Model: Exp. [Reconv.] NNLS ({n_exponentials} Exponential Components)")
    lines.append("Plotted Data Set #0 Decay:")
    lines.append(f'"{decay_name}" (0)')
    lines.append("Plotted Data Set #0 IRF:")
    lines.append(f'"{irf_name}" (0)')
    lines.append(f"X²(reduced): {reduced_chi2:.4f} ; Fitted Data Points: {fitted_points}")
    lines.append("")

    lines.append("Data Set #0")
    lines.append(
        "     Decay                 IRF                   Model Decay           Residuals"
    )
    lines.append(
        "     t[ns]     Intens.     t[ns]     Intens.     t[ns]     Intens.     t[ns]       diff."
    )

    for t, y, h, m, r in zip(time_ns, decay, irf, model, residuals):
        lines.append(
            f"{t:12.6f}"
            f"{y:12.0f}"
            f"{t:12.6f}"
            f"{h:12.0f}"
            f"{t:12.6f}"
            f"{m:12.4f}"
            f"{t:12.6f}"
            f"{r:12.4f}"
        )

    lines.append("")
    lines.append("Best Fit Parameters:")
    lines.append("")
    lines.append("Data Set #0")
    lines.append("Parameter              Value          Conf. Lower     Conf. Upper     Conf. Estimation")

    for idx, (amplitude, tau_ns) in enumerate(zip(amplitudes, taus_ns), start=1):
        lines.append(
            f"A{idx} [Cnts]        {amplitude:14.5f}              ---            ---     Fitting"
        )
        lines.append(
            f"t{idx} [ns]          {tau_ns:14.5f}              ---            ---     Fitting"
        )

    lines.append(
        f"Bkgr. Dec [Cnts] {bg_decay:14.5f}              ---            ---     <none>"
    )
    lines.append(
        f"Bkgr. IRF [Cnts] {bg_irf:14.5f}              ---            ---     <none>"
    )
    lines.append(
        f"Shift IRF [ns]   {shift_irf:14.5f}              ---            ---     Fitting"
    )

    lines.append("")
    lines.append("Average Lifetime:")
    lines.append(
        f"tAv.1={fit_result['intensity_weighted_lifetime_ns']:.5f} ns (intensity weighted)"
    )
    lines.append(
        f"tAv.2={fit_result['amplitude_weighted_lifetime_ns']:.5f} ns (amplitude weighted)"
    )
    lines.append("")

    lines.append("Fractional Intensities of the Positive Decay Components:")
    for idx, (tau_ns, fraction) in enumerate(zip(taus_ns, intensity_fractions), start=1):
        lines.append(f"t{idx} ({tau_ns:.5f} ns) : {100.0 * fraction:.2f}%")
    lines.append("")

    lines.append("Fractional Amplitudes of the Positive Decay Components:")
    for idx, (tau_ns, fraction) in enumerate(zip(taus_ns, amplitude_fractions), start=1):
        lines.append(f"t{idx} ({tau_ns:.5f} ns) : {100.0 * fraction:.2f}%")
    lines.append("")

    lines.append("Fitted Decay and Exponential Components:")
    component_header = "     Time [ns]" + "".join(
        f"      t{idx} [kCounts]" for idx in range(1, n_exponentials + 1)
    ) + "      Sum [kCounts]"
    lines.append(component_header)

    component_time = np.linspace(0.0, 22.0, 1000)
    components = [amplitude * np.exp(-component_time / tau_ns) for amplitude, tau_ns in zip(amplitudes, taus_ns)]
    component_sum = np.sum(components, axis=0)

    for row_idx, t in enumerate(component_time):
        row = f"{t:12.5f}"
        for component in components:
            row += f"{component[row_idx] / 1000.0:16.8f}"
        row += f"{component_sum[row_idx] / 1000.0:16.8f}"
        lines.append(row)

    lines.append("")
    lines.append("Confidence Intervals:")
    for idx in range(1, n_exponentials + 1):
        lines.append(f"A{idx} : [--- ; ---] Cnts")
        lines.append(f"t{idx} : [--- ; ---] ns")
    lines.append("Shift IRF : [--- ; ---] ns")
    lines.append("")

    text = "\r\n".join(lines)

    with open(output_dat_path, "w", encoding="latin-1", newline="") as f:
        f.write(text)

    return output_dat_path


# Plotting

def plot_fit_and_residuals(fit_result):
    time_ns = fit_result["time_ns"]
    decay = fit_result["decay"]
    irf = fit_result["irf"]
    model = fit_result["fit"]
    residuals = fit_result["residuals"]
    n_exponentials = fit_result["n_exponentials"]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )

    ax1.semilogy(
        time_ns, np.maximum(decay, 1), ".", markersize=2, label="Decay"
    )
    ax1.semilogy(
        time_ns, np.maximum(model, 1), "-", linewidth=1.2, label="Model fit"
    )
    ax1.semilogy(time_ns, np.maximum(irf, 1), "-", linewidth=1.0, label="IRF")
    ax1.set_ylabel("Intensity [counts]")
    ax1.set_title(
        f"{n_exponentials}-Exponential Decay Regression Fit\n"
        f"Reduced chi-square = {fit_result['reduced_chi2']:.5f}"
    )
    ax1.legend()

    param_lines = []
    for idx, (amplitude, tau_ns, amp_frac, int_frac) in enumerate(
        zip(
            fit_result["amplitudes_counts"],
            fit_result["taus_ns"],
            fit_result["amplitude_fractions"],
            fit_result["intensity_fractions"],
        ),
        start=1,
    ):
        param_lines.append(
            f"A{idx}={amplitude:.3g}, tau{idx}={tau_ns:.4g} ns, "
            f"Amp={100.0 * amp_frac:.1f}%, Int={100.0 * int_frac:.1f}%"
        )

    ax1.text(
        0.98,
        0.98,
        "\n".join(param_lines),
        transform=ax1.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )

    ax2.plot(time_ns, residuals, linewidth=0.8)
    ax2.axhline(0, linestyle="--", linewidth=0.8)
    ax2.set_xlabel("Time [ns]")
    ax2.set_ylabel("Residual")

    plt.tight_layout()
    plt.show()



# Tkinter dialogs


def ask_initial_parameter_estimates(parent, n_exponentials, default_params):
    """
    Ask the user for initial estimates before running the nonlinear fit.

    Parameter order:
        A1, tau1, A2, tau2, ..., An, taun,
        shift_irf_ns, background_decay, background_irf

    Each row has a Fix checkbox. If checked, that parameter stays at the entered
    value and the fit is done on the remaining variables.
    """
    dialog = tk.Toplevel(parent)
    dialog.title("Initial Fit Parameter Estimates")
    dialog.resizable(False, False)
    dialog.grab_set()

    default_params = np.asarray(default_params, dtype=float)
    result = {"ok": False, "initial_params": None, "fixed_params": None}
    entries = []

    intro = (
        "Enter initial guesses for the fit. Check Fix to keep a value constant; "
        "unchecked parameters are fitted or NNLS-solved."
    )
    tk.Label(dialog, text=intro, wraplength=520, justify="left").grid(
        row=0, column=0, columnspan=4, padx=12, pady=(12, 8), sticky="w"
    )

    row = 1
    tk.Label(dialog, text="Parameter", font=("TkDefaultFont", 9, "bold")).grid(
        row=row, column=0, padx=12, pady=(0, 4), sticky="w"
    )
    tk.Label(dialog, text="Initial estimate", font=("TkDefaultFont", 9, "bold")).grid(
        row=row, column=1, padx=12, pady=(0, 4), sticky="w"
    )
    tk.Label(dialog, text="Fix", font=("TkDefaultFont", 9, "bold")).grid(
        row=row, column=2, padx=12, pady=(0, 4), sticky="w"
    )
    tk.Label(dialog, text="Allowed range", font=("TkDefaultFont", 9, "bold")).grid(
        row=row, column=3, padx=12, pady=(0, 4), sticky="w"
    )

    labels = []
    allowed_ranges = []
    for idx in range(1, n_exponentials + 1):
        labels.extend([f"A{idx} [Cnts]", f"tau{idx} [ns]"])
        allowed_ranges.extend([">= 0", "0.01 to 100"])
    labels.extend(["Shift IRF [ns]", "Background decay [Cnts]", "Background IRF [Cnts]"])
    allowed_ranges.extend(["-5 to 5", ">= 0", ">= 0"])

    for label, value, allowed in zip(labels, default_params, allowed_ranges):
        row += 1
        tk.Label(dialog, text=label).grid(row=row, column=0, padx=12, pady=3, sticky="w")
        value_var = tk.StringVar(value=f"{value:.6g}")
        entry = tk.Entry(dialog, textvariable=value_var, width=16)
        entry.grid(row=row, column=1, padx=12, pady=3, sticky="w")
        fixed_var = tk.BooleanVar(value=False)
        tk.Checkbutton(dialog, variable=fixed_var).grid(row=row, column=2, padx=12, pady=3, sticky="w")
        tk.Label(dialog, text=allowed).grid(row=row, column=3, padx=12, pady=3, sticky="w")
        entries.append((label, value_var, fixed_var))

    note = (
        "Note: A and background decay are linear NNLS variables unless fixed. "
        "Background IRF is only applied when fixed; otherwise it remains 0 as in the original script."
    )
    tk.Label(dialog, text=note, wraplength=520, justify="left", fg="gray35").grid(
        row=row + 1, column=0, columnspan=4, padx=12, pady=(8, 0), sticky="w"
    )

    button_frame = tk.Frame(dialog)
    button_frame.grid(row=row + 2, column=0, columnspan=4, padx=12, pady=(10, 12), sticky="e")

    def submit():
        values = []
        fixed_flags = []
        try:
            for label, value_var, fixed_var in entries:
                text = value_var.get().strip()
                if not text:
                    raise ValueError(f"{label} is blank.")
                value = float(text)
                if not np.isfinite(value):
                    raise ValueError(f"{label} must be a finite number.")
                values.append(value)
                fixed_flags.append(bool(fixed_var.get()))

            lower_bounds, upper_bounds = make_bounds(n_exponentials)
            values_arr = np.asarray(values, dtype=float)
            if np.any(values_arr < np.asarray(lower_bounds)) or np.any(values_arr > np.asarray(upper_bounds)):
                raise ValueError(
                    "One or more values are outside the allowed ranges shown in the dialog."
                )

        except ValueError as exc:
            messagebox.showerror("Invalid initial estimate", str(exc), parent=dialog)
            return

        result["ok"] = True
        result["initial_params"] = np.asarray(values, dtype=float)
        result["fixed_params"] = np.asarray(fixed_flags, dtype=bool)
        dialog.destroy()

    def cancel():
        dialog.destroy()

    tk.Button(button_frame, text="Cancel", command=cancel).pack(side="right", padx=(6, 0))
    tk.Button(button_frame, text="Run fit", command=submit).pack(side="right")

    dialog.bind("<Return>", lambda event: submit())
    dialog.bind("<Escape>", lambda event: cancel())

    parent.wait_window(dialog)
    return result if result["ok"] else None

def ask_fit_options(parent):
    """
    Ask the user for the number of exponentials and whether to save a DAT file.
    Called after the decay and IRF files are selected.
    """
    dialog = tk.Toplevel(parent)
    dialog.title("Fit Options")
    dialog.resizable(False, False)
    dialog.grab_set()

    n_exp_var = tk.StringVar(value="1")
    save_dat_var = tk.BooleanVar(value=False)
    result = {"ok": False, "n_exponentials": None, "save_dat": False}

    tk.Label(dialog, text="Number of exponentials:").grid(
        row=0, column=0, padx=12, pady=(12, 6), sticky="w"
    )
    n_exp_entry = tk.Entry(dialog, textvariable=n_exp_var, width=10)
    n_exp_entry.grid(row=0, column=1, padx=12, pady=(12, 6), sticky="w")

    tk.Checkbutton(
        dialog,
        text="Generate .dat file",
        variable=save_dat_var,
    ).grid(row=1, column=0, columnspan=2, padx=12, pady=6, sticky="w")

    button_frame = tk.Frame(dialog)
    button_frame.grid(row=2, column=0, columnspan=2, padx=12, pady=(6, 12), sticky="e")

    def submit():
        try:
            n_exponentials = int(n_exp_var.get().strip())
            if n_exponentials < 1:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Invalid input",
                "Please enter a whole number greater than or equal to 1.",
                parent=dialog,
            )
            return

        result["ok"] = True
        result["n_exponentials"] = n_exponentials
        result["save_dat"] = bool(save_dat_var.get())
        dialog.destroy()

    def cancel():
        dialog.destroy()

    tk.Button(button_frame, text="Cancel", command=cancel).pack(side="right", padx=(6, 0))
    tk.Button(button_frame, text="OK", command=submit).pack(side="right")

    n_exp_entry.focus_set()
    n_exp_entry.selection_range(0, tk.END)
    dialog.bind("<Return>", lambda event: submit())
    dialog.bind("<Escape>", lambda event: cancel())

    parent.wait_window(dialog)
    return result if result["ok"] else None


def run_gui_workflow():
    root = tk.Tk()
    root.withdraw()

    decay_sdt_path = filedialog.askopenfilename(
        title="Select decay SDT file",
        filetypes=[
            ("SDT files", "*.sdt"),
            ("All files", "*.*"),
        ],
    )

    if not decay_sdt_path:
        messagebox.showinfo("Cancelled", "No decay SDT file selected.")
        return

    irf_sdt_path = filedialog.askopenfilename(
        title="Select IRF SDT file",
        filetypes=[
            ("SDT files", "*.sdt"),
            ("All files", "*.*"),
        ],
    )

    if not irf_sdt_path:
        messagebox.showinfo("Cancelled", "No IRF SDT file selected.")
        return

    options = ask_fit_options(root)
    if options is None:
        messagebox.showinfo("Cancelled", "Fit options were not selected.")
        return

    try:
        preview_time_ns, preview_decay = load_sdt_decay(decay_sdt_path, block=0)
        default_initial_params = make_initial_params(preview_decay, options["n_exponentials"])
    except Exception as e:
        messagebox.showerror("Error", f"Could not load decay data for initial estimates:\n{e}")
        return

    initial_options = ask_initial_parameter_estimates(
        root,
        options["n_exponentials"],
        default_initial_params,
    )
    if initial_options is None:
        messagebox.showinfo("Cancelled", "Initial fit parameter estimates were not selected.")
        return

    output_dat_path = None
    if options["save_dat"]:
        output_dat_path = filedialog.asksaveasfilename(
            title="Save output DAT file",
            defaultextension=".dat",
            filetypes=[
                ("DAT files", "*.dat"),
                ("All files", "*.*"),
            ],
            initialfile=os.path.splitext(os.path.basename(decay_sdt_path))[0]
            + f"_{options['n_exponentials']}exp_fit.dat",
        )

        if not output_dat_path:
            messagebox.showinfo("Cancelled", "No output DAT file selected.")
            return

    try:
        fit_result = fit_multi_exp_reconv(
            decay_sdt_path=decay_sdt_path,
            irf_sdt_path=irf_sdt_path,
            n_exponentials=options["n_exponentials"],
            decay_block=0,
            irf_block=0,
            fit_start_ns=None,
            fit_end_ns=None,
            initial_params=initial_options["initial_params"],
            fixed_params=initial_options["fixed_params"],
        )
        plot_fit_and_residuals(fit_result)

        if output_dat_path is not None:
            write_fluofit_style_dat(
                fit_result=fit_result,
                output_dat_path=output_dat_path,
                decay_sdt_path=decay_sdt_path,
                irf_sdt_path=irf_sdt_path,
            )

        msg_lines = [
            "Fit complete.",
            "",
            f"Number of exponentials: {fit_result['n_exponentials']}",
        ]

        if output_dat_path is not None:
            msg_lines.extend(["", f"Output DAT file:\n{output_dat_path}"])
        else:
            msg_lines.extend(["", "DAT file was not generated."])

        msg_lines.append("")
        msg_lines.append("Initial estimates used:")
        initial_amplitudes, initial_taus, initial_shift, initial_bg_decay, initial_bg_irf = unpack_multi_exp_params(
            fit_result["initial_params"], fit_result["n_exponentials"]
        )
        fixed_params = np.asarray(fit_result["fixed_params"], dtype=bool)
        fixed_amplitudes = fixed_params[0 : 2 * fit_result["n_exponentials"] : 2]
        fixed_taus = fixed_params[1 : 2 * fit_result["n_exponentials"] : 2]
        for idx, (initial_amplitude, initial_tau, fixed_a, fixed_tau) in enumerate(
            zip(initial_amplitudes, initial_taus, fixed_amplitudes, fixed_taus), start=1
        ):
            msg_lines.append(f"Initial A{idx} [Cnts]: {initial_amplitude:.5f}" + (" (fixed)" if fixed_a else ""))
            msg_lines.append(f"Initial tau{idx} [ns]: {initial_tau:.5f}" + (" (fixed)" if fixed_tau else ""))
        msg_lines.append(f"Initial Shift IRF [ns]: {initial_shift:.5f}" + (" (fixed)" if fixed_params[2 * fit_result["n_exponentials"]] else ""))
        msg_lines.append(f"Initial background decay [Cnts]: {initial_bg_decay:.5f}" + (" (fixed)" if fixed_params[2 * fit_result["n_exponentials"] + 1] else ""))
        msg_lines.append(f"Initial background IRF [Cnts]: {initial_bg_irf:.5f}" + (" (fixed)" if fixed_params[2 * fit_result["n_exponentials"] + 2] else ""))

        msg_lines.append("")
        msg_lines.append("Fitted parameters (LifeFit-style NNLS amplitudes/background):")
        for idx, (amplitude, tau_ns, amp_frac, int_frac) in enumerate(
            zip(
                fit_result["amplitudes_counts"],
                fit_result["taus_ns"],
                fit_result["amplitude_fractions"],
                fit_result["intensity_fractions"],
            ),
            start=1,
        ):
            msg_lines.append(f"A{idx} [Cnts]: {amplitude:.5f}")
            msg_lines.append(f"tau{idx} [ns]: {tau_ns:.5f}")
            msg_lines.append(f"Amplitude fraction {idx}: {100.0 * amp_frac:.2f}%")
            msg_lines.append(f"Intensity fraction {idx}: {100.0 * int_frac:.2f}%")

        msg_lines.extend(
            [
                "",
                f"Shift IRF [ns]: {fit_result['shift_irf_ns']:.5f}",
                #f"IRF offset subtracted: {fit_result['irf_offset_subtracted']:.5f}",
                f"Reduced chi-square: {fit_result['reduced_chi2']:.5f}",
                f"Fitted points: {fit_result['fitted_points']}",
                f"Decay regularization: {fit_result['decay_regularization_info']}",
            ]
        )

        messagebox.showinfo("Done", "\n".join(msg_lines))

    except Exception as e:
        messagebox.showerror("Error", str(e))


if __name__ == "__main__":
    run_gui_workflow()
