import os
from datetime import datetime

import numpy as np
from scipy.optimize import least_squares
from sdtfile import SdtFile

import tkinter as tk
from tkinter import filedialog, messagebox


# Preview of data then select exponentials, see if there's a one size fits all exponential based on a numerical input
# Giving one exponential or fit until chi square is below 2
# Give guesses for exponential values for fitting
# Background: ignore the baselines 
# Shows the life time and amplitude ratios in the plot, and chi squared as a part of the plot 
# fancy error functions, A1 A2, amplitutde and life time errors 
# Amplitude ratios, what is the intensity ratio between A1, A2, and A3 



# ============================================================
# SDT loading
# ============================================================

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

# ============================================================
# Model functions
# ============================================================

def shift_curve(y, shift_bins):
    """
    Fractionally shift a curve using interpolation.

    Positive shift_bins shifts the curve to the right.
    Negative shift_bins shifts the curve to the left.
    """
    x = np.arange(len(y))
    return np.interp(x - shift_bins, x, y, left=0.0, right=0.0)


def reconvolved_single_exp(params, time_ns, irf, dt):
    """
    Single-exponential reconvolution model.

    params:
        A1, tau1_ns, shift_irf_ns, background_decay, background_irf
    """
    A1, tau1_ns, shift_irf_ns, background_decay, background_irf = params

    shift_bins = shift_irf_ns / dt
    shifted_irf = shift_curve(irf + background_irf, shift_bins)

    t0 = time_ns[0]
    pure_decay = A1 * np.exp(-(time_ns - t0) / tau1_ns)
    pure_decay[time_ns < t0] = 0.0

    model = np.convolve(shifted_irf, pure_decay, mode="full")[: len(time_ns)]

    irf_area = shifted_irf.sum()
    if irf_area > 0:
        model = model / irf_area

    return model + background_decay


def poisson_residuals(params, time_ns, decay, irf, dt, fit_mask):
    model = reconvolved_single_exp(params, time_ns, irf, dt)

    # Poisson weighting, common for TCSPC counts.
    sigma = np.sqrt(np.maximum(decay, 1.0))

    return (decay[fit_mask] - model[fit_mask]) / sigma[fit_mask]



# Fitting workflow


def fit_single_exp_reconv(
    decay_sdt_path,
    irf_sdt_path,
    decay_block=0,
    irf_block=0,
    fit_start_ns=5,
    fit_end_ns=48,
):
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

    dt = float(np.median(np.diff(time_ns)))

    if fit_start_ns is None:
        fit_start_ns = float(time_ns[np.argmax(decay)])

    if fit_end_ns is None:
        fit_end_ns = float(time_ns[-1])

    fit_mask = (time_ns >= fit_start_ns) & (time_ns <= fit_end_ns)

    # Initial guesses.
    p0 = np.array([
        np.max(decay),  # A1
        5.5,            # tau1, ns
        -0.04,          # IRF shift, ns
        0.0,            # decay background
        0.0,            # IRF background
    ], dtype=float)

    lower_bounds = [0.0, 0.01, -5.0, 0.0, 0.0]
    upper_bounds = [np.inf, 100.0, 5.0, np.inf, np.inf]

    result = least_squares(
        poisson_residuals,
        p0,
        bounds=(lower_bounds, upper_bounds),
        args=(time_ns, decay, irf, dt, fit_mask),
        max_nfev=20000,
    )

    params = result.x

    model = reconvolved_single_exp(params, time_ns, irf, dt)

    full_mask = np.ones_like(time_ns, dtype=bool)
    residuals = poisson_residuals(params, time_ns, decay, irf, dt, full_mask)

    fit_residuals = poisson_residuals(params, time_ns, decay, irf, dt, fit_mask)
    n_fit_points = int(fit_mask.sum())
    n_params = len(params)

    reduced_chi2 = np.sum(fit_residuals ** 2) / max(n_fit_points - n_params, 1)

    return {
        "time_ns": time_ns,
        "decay": decay,
        "irf": irf,
        "fit": model,
        "residuals": residuals,
        "fit_mask": fit_mask,
        "A1_counts": params[0],
        "tau1_ns": params[1],
        "shift_irf_ns": params[2],
        "background_decay_counts": params[3],
        "background_irf_counts": params[4],
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
    title="Python single-exponential reconvolution fit",
):
    time_ns = np.asarray(fit_result["time_ns"])
    decay = np.asarray(fit_result["decay"])
    irf = np.asarray(fit_result["irf"])
    model = np.asarray(fit_result["fit"])
    residuals = np.asarray(fit_result["residuals"])

    A1 = fit_result["A1_counts"]
    tau1 = fit_result["tau1_ns"]
    shift_irf = fit_result["shift_irf_ns"]
    bg_decay = fit_result["background_decay_counts"]
    bg_irf = fit_result["background_irf_counts"]
    reduced_chi2 = fit_result["reduced_chi2"]
    fitted_points = fit_result["fitted_points"]

    decay_name = os.path.basename(decay_sdt_path)
    irf_name = os.path.basename(irf_sdt_path)

    now = datetime.now().strftime("%m/%d/%Y %I:%M:%S %p")

    lines = []

    lines.append("PicoQuant FluoFit-style Python Export")
    lines.append(f"Saved : {now}")
    lines.append("")
    lines.append(title)
    lines.append("")
    lines.append("Model: Exp. [Reconv.] (Exponential)")
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

    lines.append(
        f"A1 [Cnts]        {A1:14.5f}              ---            ---     Fitting"
    )
    lines.append(
        f"t1 [ns]          {tau1:14.5f}              ---            ---     Fitting"
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
    lines.append(f"tAv.1={tau1:.5f} ns (intensity weighted)")
    lines.append(f"tAv.2={tau1:.5f} ns (amplitude weighted)")
    lines.append("")

    lines.append("Fractional Intensities of the Positive Decay Components:")
    lines.append(f"t1 ({tau1:.5f} ns) : 100.00%")
    lines.append("")

    lines.append("Fractional Amplitudes of the Positive Decay Components:")
    lines.append(f"t1 ({tau1:.5f} ns) : 100.00%")
    lines.append("")

    lines.append("Fitted Decay and Exponential Components:")
    lines.append("     Time [ns]      t1 [kCounts]      Sum [kCounts]")

    component_time = np.linspace(0.0, 22.0, 1000)
    component = A1 * np.exp(-component_time / tau1)

    for t, c in zip(component_time, component):
        lines.append(f"{t:12.5f}{c / 1000.0:16.8f}{c / 1000.0:16.8f}")

    lines.append("")
    lines.append("Confidence Intervals:")
    lines.append("A1 : [--- ; ---] Cnts")
    lines.append("t1 : [--- ; ---] ns")
    lines.append("Shift IRF : [--- ; ---] ns")
    lines.append("")

    text = "\r\n".join(lines)

    with open(output_dat_path, "w", encoding="latin-1", newline="") as f:
        f.write(text)

    return output_dat_path

import matplotlib.pyplot as plt
import numpy as np


def plot_fit_and_residuals(fit_result):
    time_ns = fit_result["time_ns"]
    decay = fit_result["decay"]
    irf = fit_result["irf"]
    model = fit_result["fit"]
    residuals = fit_result["residuals"]

    # Create a figure with 2 subplots stacked vertically
    # gridspec_kw={'height_ratios': [3, 1]} makes the fit plot larger than the residual plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )

    # 1. Top Subplot: Regression / fit plot (Log scale)
    ax1.semilogy(
        time_ns, np.maximum(decay, 1), ".", markersize=2, label="Decay"
    )
    ax1.semilogy(
        time_ns, np.maximum(model, 1), "-", linewidth=1.2, label="Model fit"
    )
    ax1.semilogy(time_ns, np.maximum(irf, 1), "-", linewidth=1.0, label="IRF")
    ax1.set_ylabel("Intensity [counts]")
    ax1.set_title("Decay Regression Fit and Residuals")
    ax1.legend()

    # 2. Bottom Subplot: Residual plot (Linear scale)
    ax2.plot(time_ns, residuals, linewidth=0.8, color="purple")
    ax2.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax2.set_xlabel("Time [ns]")
    ax2.set_ylabel("Residual")

    # Clean up layout and display
    plt.tight_layout()
    plt.show()
# ============================================================
# Tkinter file dialogs
# ============================================================

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

    output_dat_path = filedialog.asksaveasfilename(
        title="Save output DAT file",
        defaultextension=".dat",
        filetypes=[
            ("DAT files", "*.dat"),
            ("All files", "*.*"),
        ],
        initialfile=os.path.splitext(os.path.basename(decay_sdt_path))[0] + "_fit.dat",
    )

    if not output_dat_path:
        messagebox.showinfo("Cancelled", "No output DAT file selected.")
        return

    try:
        fit_result = fit_single_exp_reconv(
            decay_sdt_path=decay_sdt_path,
            irf_sdt_path=irf_sdt_path,
            decay_block=0,
            irf_block=0,
            fit_start_ns=None,
            fit_end_ns=None,
        )
        plot_fit_and_residuals(fit_result)
        write_fluofit_style_dat(
            fit_result=fit_result,
            output_dat_path=output_dat_path,
            decay_sdt_path=decay_sdt_path,
            irf_sdt_path=irf_sdt_path,
        )

        msg = (
            "Fit complete.\n\n"
            f"Output file:\n{output_dat_path}\n\n"
            f"A1 [Cnts]: {fit_result['A1_counts']:.5f}\n"
            f"tau1 [ns]: {fit_result['tau1_ns']:.5f}\n"
            f"Shift IRF [ns]: {fit_result['shift_irf_ns']:.5f}\n"
            f"Reduced chi-square: {fit_result['reduced_chi2']:.5f}\n"
            f"Fitted points: {fit_result['fitted_points']}"
        )

        messagebox.showinfo("Done", msg)

    except Exception as e:
        messagebox.showerror("Error", str(e))
import matplotlib.pyplot as plt
import numpy as np



if __name__ == "__main__":
    run_gui_workflow()