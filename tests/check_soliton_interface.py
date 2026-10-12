"""Exercise hidden Tk dialogs and their real callbacks without manual input."""
from pathlib import Path
import runpy
import sys
import tempfile
from unittest.mock import patch

import numpy as np
from matplotlib.backend_bases import MouseEvent
import tkinter as tk
from tkinter import ttk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def descendants(widget):
    for child in widget.winfo_children():
        yield child
        yield from descendants(child)


for relative in ('Lab Interface.py', 'LabInterfaceApp/lab_interface.py'):
    ns = runpy.run_path(str(ROOT / relative), run_name='interface_check')
    cls = ns['PsTAAnalysisApp']
    cls._maximize_main_window = lambda self: self.withdraw()
    app = cls()
    original_toplevel = tk.Toplevel

    def hidden_dialog(*args, **kwargs):
        dlg = original_toplevel(*args, **kwargs)
        dlg.withdraw()
        return dlg

    try:
        with tempfile.TemporaryDirectory() as folder:
            times = np.arange(-3., 4.)
            matrix = np.column_stack((times, np.arange(28.).reshape(7, 4)))
            np.savetxt(Path(folder) / 'scan1.txt', matrix)
            dat = ns['SolitonTAExperimentSelfContained'](folder, 'scan1.txt', 1)
            dat.mainSoliton_SC([0., 2., 3., 4.], False, -2.)
            app.raw_dat = app.proc_dat = dat
            app.pipeline_type = 'soliton'
            app.soliton_has_energy_calibration = True
            app._reset_spectral_mode_on_pull = True
            app._pull_current_dat_arrays()
            app.plot_time_slices()
            app.plot_wavelength_slices()
            with patch.object(tk, 'Toplevel', hidden_dialog), patch.object(ns['messagebox'], 'showinfo'), patch.object(ns['messagebox'], 'showerror') as errors:
                app.dispersion_correction()
                chooser = next(w for w in app.winfo_children() if isinstance(w, original_toplevel))
                assert any(isinstance(w, ttk.Button) and w.cget('text') == 'Interactive Point-Selection Correction' for w in descendants(chooser))
                chooser.destroy()
                app._open_soliton_chirp_dialog()
                dlg = next(w for w in app.winfo_children() if isinstance(w, original_toplevel))
                assert dlg.minsize() == (650, 240)
                entry = next(w for w in descendants(dlg) if isinstance(w, ttk.Entry))
                apply = next(w for w in descendants(dlg) if isinstance(w, ttk.Button) and w.cget('text') == 'Apply')
                entry.delete(0, 'end')
                entry.insert(0, 'bad input')
                apply.invoke()
                assert dlg.winfo_exists()
                entry.delete(0, 'end')
                entry.insert(0, '0, 0, 0, 0, 0.5')
                apply.invoke()
                assert not dlg.winfo_exists()
                np.testing.assert_allclose(dat.chirpCorrectionCoefficients, [0., 0., 0., 0., .5])
                np.testing.assert_allclose(app.data, dat.TAMeanSortedBackgroundSub_CC[:, 1:])

                canvases = []
                globals_dict = cls._open_interactive_dispersion_dialog.__globals__
                original_canvas = globals_dict['FigureCanvasTkAgg']

                def capture_canvas(*args, **kwargs):
                    canvas = original_canvas(*args, **kwargs)
                    canvases.append(canvas)
                    return canvas

                with patch.dict(globals_dict, FigureCanvasTkAgg=capture_canvas):
                    app._open_interactive_dispersion_dialog()
                dlg = next(w for w in app.winfo_children() if isinstance(w, original_toplevel))
                canvas = canvases[0]
                canvas.draw()
                axes = canvas.figure.axes[0]
                # Channel 1 is masked: preserve original pixels 2, 3, 4.
                for pixel in (2., 3., 4.):
                    x, y = axes.transData.transform((pixel, .1 * pixel))
                    event = MouseEvent('button_press_event', canvas, x, y, button=1)
                    canvas.callbacks.process('button_press_event', event)
                apply = next(w for w in descendants(dlg) if isinstance(w, ttk.Button) and w.cget('text') == 'Fit and Apply')
                apply.invoke()
                assert not dlg.winfo_exists()
                np.testing.assert_allclose(np.polyval(dat.chirpCorrectionCoefficients, [2., 3., 4.]), [.2, .3, .4], atol=.02)
                assert dat.chirpCorrectionCoefficients.size == 5
                np.testing.assert_allclose(app.data, dat.TAMeanSortedBackgroundSub_CC[:, 1:])
                with patch.object(ns['simpledialog'], 'askfloat', return_value=-1.):
                    app.background_correction()
                np.testing.assert_allclose(app.data, dat.TAMeanSortedBackgroundSub_CC[:, 1:])
                spectral_line = app.ax_time.lines[1]
                target_time = float(app.time_slice_var.get().split(',')[0])
                row = np.argmin(np.abs(app.times - target_time))
                np.testing.assert_allclose(spectral_line.get_ydata(), app.data[row])
                kinetic_line = app.ax_wl.lines[1]
                target_energy = float(app.wl_slice_var.get().split(',')[0])
                column = np.argmin(np.abs(app.wavelength - target_energy))
                np.testing.assert_allclose(kinetic_line.get_ydata(), app.data[:, column])
                assert '-1 ps' in app.status_label.cget('text')
                errors.assert_not_called()
        print(relative + ': dialogs, coefficient validation, interactive fit, and background refresh passed')
    finally:
        app.destroy()
