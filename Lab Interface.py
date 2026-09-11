import tkinter as tk
from tkinter import ttk
from tkinter import messagebox  
from tkinter import filedialog
from tkinter import simpledialog
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import matlab.engine
import threading
from pathlib import Path
import os
from datetime import datetime

try:
    from scipy.optimize import least_squares
except ImportError:
    least_squares = None

try:
    from sdtfile import SdtFile
except ImportError:
    SdtFile = None

# Wavelength and time slices based on comma separated list of values
# implement min and max fitting for the heatmap where any value above the maximum or below the minimum becomes it.
# Ideally use the bar on the heatmap  

# 03/26 1. Call functions from new matlab class
# 2. Wavelength to energy space - Note: When you go from energy to wavelength for absorption it does not matter, but for counts you need apply a jacobian. 
##Low prio Getting the basics right -> Anyything that can be seen as a historgam florescence and phosphoflorescence needs adjustments 
# 1 - 10% excitation fraction 
# 3. Plot the ground state spectrum, the florescnece as well 

# Time-Domain is actually spectral slices -> Pick range + plot ground state and florescence 
# Spectral slices average between some range of delays 1 dropdown menu on Time(ps), comma to implement table for averaging ranges of times
### Table: top is top bounds, bottom is bottom bounds 
# Wavelength-Domain is Kinetics trace 

# Normalization -> kinetics traces point divided by max possible val 
# waterfall: shift apart the different spectral slices 
# Compare across objects on the same graph (time traces: solvent vs experiment)
# Ideally in slices and traces have both the option to analyze raw and preprocessed 



# When plotting Ground state and Fluorescence have the user input scaling factor
# Plot: scaling factor * absorbance at any given wavelength. 
# After the buttons have been clicked have a box for show or don't show for fluorescece and ground state

# For wavelength average within ranges to smooth it out, binning 
# binning should be optional for both spectral and time 
# select ranges of time
# Ranges for the 2D map view 
# Remove pubchem search 
# And declasssify
# Comma separated lists or a table that averages wavelength ranges instead of discrete individual wavelengths. 
# [Range , Range] average to obtain value at time.


# Packaging the software to be readily sendable. 
# .exe is ideal 

# TCSPC
# .Sdt files
# https://github.com/cgohlke/sdtfile/
# https://github.com/glotaran/pyglotaran
# Instrument response function 
# Convolution of IRF(instrument response function ) with exp^-t/lambda 
# scaling factor exists 
# Select 1 2 or 3 exponentials
# Linear vs log scaling linear on x log on y 


# 1. Cut pubchem search out
# 2. Implement new interface features
# 3. TCSPC functionality reconstruction via python
# 4. Later Matlab translated into python 

# Use RdBu as color map for matplotlib 
# Plotly 
class PsTAAnalysisApp(tk.Tk):
    """
    A Python tkinter application replicating the psTA Analysis Suite layout
    described for MATLAB App Designer.
    """ 
    def __init__(self):
        super().__init__()

        # Main App Window Setup
        self.title("psTA Analysis Suite (Python)")

        # This interface is designed primarily for a full-screen/maximized lab
        # workstation.  Keep a sensible minimum size for smaller monitors, then
        # maximize using the platform-specific mechanism where available.
        self.minsize(1280, 720)
        self._maximize_main_window()
        
        # Main Grid Layout (3 Columns)
        # Column 0: analysis/data controls
        self.grid_columnconfigure(0, weight=2, minsize=280)
        # Column 1: plots/visualization; this remains the largest area.
        self.grid_columnconfigure(1, weight=5, minsize=680)
        # Column 2: fitting/range controls.  TCSPC parameter entry needs room
        # for parameter names, editable estimates, fixed flags, and bounds.
        self.grid_columnconfigure(2, weight=3, minsize=440)
        # Row 0: Main content
        self.grid_rowconfigure(0, weight=1)
        # Row 1: Status Bar
        self.grid_rowconfigure(1, weight=0)

        # Create the Three Major Panels
        self.left_panel = ttk.Frame(self, padding="10")
        self.middle_panel = ttk.Frame(self, padding="10")
        self.right_panel = ttk.Frame(self, padding="10")

        self.left_panel.pack_propagate(False)
        self.middle_panel.pack_propagate(False)
        self.right_panel.pack_propagate(False)

        # Place panels on the grid
        self.left_panel.grid(row=0, column=0, sticky="nsew")
        self.middle_panel.grid(row=0, column=1, sticky="nsew")
        self.right_panel.grid(row=0, column=2, sticky="nsew")

        # Populate Each Panel 
        self._create_left_panel()
        self._create_middle_panel()
        self._create_right_panel()

        # Bottom Status Bar 
        self._create_status_bar()

        # Data Containers 
        # Data is the container for the primary data - ie experiment being analzyed 
        # data 2 is a place holder for TCSPC for future use
        # 1 set of wavelength and times place holder for now, if need be duplica
        self.data = None
        self.data2 = None
        self.wavelength = None
        self.times = None

        self.raw_data = None
        self.proc_data = None
        self.raw_dat = None
        self.proc_dat = None


        self.pipeline_type = None 

        # Optional spectral overlays shown on the Time-Domain Traces tab.
        # Each overlay stores the source CSV, scale factor, plotted line, and visibility state.
        self.time_trace_overlays = {
            "ground_state": {"label": "Ground State", "line": None, "data": None, "scale": None},
            "fluorescence": {"label": "Fluorescence", "line": None, "data": None, "scale": None},
        }

        # Matlab Environment Start
        self.eng = matlab.engine.start_matlab()
        # Add path to MATLAB script
        self.eng.addpath(r"MatlabItems", nargout=0)
        self.eng.addpath(r"MatlabItems", nargout=0)
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.map_colorbar = None

    def _maximize_main_window(self):
        """Open the lab interface maximized while retaining normal window controls."""
        # Windows/Tk commonly supports state("zoomed").  Many Linux Tk builds
        # use the -zoomed attribute instead.  The screen-size fallback still
        # gives a full-screen-oriented layout on other platforms.
        try:
            self.state("zoomed")
            return
        except tk.TclError:
            pass

        try:
            self.attributes("-zoomed", True)
            return
        except tk.TclError:
            pass

        try:
            width = self.winfo_screenwidth()
            height = self.winfo_screenheight()
            self.geometry(f"{width}x{height}+0+0")
        except tk.TclError:
            self.geometry("1600x900")

    

    def _create_left_panel(self):
        """Populates the Left Panel (Data Control & Processing)"""
        frame = self.left_panel
        
        # Section: Data I/O 
        data_io_frame = ttk.LabelFrame(frame, text="Data I/O", padding="10")
        data_io_frame.pack(fill='x', pady=5)

        ttk.Button(
            data_io_frame, text="Load psTA Data", command=self.load_data
        ).pack(fill='x', pady=2)

        ttk.Button(
            data_io_frame, text="Load TCSPC Decay", command=self.import_tcspc_decay).pack(fill='x', pady=2)
        
        
        ttk.Button(
            data_io_frame, text="Load nsTA Data (optional)", command=self.placeholder_command
        ).pack(fill='x', pady=2)
        
        self.combine_zeros_var = tk.BooleanVar()
        ttk.Checkbutton(
            data_io_frame, 
            text="Include zeros (1 ns–50 ns)", 
            variable=self.combine_zeros_var
        ).pack(fill='x', pady=5)

        #  Preprocessing 
        preproc_frame = ttk.LabelFrame(frame, text="Preprocessing", padding="10")
        preproc_frame.pack(fill='x', pady=5)

        # Background Correction controls (user input + run)
        ttk.Button(
            preproc_frame, text="Background Correction", command=self.background_correction
        ).pack(fill='x', pady=2)
        
        ttk.Button(
            preproc_frame, text="Dispersion Correction", command=self.dispersion_correction
        ).pack(fill='x', pady=2)

        # Sub-section: Outlier Removal
        outlier_frame = ttk.LabelFrame(preproc_frame, text="Outlier Removal", padding="10")
        outlier_frame.pack(fill='x', pady=10)
        
        ttk.Label(outlier_frame, text="Sensitivity (0-1):").pack()
        self.sensitivity_slider = ttk.Scale(outlier_frame, from_=0, to=1, orient='horizontal')
        self.sensitivity_slider.pack(fill='x', pady=2)

        self.show_before_after_var = tk.BooleanVar()
        ttk.Checkbutton(
            outlier_frame, 
            text="Show before/after plots", 
            variable=self.show_before_after_var
        ).pack(fill='x', pady=2)

        ttk.Label(outlier_frame, text="Select wavelength (nm):").pack()
        self.wavelength_combo = ttk.Combobox(
            outlier_frame, 
            values=["(no data loaded)", "500", "520", "550"]
        )
        self.wavelength_combo.current(0)
        self.wavelength_combo.pack(fill='x', pady=2)

        ttk.Button(
            outlier_frame, text="Clean Outliers", command=self.placeholder_command
        ).pack(fill='x', pady=5)
        # End Outlier Sub-section

        ttk.Button(
            preproc_frame, text="Recompute Average", command=self.placeholder_command
        ).pack(fill='x', pady=5)
        

    def _create_middle_panel(self):
        """Populates the Middle Panel (Visualization & Diagnostics)"""
        frame = self.middle_panel
        
        # Create the Tab Group (Notebook)
        self.tab_control = ttk.Notebook(frame)
        tab_control = self.tab_control
        
        # Tab Frames
        tab1 = ttk.Frame(tab_control, padding="10")
        tab2 = ttk.Frame(tab_control, padding="10")
        tab3 = ttk.Frame(tab_control, padding="10")
        tab4 = ttk.Frame(tab_control, padding="10")
        tab5 = ttk.Frame(tab_control, padding="10")
        self.tcspc_tab = tab5
        
        # Add tabs to the notebook
        tab_control.add(tab1, text='Time-Domain Traces')
        tab_control.add(tab4, text='Wavelength-Domain Traces')
        tab_control.add(tab2, text='2D Map View')
        tab_control.add(tab3, text='Fitting')
        tab_control.add(tab5, text='TCSPC Exponential Fit')
        
        # Make the notebook fill the middle panel
        tab_control.pack(expand=1, fill='both')
        tab_control.bind('<<NotebookTabChanged>>', self._on_middle_tab_changed)

        # Populate Tab 1: Time-Domain Traces 
        tab1_controls = ttk.Frame(tab1)
        tab1_controls.pack(fill='x', pady=5)

        ttk.Label(tab1_controls, text="Times (ps), comma-separated:").grid(row=0, column=0, sticky='w', padx=5, pady=2)
        self.time_slice_var = tk.StringVar(value="0.1, 1, 10, 100")
        ttk.Entry(tab1_controls, textvariable=self.time_slice_var, width=35).grid(row=0, column=1, sticky='w', padx=5, pady=2)
        ttk.Button(tab1_controls, text="Plot", command=self.plot_time_slices).grid(row=0, column=2, padx=2, pady=2)
        ttk.Button(
            tab1_controls, text="Clear",
            command=lambda: self._clear_slice_axis(
                self.ax_time, self.canvas_time,
                "Wavelength (nm)", "ΔA", "Time-Domain Traces"
            )
        ).grid(row=0, column=3, padx=2, pady=2)

        # New wavelength range controls for the time-domain traces plot only
        ttk.Label(tab1_controls, text="Visible wavelength range (nm):").grid(row=1, column=0, sticky='w', padx=5, pady=2)
        self.time_trace_wl_range_var = tk.StringVar(value="")
        ttk.Entry(tab1_controls, textvariable=self.time_trace_wl_range_var, width=35).grid(row=1, column=1, sticky='w', padx=5, pady=2)
        ttk.Button(tab1_controls, text="Apply Range", command=self.apply_time_trace_wavelength_range).grid(row=1, column=2, padx=2, pady=2)
        ttk.Button(tab1_controls, text="Reset Range", command=self.reset_time_trace_wavelength_range).grid(row=1, column=3, padx=2, pady=2)

        # Optional CSV overlays for ground-state and fluorescence spectra.
        # CSV headers expected: Wavelength (nm), Absorbance (AU), optional Std.Dev.
        ttk.Button(
            tab1_controls, text="Plot Ground State",
            command=lambda: self.plot_scaled_absorbance_overlay("ground_state")
        ).grid(row=2, column=0, padx=5, pady=6, sticky='w')
        ttk.Button(
            tab1_controls, text="Plot Fluorescence",
            command=lambda: self.plot_scaled_absorbance_overlay("fluorescence")
        ).grid(row=2, column=1, padx=5, pady=6, sticky='w')

        self.ground_state_visible_var = tk.BooleanVar(value=True)
        self.fluorescence_visible_var = tk.BooleanVar(value=True)
        self.ground_state_toggle = ttk.Checkbutton(
            tab1_controls, text="Show Ground State",
            variable=self.ground_state_visible_var,
            command=lambda: self.toggle_time_trace_overlay("ground_state"),
            state="disabled"
        )
        self.ground_state_toggle.grid(row=3, column=0, padx=5, pady=2, sticky='w')
        self.fluorescence_toggle = ttk.Checkbutton(
            tab1_controls, text="Show Fluorescence",
            variable=self.fluorescence_visible_var,
            command=lambda: self.toggle_time_trace_overlay("fluorescence"),
            state="disabled"
        )
        self.fluorescence_toggle.grid(row=3, column=1, padx=5, pady=2, sticky='w')

        tab1_controls.grid_columnconfigure(1, weight=1)

        self.fig_time = Figure(figsize=(5, 4), dpi=100)
        self.ax_time = self.fig_time.add_subplot(111)
        self.ax_time.set_title("Time-Domain Traces")
        self.ax_time.set_xlabel("Wavelength (nm)")
        self.ax_time.set_ylabel("ΔA")
        self.ax_time.grid(True)
        self.ax_time.text(
            0.5, 0.5, "Load data, then enter times above",
            ha='center', va='center', transform=self.ax_time.transAxes,
            fontsize=12, color='gray', alpha=0.5
        )
        self.canvas_time = FigureCanvasTkAgg(self.fig_time, master=tab1)
        self.canvas_time.draw()
        self.canvas_time.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


        # Populate Tab 2: Wavelength-domain Traces
        # Add controls
        tab4_controls = ttk.Frame(tab4)
        tab4_controls.pack(fill='x', pady=5)

        # --- Row 1: Wavelength Input ---
        row1 = ttk.Frame(tab4_controls)
        row1.pack(side='top', fill='x', pady=2)

        ttk.Label(row1, text="Wavelengths (nm), comma-separated:").pack(side='left', padx=5)
        self.wl_slice_var = tk.StringVar(value="500, 550, 600, 650")
        ttk.Entry(row1, textvariable=self.wl_slice_var, width=35).pack(side='left', padx=5)

        # --- Row 2: Average Range & Buttons ---
        row2 = ttk.Frame(tab4_controls)
        row2.pack(side='top', fill='x', pady=2)

        ttk.Label(row2, text="Average range (+/- nm):").pack(side='left', padx=5)
        self.wl_avg_range_var = tk.StringVar(value="0")
        ttk.Entry(row2, textvariable=self.wl_avg_range_var, width=8).pack(side='left', padx=5)
        
        ttk.Button(
            row2, text="Plot", command=self.plot_wavelength_slices
        ).pack(side='left', padx=10)
        
        ttk.Button(
            row2, text="Clear",
            command=lambda: self._clear_slice_axis(
                self.ax_wl, self.canvas_wl,
                "Time (ps)", "ΔA", "Wavelength-Domain Traces"
            )
        ).pack(side='left', padx=2)

        self.fig_wl = Figure(figsize=(5, 4), dpi=100)
        self.ax_wl = self.fig_wl.add_subplot(111)
        self.ax_wl.set_title("Wavelength-Domain Traces")
        self.ax_wl.set_xlabel("Time (ps)")
        self.ax_wl.set_ylabel("ΔA")
        self.ax_wl.grid(True)
        self.ax_wl.text(
            0.5, 0.5, "Load data, then enter wavelengths above",
            ha='center', va='center', transform=self.ax_wl.transAxes,
            fontsize=12, color='gray', alpha=0.5
        )
        self.canvas_wl = FigureCanvasTkAgg(self.fig_wl, master=tab4)
        self.canvas_wl.draw()
        self.canvas_wl.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


        # Populate Tab 3: 2D Map View
    
        tab2_controls = ttk.Frame(tab2)
        tab2_controls.pack(fill='x', pady=5)

        # Min/Max clamp controls
        ttk.Label(tab2_controls, text="Min:").pack(side='left', padx=(5, 2))
        self.map_vmin_var = tk.StringVar(value="")
        self.map_vmin_entry = ttk.Entry(tab2_controls, textvariable=self.map_vmin_var, width=8)
        self.map_vmin_entry.pack(side='left', padx=(0, 5))

        ttk.Label(tab2_controls, text="Max:").pack(side='left', padx=(5, 2))
        self.map_vmax_var = tk.StringVar(value="")
        self.map_vmax_entry = ttk.Entry(tab2_controls, textvariable=self.map_vmax_var, width=8)
        self.map_vmax_entry.pack(side='left', padx=(0, 5))

        ttk.Button(
            tab2_controls, text="Apply", command=self.update_2d_map
        ).pack(side='left', padx=5)

        ttk.Button(
            tab2_controls, text="Reset", command=self.reset_map_clamp
        ).pack(side='left', padx=5)

        self.colorbar_toggle_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            tab2_controls, text="Show Colorbar", variable=self.colorbar_toggle_var
        ).pack(side='left', padx=5)


        # Add Matplotlib Plot Canvas (UIAxes_Map)
        self.map_fig = Figure(figsize=(5, 4), dpi=100)
        self.map_ax = self.map_fig.add_subplot(111)
        self.map_ax.set_title("2D Map View (dat.TAmean)")
        self.map_ax.set_xlabel("Wavelength (nm)")
        self.map_ax.set_ylabel("Time (log scale)")
        self.map_ax.text(0.5, 0.5, "UIAxes_Map", horizontalalignment='center', verticalalignment='center', transform=self.map_ax.transAxes, fontsize=16, color='gray', alpha=0.5)

        self.canvas_tab2 = FigureCanvasTkAgg(self.map_fig, master=tab2)
        self.canvas_tab2.draw()
        self.canvas_tab2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Populate Tab 4: Residuals / Cleaned Comparison
        # This tab will have two plots, so we use a main frame
        tab3_main_frame = ttk.Frame(tab3)
        tab3_main_frame.pack(fill=tk.BOTH, expand=True)
        tab3_main_frame.grid_rowconfigure(0, weight=1)
        tab3_main_frame.grid_rowconfigure(1, weight=1)
        tab3_main_frame.grid_columnconfigure(0, weight=1)
        
        # Top Plot: TA data
        before_fig = Figure(dpi=100)
        before_ax = before_fig.add_subplot(111)
        before_ax.set_title("psTA")
        before_ax.set_ylabel("ΔA")
        before_ax.grid(True)
        before_ax.text(0.5, 0.5, "psTA", horizontalalignment='center', verticalalignment='center', transform=before_ax.transAxes, fontsize=14, color='gray', alpha=0.5)
        
        self.canvas_tab3_top = FigureCanvasTkAgg(before_fig, master=tab3_main_frame)
        self.canvas_tab3_top.draw()
        self.canvas_tab3_top.get_tk_widget().grid(row=0, column=0, sticky="nsew", pady=2)

        # Bottom Plot: TCSPC
        after_fig = Figure(dpi=100)
        after_ax = after_fig.add_subplot(111)
        after_ax.set_title("TCSPC")
        after_ax.set_xlabel("Time (ps/ns)")
        after_ax.set_ylabel("Emission Intensity")
        after_ax.grid(True)
        after_ax.text(0.5, 0.5, "TCSPC", horizontalalignment='center', verticalalignment='center', transform=after_ax.transAxes, fontsize=14, color='gray', alpha=0.5)
        
        self.canvas_tab3_bot = FigureCanvasTkAgg(after_fig, master=tab3_main_frame)
        self.canvas_tab3_bot.draw()
        self.canvas_tab3_bot.get_tk_widget().grid(row=1, column=0, sticky="nsew", pady=2)


        # Populate Tab 5: TCSPC Exponential Fit
        # The standalone fitting workflow is embedded here so the user can load,
        # preview, fit, inspect residuals, and optionally save a report without
        # leaving the main lab interface.
        self.tcspc_decay_path = None
        self.tcspc_irf_path = None
        self.tcspc_decay_time_ns = None
        self.tcspc_decay_raw = None
        self.tcspc_irf_time_ns = None
        self.tcspc_irf_raw = None
        self.tcspc_fit_result = None
        self.tcspc_last_summary = "No fit has been run yet."

        self.tcspc_n_exp_var = tk.StringVar(value="1")
        self.tcspc_fit_start_var = tk.StringVar(value="")
        self.tcspc_fit_end_var = tk.StringVar(value="")
        self.tcspc_report_var = tk.BooleanVar(value=False)
        self.tcspc_estimate_vars = {}
        self.tcspc_fixed_vars = {}

        tcspc_controls = ttk.Frame(tab5)
        tcspc_controls.pack(fill='x', pady=(0, 6))
        tcspc_controls.grid_columnconfigure(0, weight=1)
        tcspc_controls.grid_columnconfigure(1, weight=1)

        ttk.Button(
            tcspc_controls,
            text="Import TCSPC Lifetime Decay",
            command=self.import_tcspc_decay,
        ).grid(row=0, column=0, sticky="ew", padx=(0, 4), pady=2)

        ttk.Button(
            tcspc_controls,
            text="Import TCSPC IRF",
            command=self.import_tcspc_irf,
        ).grid(row=0, column=1, sticky="ew", padx=(4, 0), pady=2)

        self.tcspc_decay_file_label = ttk.Label(
            tcspc_controls, text="Decay: not loaded", anchor="w"
        )
        self.tcspc_decay_file_label.grid(row=1, column=0, sticky="ew", padx=(0, 4), pady=2)

        self.tcspc_irf_file_label = ttk.Label(
            tcspc_controls, text="IRF: not loaded", anchor="w"
        )
        self.tcspc_irf_file_label.grid(row=1, column=1, sticky="ew", padx=(4, 0), pady=2)

        self.tcspc_fig = Figure(figsize=(6, 5), dpi=100)
        self.tcspc_ax = self.tcspc_fig.add_subplot(111)
        self._draw_tcspc_empty_plot()

        self.tcspc_canvas = FigureCanvasTkAgg(self.tcspc_fig, master=tab5)
        self.tcspc_canvas.draw()
        self.tcspc_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


    def _create_right_panel(self):
        """Populates the Right Panel with a dynamic tab-specific top area and shared exports."""
        frame = self.right_panel

        # The top portion changes when the user switches tabs in the middle panel.
        self.right_dynamic_frame = ttk.Frame(frame)
        # Let tab-specific controls use the full height and width of the expanded
        # right panel.  This is especially useful for multi-exponential TCSPC fits.
        self.right_dynamic_frame.pack(fill='both', expand=True, pady=5, side='top')

        # Section: Export Options stays available for every middle-panel tab.
        export_frame = ttk.LabelFrame(frame, text="Export Options", padding="10")
        export_frame.pack(fill='x', pady=5, side='bottom')

        ttk.Button(
            export_frame, text="Save Cleaned Data", command=self.placeholder_command
        ).pack(fill='x', pady=2)

        ttk.Label(export_frame, text="Export Plots:").pack()
        self.export_plot_combo = ttk.Combobox(
            export_frame,
            values=["PNG", "PDF", "SVG"]
        )
        self.export_plot_combo.current(0)
        self.export_plot_combo.pack(fill='x', pady=2)

        ttk.Button(
            export_frame, text="Save Report", command=self.placeholder_command
        ).pack(fill='x', pady=2)

        self._update_right_panel_for_active_tab()

    def _on_middle_tab_changed(self, event=None):
        if hasattr(self, "right_dynamic_frame"):
            self._update_right_panel_for_active_tab()

        self.after_idle(self._normalize_active_tab_layout)

    def _normalize_active_tab_layout(self):
        self.update_idletasks()

        active_tab = self._active_middle_tab_text()
        canvas_groups = {
            "Time-Domain Traces": (self.canvas_time,),
            "Wavelength-Domain Traces": (self.canvas_wl,),
            "2D Map View": (self.canvas_tab2,),
            "Fitting": (self.canvas_tab3_top, self.canvas_tab3_bot),
            "TCSPC Exponential Fit": (self.tcspc_canvas,)
        }

        for canvas in canvas_groups.get(active_tab, ()):
            widget = canvas.get_tk_widget()
            width = max(widget.winfo_width(), 1)
            height = max(widget.winfo_height(), 1)
            widget.event_generate("<Configure>", width=width, height=height)
            canvas.draw_idle()

    def _active_middle_tab_text(self):
        """Return the current middle notebook tab label."""
        if not hasattr(self, "tab_control"):
            return "Fitting"
        selected = self.tab_control.select()
        if not selected:
            return "Fitting"
        return self.tab_control.tab(selected, "text")

    def _clear_right_dynamic_frame(self):
        """Remove the tab-specific widgets from the dynamic right-panel area."""
        for child in self.right_dynamic_frame.winfo_children():
            child.destroy()

    def _update_right_panel_for_active_tab(self):
        """Show fitting controls only on Fitting; otherwise show tab-specific graph range controls."""
        self._clear_right_dynamic_frame()
        active_tab = self._active_middle_tab_text()

        if active_tab == "TCSPC Exponential Fit":
            self._create_tcspc_right_controls(self.right_dynamic_frame)
        elif active_tab == "Fitting":
            self._create_fitting_right_controls(self.right_dynamic_frame)
        else:
            self._create_range_adjuster_right_controls(self.right_dynamic_frame, active_tab)

    def _create_fitting_right_controls(self, parent):
        """Original Kinetic Fitting layout, preserved for the Fitting tab."""
        fit_frame = ttk.LabelFrame(parent, text="Kinetic Fitting", padding="10")
        fit_frame.pack(fill='x', pady=5)

        ttk.Label(fit_frame, text="Fit Mode:").pack()
        self.fit_mode_combo = ttk.Combobox(
            fit_frame,
            values=["single exp", "bi exp", "tri exp"]
        )
        self.fit_mode_combo.current(0)
        self.fit_mode_combo.pack(fill='x', pady=2)

        fit_range_frame = ttk.Frame(fit_frame)
        fit_range_frame.pack(fill='x', pady=5)

        ttk.Label(fit_range_frame, text="Fit Range:").pack(side='left')
        self.fit_start_entry = ttk.Entry(fit_range_frame, width=7)
        self.fit_start_entry.insert(0, "start ps")
        self.fit_start_entry.pack(side='left', padx=2)

        self.fit_end_entry = ttk.Entry(fit_range_frame, width=7)
        self.fit_end_entry.insert(0, "end ns")
        self.fit_end_entry.pack(side='left', padx=2)

        ttk.Button(
            fit_frame, text="Perform Fit", command=self.placeholder_command
        ).pack(fill='x', pady=5)

        ttk.Label(fit_frame, text="Fit Results:").pack()
        self.results_text = tk.Text(fit_frame, height=8, width=30)
        self.results_text.insert(
            '1.0', "τ₁, τ₂, amplitudes, χ²...\n"
        )
        self.results_text.config(state='disabled', bg='#f0f0f0')
        self.results_text.pack(fill='x', expand=True, pady=5)


    # ------------------------------------------------------------------
    # TCSPC exponential fitting interface and calculations
    # ------------------------------------------------------------------
    def _create_tcspc_right_controls(self, parent):
        """Build controls used by the TCSPC Exponential Fit tab."""
        fit_frame = ttk.LabelFrame(parent, text="TCSPC Exponential Fit", padding="12")
        fit_frame.pack(fill='both', expand=True, pady=5)

        settings = ttk.Frame(fit_frame)
        settings.pack(fill='x', pady=(0, 6))
        settings.grid_columnconfigure(1, weight=1)

        ttk.Label(settings, text="Exponentials:").grid(row=0, column=0, sticky='w', pady=2)
        n_exp_spin = tk.Spinbox(
            settings,
            from_=1,
            to=6,
            width=6,
            textvariable=self.tcspc_n_exp_var,
            command=self._tcspc_n_exp_changed,
        )
        n_exp_spin.grid(row=0, column=1, sticky='w', padx=(6, 0), pady=2)
        n_exp_spin.bind("<Return>", self._tcspc_n_exp_changed)
        n_exp_spin.bind("<FocusOut>", self._tcspc_n_exp_changed)

        ttk.Label(settings, text="Fit start (ns):").grid(row=1, column=0, sticky='w', pady=2)
        ttk.Entry(settings, textvariable=self.tcspc_fit_start_var, width=12).grid(
            row=1, column=1, sticky='ew', padx=(6, 0), pady=2
        )
        ttk.Label(settings, text="Fit end (ns):").grid(row=2, column=0, sticky='w', pady=2)
        ttk.Entry(settings, textvariable=self.tcspc_fit_end_var, width=12).grid(
            row=2, column=1, sticky='ew', padx=(6, 0), pady=2
        )
        ttk.Label(
            settings,
            text="Blank start = first non-zero decay; blank end = final point.",
            wraplength=390,
        ).grid(row=3, column=0, columnspan=2, sticky='w', pady=(2, 6))

        try:
            n_exp = self._get_tcspc_n_exponentials()
        except ValueError:
            n_exp = 1
            self.tcspc_n_exp_var.set("1")

        self._ensure_tcspc_parameter_vars(n_exp)

        estimates_frame = ttk.LabelFrame(
            fit_frame, text="Initial Numeric Estimates", padding="8"
        )
        estimates_frame.pack(fill='x', pady=6)
        estimates_frame.grid_columnconfigure(0, weight=2, minsize=135)
        estimates_frame.grid_columnconfigure(1, weight=2, minsize=115)
        estimates_frame.grid_columnconfigure(2, weight=0, minsize=55)
        estimates_frame.grid_columnconfigure(3, weight=1, minsize=75)

        ttk.Label(
            estimates_frame,
            text=(
                "Enter the starting numeric value for every fit parameter below. "
                "These values are passed directly to the optimizer when you click "
                "Perform TCSPC Fit. Check Fixed to hold a value constant."
            ),
            wraplength=400,
            justify='left',
        ).grid(row=0, column=0, columnspan=4, sticky='ew', pady=(0, 8))

        ttk.Label(estimates_frame, text="Parameter").grid(row=1, column=0, sticky='w')
        ttk.Label(estimates_frame, text="Starting value").grid(row=1, column=1, sticky='w')
        ttk.Label(estimates_frame, text="Fixed").grid(row=1, column=2, sticky='w')
        ttk.Label(estimates_frame, text="Allowed").grid(row=1, column=3, sticky='w')

        for row, spec in enumerate(self._tcspc_parameter_specs(n_exp), start=2):
            key, label, allowed = spec
            ttk.Label(estimates_frame, text=label).grid(row=row, column=0, sticky='w', pady=3)
            entry = ttk.Entry(
                estimates_frame,
                textvariable=self.tcspc_estimate_vars[key],
                width=14,
                justify='right',
            )
            entry.grid(row=row, column=1, sticky='ew', padx=(6, 8), pady=3)
            ttk.Checkbutton(
                estimates_frame,
                variable=self.tcspc_fixed_vars[key],
            ).grid(row=row, column=2, sticky='w', padx=(2, 8), pady=3)
            ttk.Label(estimates_frame, text=allowed).grid(
                row=row, column=3, sticky='w', padx=(3, 0), pady=3
            )

        ttk.Button(
            fit_frame,
            text="Reset Estimates",
            command=self._reset_tcspc_initial_estimates,
        ).pack(fill='x', pady=(3, 2))

        ttk.Checkbutton(
            fit_frame,
            text="Generate/print FluoFit-style .dat report",
            variable=self.tcspc_report_var,
        ).pack(anchor='w', pady=(5, 2))

        ttk.Button(
            fit_frame,
            text="Perform TCSPC Fit",
            command=self.perform_tcspc_fit,
        ).pack(fill='x', pady=(6, 5))

        ttk.Label(fit_frame, text="Fit Results:").pack(anchor='w')
        self.tcspc_results_text = tk.Text(fit_frame, height=12, width=45, wrap='word')
        self.tcspc_results_text.insert('1.0', self.tcspc_last_summary)
        self.tcspc_results_text.config(state='disabled', bg='#f0f0f0')
        self.tcspc_results_text.pack(fill='x', expand=True, pady=5)

    def _tcspc_n_exp_changed(self, event=None):
        """Refresh parameter rows after the requested exponential count changes."""
        try:
            self._get_tcspc_n_exponentials()
        except ValueError:
            self.tcspc_n_exp_var.set("1")
        if self._active_middle_tab_text() == "TCSPC Exponential Fit":
            self.after_idle(self._update_right_panel_for_active_tab)

    def _get_tcspc_n_exponentials(self):
        try:
            n_exp = int(self.tcspc_n_exp_var.get().strip())
        except (ValueError, AttributeError):
            raise ValueError("Number of exponentials must be a whole number.")
        if n_exp < 1 or n_exp > 6:
            raise ValueError("Number of exponentials must be between 1 and 6.")
        return n_exp

    def _tcspc_parameter_specs(self, n_exponentials):
        specs = []
        for idx in range(1, n_exponentials + 1):
            specs.append((f"A{idx}", f"A{idx} [Cnts]", ">= 0"))
            specs.append((f"tau{idx}", f"tau{idx} [ns]", "0.01-100"))
        specs.extend([
            ("shift_irf_ns", "Shift IRF [ns]", "-5 to 5"),
            ("background_decay", "Decay bg [Cnts]", ">= 0"),
            ("background_irf", "IRF bg [Cnts]", ">= 0"),
        ])
        return specs

    def _tcspc_parameter_keys(self, n_exponentials):
        return [spec[0] for spec in self._tcspc_parameter_specs(n_exponentials)]

    def _ensure_tcspc_parameter_vars(self, n_exponentials):
        defaults = self._tcspc_default_initial_params(n_exponentials)
        keys = self._tcspc_parameter_keys(n_exponentials)
        for key, default in zip(keys, defaults):
            if key not in self.tcspc_estimate_vars:
                self.tcspc_estimate_vars[key] = tk.StringVar(value=f"{default:.6g}")
            if key not in self.tcspc_fixed_vars:
                self.tcspc_fixed_vars[key] = tk.BooleanVar(value=False)

    def _reset_tcspc_initial_estimates(self):
        try:
            n_exp = self._get_tcspc_n_exponentials()
        except ValueError as exc:
            messagebox.showerror("Invalid exponentials", str(exc), parent=self)
            return
        defaults = self._tcspc_default_initial_params(n_exp)
        self._ensure_tcspc_parameter_vars(n_exp)
        for key, value in zip(self._tcspc_parameter_keys(n_exp), defaults):
            self.tcspc_estimate_vars[key].set(f"{value:.6g}")
            self.tcspc_fixed_vars[key].set(False)
        if self._active_middle_tab_text() == "TCSPC Exponential Fit":
            self._update_right_panel_for_active_tab()

    def _tcspc_default_initial_params(self, n_exponentials):
        if self.tcspc_decay_raw is not None:
            try:
                decay, _ = self._tcspc_divide_nonzero_by_first_nonzero(
                    self.tcspc_decay_raw, "Decay"
                )
            except ValueError:
                decay = np.asarray([1.0])
        else:
            decay = np.asarray([1.0])
        return self._tcspc_make_initial_params(decay, n_exponentials)

    def import_tcspc_decay(self):
        self._import_tcspc_curve("decay")

    def import_tcspc_irf(self):
        self._import_tcspc_curve("irf")

    def _import_tcspc_curve(self, curve_type):
        """Load an SDT decay or IRF into the fitting tab and immediately preview it."""
        if SdtFile is None:
            messagebox.showerror(
                "Missing dependency",
                "TCSPC loading requires the 'sdtfile' package. Install it with: pip install sdtfile",
                parent=self,
            )
            return

        is_decay = curve_type == "decay"
        title = "Select TCSPC lifetime decay SDT file" if is_decay else "Select TCSPC IRF SDT file"
        path = filedialog.askopenfilename(
            title=title,
            filetypes=[("SDT files", "*.sdt"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            time_ns, values = self._tcspc_load_sdt_decay(path, block=0)
            self._tcspc_divide_nonzero_by_first_nonzero(
                values, "Decay" if is_decay else "IRF"
            )
        except Exception as exc:
            messagebox.showerror("TCSPC import error", str(exc), parent=self)
            return

        if is_decay:
            self.tcspc_decay_path = path
            self.tcspc_decay_time_ns = time_ns
            self.tcspc_decay_raw = values
            self.tcspc_decay_file_label.config(text=f"Decay: {Path(path).name}")
            # Update generic amplitude guesses to the scale of the newly loaded decay.
            self._reset_tcspc_initial_estimates()
        else:
            self.tcspc_irf_path = path
            self.tcspc_irf_time_ns = time_ns
            self.tcspc_irf_raw = values
            self.tcspc_irf_file_label.config(text=f"IRF: {Path(path).name}")

        self.tcspc_fit_result = None
        self.tcspc_last_summary = "Data loaded. Configure estimates on the right, then perform the fit."
        self._set_tcspc_results_text(self.tcspc_last_summary)
        self._plot_tcspc_preview()
        self.status_label.config(
            text=f"Status: TCSPC {'decay' if is_decay else 'IRF'} loaded."
        )

    def _draw_tcspc_empty_plot(self):
        self.tcspc_fig.clear()
        self.tcspc_ax = self.tcspc_fig.add_subplot(111)
        self.tcspc_ax.set_title("TCSPC lifetime decay / IRF preview")
        self.tcspc_ax.set_xlabel("Time [ns]")
        self.tcspc_ax.set_ylabel("Normalized intensity")
        self.tcspc_ax.grid(True)
        self.tcspc_ax.text(
            0.5, 0.5,
            "Import a TCSPC lifetime decay and IRF above",
            ha='center', va='center', transform=self.tcspc_ax.transAxes,
            fontsize=11, color='gray', alpha=0.6,
        )
        self.tcspc_fig.tight_layout()

    def _plot_tcspc_preview(self):
        """Plot loaded data before fitting; residuals are added only after a fit."""
        self.tcspc_fig.clear()
        ax = self.tcspc_fig.add_subplot(111)
        self.tcspc_ax = ax

        plotted = False
        if self.tcspc_decay_raw is not None:
            decay, _ = self._tcspc_divide_nonzero_by_first_nonzero(
                self.tcspc_decay_raw, "Decay"
            )
            plotted |= self._plot_tcspc_normalized_trace(
                ax,
                self.tcspc_decay_time_ns,
                decay,
                '.', markersize=2.5, label="Lifetime decay",
            )

        if self.tcspc_irf_raw is not None:
            irf, _ = self._tcspc_divide_nonzero_by_first_nonzero(
                self.tcspc_irf_raw, "IRF"
            )
            plotted |= self._plot_tcspc_normalized_trace(
                ax,
                self.tcspc_irf_time_ns,
                irf,
                '-', linewidth=1.0, label="IRF",
            )

        ax.set_title("TCSPC data before fit")
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel("Normalized intensity (≥ 1)")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1.0)
        ax.grid(True)
        if plotted:
            ax.legend(fontsize=8)
        else:
            ax.text(
                0.5, 0.5, "Import TCSPC data above",
                ha='center', va='center', transform=ax.transAxes,
                color='gray', alpha=0.6,
            )
        self.tcspc_fig.tight_layout()
        self.tcspc_canvas.draw()

    def _plot_tcspc_fit_and_residuals(self, fit_result):
        """Add the fitted curve/chi-square and a residual plot below the main graph."""
        self.tcspc_fig.clear()
        grid = self.tcspc_fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
        ax_fit = self.tcspc_fig.add_subplot(grid[0, 0])
        ax_res = self.tcspc_fig.add_subplot(grid[1, 0], sharex=ax_fit)
        self.tcspc_ax = ax_fit
        self.tcspc_residual_ax = ax_res

        time_ns = fit_result["time_ns"]
        decay = fit_result["decay"]
        irf = fit_result["irf"]
        model = fit_result["fit"]
        residuals = fit_result["residuals"]

        self._plot_tcspc_normalized_trace(
            ax_fit, time_ns, decay, '.', markersize=2.5, label="Decay"
        )
        self._plot_tcspc_normalized_trace(
            ax_fit, time_ns, irf, '-', linewidth=1.0, label="IRF"
        )
        self._plot_tcspc_normalized_trace(
            ax_fit, time_ns, model, '-', linewidth=1.3, label="Fitted model"
        )
        ax_fit.set_ylabel("Normalized intensity (≥ 1)")
        ax_fit.set_yscale("log")
        ax_fit.set_ylim(bottom=1.0)
        ax_fit.set_title(
            f"{fit_result['n_exponentials']}-exponential reconvolution fit | "
            f"Reduced chi-square = {fit_result['reduced_chi2']:.5f}"
        )
        ax_fit.grid(True)
        ax_fit.legend(fontsize=8)
        ax_fit.tick_params(labelbottom=False)

        parameter_lines = []
        for idx, (amplitude, tau_ns) in enumerate(
            zip(fit_result["amplitudes_counts"], fit_result["taus_ns"]), start=1
        ):
            parameter_lines.append(f"A{idx}={amplitude:.3g}, tau{idx}={tau_ns:.4g} ns")
        ax_fit.text(
            0.98, 0.97, "\n".join(parameter_lines),
            transform=ax_fit.transAxes, ha='right', va='top', fontsize=8,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
        )

        ax_res.plot(time_ns, residuals, linewidth=0.8)
        ax_res.axhline(0, linestyle='--', linewidth=0.8)
        ax_res.set_xlabel("Time [ns]")
        ax_res.set_ylabel("Residual")
        ax_res.grid(True)

        self.tcspc_fig.tight_layout()
        self.tcspc_canvas.draw()

    def perform_tcspc_fit(self):
        """Validate right-panel settings, perform the reconvolution fit, and update plots."""
        if least_squares is None:
            messagebox.showerror(
                "Missing dependency",
                "TCSPC fitting requires SciPy. Install it with: pip install scipy",
                parent=self,
            )
            return
        if self.tcspc_decay_raw is None or self.tcspc_irf_raw is None:
            messagebox.showwarning(
                "Missing TCSPC data",
                "Import both the lifetime decay and the IRF before fitting.",
                parent=self,
            )
            return

        try:
            n_exp = self._get_tcspc_n_exponentials()
            self._ensure_tcspc_parameter_vars(n_exp)
            initial_params, fixed_params = self._read_tcspc_parameter_controls(n_exp)
            fit_start = self._parse_optional_float(self.tcspc_fit_start_var.get(), "Fit start")
            fit_end = self._parse_optional_float(self.tcspc_fit_end_var.get(), "Fit end")
            if fit_start is not None and fit_end is not None and fit_start > fit_end:
                raise ValueError("Fit start must be less than or equal to fit end.")
        except ValueError as exc:
            messagebox.showerror("TCSPC fit settings", str(exc), parent=self)
            return

        self.status_label.config(text="Status: Running TCSPC exponential fit...")
        self.update_idletasks()

        try:
            fit_result = self._fit_tcspc_multi_exp_reconv(
                n_exponentials=n_exp,
                fit_start_ns=fit_start,
                fit_end_ns=fit_end,
                initial_params=initial_params,
                fixed_params=fixed_params,
            )
            self.tcspc_fit_result = fit_result
            self.tcspc_last_summary = self._format_tcspc_fit_summary(fit_result)
            self._set_tcspc_results_text(self.tcspc_last_summary)
            self._plot_tcspc_fit_and_residuals(fit_result)

            report_path = None
            if self.tcspc_report_var.get():
                report_path = filedialog.asksaveasfilename(
                    title="Save TCSPC fit report",
                    defaultextension=".dat",
                    filetypes=[("DAT files", "*.dat"), ("All files", "*.*")],
                    initialfile=(
                        Path(self.tcspc_decay_path).stem
                        + f"_{n_exp}exp_fit.dat"
                    ),
                )
                if report_path:
                    self._write_tcspc_fluofit_style_dat(fit_result, report_path)

            status = f"Status: TCSPC fit complete (reduced chi-square {fit_result['reduced_chi2']:.5f})."
            if report_path:
                status += f" Report saved: {Path(report_path).name}."
            self.status_label.config(text=status)

        except Exception as exc:
            messagebox.showerror("TCSPC fit error", str(exc), parent=self)
            self.status_label.config(text="Status: TCSPC fit failed.")

    def _set_tcspc_results_text(self, text):
        widget = getattr(self, "tcspc_results_text", None)
        if widget is None or not widget.winfo_exists():
            return
        widget.config(state='normal')
        widget.delete('1.0', tk.END)
        widget.insert('1.0', text)
        widget.config(state='disabled')

    def _parse_optional_float(self, text, label):
        text = str(text).strip()
        if not text:
            return None
        try:
            value = float(text)
        except ValueError:
            raise ValueError(f"{label} must be numeric or blank.")
        if not np.isfinite(value):
            raise ValueError(f"{label} must be finite.")
        return value

    def _read_tcspc_parameter_controls(self, n_exponentials):
        keys = self._tcspc_parameter_keys(n_exponentials)
        values = []
        fixed = []
        for key in keys:
            text = self.tcspc_estimate_vars[key].get().strip()
            if not text:
                raise ValueError(f"Initial estimate for {key} is blank.")
            try:
                value = float(text)
            except ValueError:
                raise ValueError(f"Initial estimate for {key} must be numeric.")
            if not np.isfinite(value):
                raise ValueError(f"Initial estimate for {key} must be finite.")
            values.append(value)
            fixed.append(bool(self.tcspc_fixed_vars[key].get()))

        values = np.asarray(values, dtype=float)
        lower, upper = self._tcspc_make_bounds(n_exponentials)
        lower = np.asarray(lower, dtype=float)
        upper = np.asarray(upper, dtype=float)
        if np.any(values < lower) or np.any(values > upper):
            raise ValueError(
                "One or more initial estimates are outside the allowed ranges shown on the right."
            )
        return values, np.asarray(fixed, dtype=bool)

    def _tcspc_load_sdt_decay(self, path, block=0):
        """Load a decay and time axis from an SDT file, converting seconds to ns."""
        if SdtFile is None:
            raise RuntimeError("The sdtfile package is not installed.")
        with SdtFile(path) as sdt:
            data = np.asarray(sdt.data[block], dtype=float)
            time = np.asarray(sdt.times[block], dtype=float)

        if data.ndim == 1:
            decay = data
        else:
            decay = data.reshape(-1, data.shape[-1]).sum(axis=0)

        time = np.squeeze(time)
        if time.ndim != 1:
            raise ValueError("SDT time axis could not be reduced to one dimension.")
        if decay.size != time.size:
            raise ValueError("SDT decay length does not match its time axis.")
        if time.size < 2:
            raise ValueError("SDT data needs at least two time points.")

        if np.nanmax(time) < 1e-3:
            time_ns = time * 1e9
        else:
            time_ns = time
        return np.asarray(time_ns, dtype=float), np.asarray(decay, dtype=float)

    def _tcspc_shift_curve(self, values, shift_bins):
        x = np.arange(len(values), dtype=float)
        return np.interp(x - shift_bins, x, values, left=0.0, right=0.0)

    def _tcspc_unpack_multi_exp_params(self, params, n_exponentials):
        params = np.asarray(params, dtype=float)
        amplitudes = params[0 : 2 * n_exponentials : 2]
        taus_ns = params[1 : 2 * n_exponentials : 2]
        shift_irf_ns = params[2 * n_exponentials]
        background_decay = params[2 * n_exponentials + 1]
        background_irf = params[2 * n_exponentials + 2]
        return amplitudes, taus_ns, shift_irf_ns, background_decay, background_irf

    def _tcspc_reconvolved_multi_exp(self, params, time_ns, irf, dt, n_exponentials):
        amplitudes, taus_ns, shift_irf_ns, background_decay, background_irf = (
            self._tcspc_unpack_multi_exp_params(params, n_exponentials)
        )
        shift_bins = shift_irf_ns / dt
        shifted_irf = self._tcspc_shift_curve(irf + background_irf, shift_bins)

        t0 = time_ns[0]
        time_from_start = time_ns - t0
        pure_decay = np.zeros_like(time_ns, dtype=float)
        for amplitude, tau_ns in zip(amplitudes, taus_ns):
            pure_decay += amplitude * np.exp(-time_from_start / tau_ns)
        pure_decay[time_ns < t0] = 0.0

        model = np.convolve(shifted_irf, pure_decay, mode="full")[: len(time_ns)]
        irf_area = shifted_irf.sum()
        if irf_area > 0:
            model = model / irf_area
        return model + background_decay

    def _tcspc_poisson_residuals(self, params, time_ns, decay, irf, dt, fit_mask, n_exponentials):
        model = self._tcspc_reconvolved_multi_exp(params, time_ns, irf, dt, n_exponentials)
        sigma = np.sqrt(np.maximum(decay, 1.0))
        return (decay[fit_mask] - model[fit_mask]) / sigma[fit_mask]

    def _tcspc_expand_free_params(self, free_params, base_params, fixed_params):
        full_params = np.asarray(base_params, dtype=float).copy()
        fixed_params = np.asarray(fixed_params, dtype=bool)
        full_params[~fixed_params] = np.asarray(free_params, dtype=float)
        return full_params

    def _tcspc_masked_poisson_residuals(
        self, free_params, base_params, fixed_params,
        time_ns, decay, irf, dt, fit_mask, n_exponentials,
    ):
        full_params = self._tcspc_expand_free_params(free_params, base_params, fixed_params)
        return self._tcspc_poisson_residuals(
            full_params, time_ns, decay, irf, dt, fit_mask, n_exponentials
        )

    def _tcspc_make_initial_params(self, decay, n_exponentials):
        max_decay = float(np.nanmax(decay)) if np.size(decay) else 1.0
        if not np.isfinite(max_decay) or max_decay <= 0:
            max_decay = 1.0
        amplitude_guess = max_decay / max(n_exponentials, 1)
        if n_exponentials == 1:
            tau_guesses = np.asarray([5.5], dtype=float)
        else:
            tau_guesses = np.geomspace(0.5, 12.0, n_exponentials)

        params = []
        for tau_guess in tau_guesses:
            params.extend([amplitude_guess, float(tau_guess)])
        params.extend([-0.04, 0.0, 0.0])
        return np.asarray(params, dtype=float)

    def _tcspc_make_bounds(self, n_exponentials):
        lower_bounds = []
        upper_bounds = []
        for _ in range(n_exponentials):
            lower_bounds.extend([0.0, 0.01])
            upper_bounds.extend([np.inf, 100.0])
        lower_bounds.extend([-5.0, 0.0, 0.0])
        upper_bounds.extend([5.0, np.inf, np.inf])
        return lower_bounds, upper_bounds

    def _tcspc_divide_nonzero_by_first_nonzero(self, values, data_name="data"):
        values = np.asarray(values, dtype=float)
        nonzero_indices = np.flatnonzero(values != 0)
        if nonzero_indices.size == 0:
            raise ValueError(f"{data_name} contains no non-zero values.")
        first_nonzero = float(values[nonzero_indices[0]])
        normalized = np.where(values != 0, values / first_nonzero, values)
        return normalized, first_nonzero

    def _plot_tcspc_normalized_trace(self, ax, time_ns, intensity, *args, **kwargs):
        """Plot only finite TCSPC points at or above the log-scale minimum."""
        time_ns = np.asarray(time_ns, dtype=float)
        intensity = np.asarray(intensity, dtype=float)
        mask = np.isfinite(time_ns) & np.isfinite(intensity) & (intensity >= 1.0)
        if not np.any(mask):
            return False
        ax.plot(time_ns[mask], intensity[mask], *args, **kwargs)
        return True

    def _fit_tcspc_multi_exp_reconv(
        self,
        n_exponentials=1,
        fit_start_ns=None,
        fit_end_ns=None,
        initial_params=None,
        fixed_params=None,
    ):
        n_exponentials = int(n_exponentials)
        if n_exponentials < 1:
            raise ValueError("Number of exponentials must be at least 1.")

        time_ns = np.asarray(self.tcspc_decay_time_ns, dtype=float).copy()
        decay = np.asarray(self.tcspc_decay_raw, dtype=float).copy()
        irf_time_ns = np.asarray(self.tcspc_irf_time_ns, dtype=float).copy()
        irf = np.asarray(self.tcspc_irf_raw, dtype=float).copy()

        decay, decay_divisor = self._tcspc_divide_nonzero_by_first_nonzero(decay, "Decay")
        irf, irf_divisor = self._tcspc_divide_nonzero_by_first_nonzero(irf, "IRF")

        if len(irf_time_ns) != len(time_ns) or not np.allclose(irf_time_ns, time_ns):
            irf = np.interp(time_ns, irf_time_ns, irf, left=0.0, right=0.0)

        dt = float(np.median(np.diff(time_ns)))
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("TCSPC time axis must be strictly increasing with a positive spacing.")

        nonzero_decay_indices = np.flatnonzero(decay != 0)
        if nonzero_decay_indices.size == 0:
            raise ValueError("Decay contains no non-zero values, so a fit start cannot be determined.")
        first_nonzero_idx = int(nonzero_decay_indices[0])
        first_nonzero_ns = float(time_ns[first_nonzero_idx])

        if fit_start_ns is None:
            fit_start_ns = first_nonzero_ns
        if fit_end_ns is None:
            fit_end_ns = float(time_ns[-1])

        # Ignore normalized values below 10^0 so zero/low points do not enter the
        # logarithmic-display range or the nonlinear fit.
        fit_mask = (
            (time_ns >= fit_start_ns)
            & (time_ns <= fit_end_ns)
            & np.isfinite(decay)
            & (decay >= 1.0)
        )
        chi_square_mask = fit_mask & (time_ns >= first_nonzero_ns)
        if not np.any(fit_mask):
            raise ValueError("The selected fit range contains no TCSPC data points.")

        if initial_params is None:
            p0 = self._tcspc_make_initial_params(decay, n_exponentials)
        else:
            p0 = np.asarray(initial_params, dtype=float)
            expected = 2 * n_exponentials + 3
            if p0.size != expected:
                raise ValueError(
                    f"Expected {expected} initial parameters for {n_exponentials} exponentials, "
                    f"but received {p0.size}."
                )

        lower_bounds, upper_bounds = self._tcspc_make_bounds(n_exponentials)
        lower_bounds = np.asarray(lower_bounds, dtype=float)
        upper_bounds = np.asarray(upper_bounds, dtype=float)
        if np.any(p0 < lower_bounds) or np.any(p0 > upper_bounds):
            raise ValueError("One or more initial estimates are outside the allowed bounds.")

        if fixed_params is None:
            fixed_params = np.zeros_like(p0, dtype=bool)
        else:
            fixed_params = np.asarray(fixed_params, dtype=bool)
            if fixed_params.size != p0.size:
                raise ValueError(
                    f"Expected {p0.size} fixed/free flags, but received {fixed_params.size}."
                )

        free_mask = ~fixed_params
        free_p0 = p0[free_mask]
        free_lower = lower_bounds[free_mask]
        free_upper = upper_bounds[free_mask]

        if free_p0.size == 0:
            params = p0.copy()
            fit_success = True
            fit_message = "All parameters were fixed; no nonlinear optimization was run."
        else:
            result = least_squares(
                self._tcspc_masked_poisson_residuals,
                free_p0,
                bounds=(free_lower, free_upper),
                args=(p0, fixed_params, time_ns, decay, irf, dt, fit_mask, n_exponentials),
                loss="linear",
                max_nfev=2000,
                ftol=1e-6,
                xtol=1e-6,
            )
            params = self._tcspc_expand_free_params(result.x, p0, fixed_params)
            fit_success = bool(result.success)
            fit_message = str(result.message)

        amplitudes, taus_ns, shift_irf_ns, bg_decay, bg_irf = (
            self._tcspc_unpack_multi_exp_params(params, n_exponentials)
        )
        model = self._tcspc_reconvolved_multi_exp(params, time_ns, irf, dt, n_exponentials)

        full_mask = np.ones_like(time_ns, dtype=bool)
        residuals = self._tcspc_poisson_residuals(
            params, time_ns, decay, irf, dt, full_mask, n_exponentials
        )
        chi_residuals = self._tcspc_poisson_residuals(
            params, time_ns, decay, irf, dt, chi_square_mask, n_exponentials
        )
        n_fit_points = int(chi_square_mask.sum())
        n_params = int(np.sum(~fixed_params))
        reduced_chi2 = np.sum(chi_residuals ** 2) / max(n_fit_points - n_params, 1)

        amplitude_sum = float(np.sum(amplitudes))
        amplitude_fractions = (
            amplitudes / amplitude_sum if amplitude_sum > 0 else np.zeros_like(amplitudes)
        )
        intensity_weights_raw = amplitudes * taus_ns
        intensity_sum = float(np.sum(intensity_weights_raw))
        intensity_fractions = (
            intensity_weights_raw / intensity_sum
            if intensity_sum > 0 else np.zeros_like(intensity_weights_raw)
        )
        amplitude_weighted_lifetime = (
            float(np.sum(amplitudes * taus_ns) / amplitude_sum)
            if amplitude_sum > 0 else np.nan
        )
        intensity_weighted_lifetime = (
            float(np.sum(amplitudes * taus_ns ** 2) / np.sum(amplitudes * taus_ns))
            if np.sum(amplitudes * taus_ns) > 0 else np.nan
        )

        return {
            "time_ns": time_ns,
            "decay": decay,
            "irf": irf,
            "fit": model,
            "residuals": residuals,
            "fit_mask": fit_mask,
            "n_exponentials": n_exponentials,
            "initial_params": p0.copy(),
            "fixed_params": fixed_params.copy(),
            "free_params_count": int(np.sum(~fixed_params)),
            "amplitudes_counts": amplitudes,
            "taus_ns": taus_ns,
            "amplitude_fractions": amplitude_fractions,
            "intensity_fractions": intensity_fractions,
            "amplitude_weighted_lifetime_ns": amplitude_weighted_lifetime,
            "intensity_weighted_lifetime_ns": intensity_weighted_lifetime,
            "shift_irf_ns": shift_irf_ns,
            "background_decay_counts": bg_decay,
            "background_irf_counts": bg_irf,
            "reduced_chi2": float(reduced_chi2),
            "fitted_points": n_fit_points,
            "success": fit_success,
            "fit_message": fit_message,
            "decay_divisor": decay_divisor,
            "irf_divisor": irf_divisor,
        }

    def _format_tcspc_fit_summary(self, fit_result):
        lines = [
            f"{fit_result['n_exponentials']}-exponential fit",
            f"Reduced chi-square: {fit_result['reduced_chi2']:.5f}",
            f"Fitted points: {fit_result['fitted_points']}",
            f"Free parameters: {fit_result['free_params_count']}",
            "",
        ]
        for idx, (amplitude, tau_ns, amp_frac, int_frac) in enumerate(
            zip(
                fit_result["amplitudes_counts"],
                fit_result["taus_ns"],
                fit_result["amplitude_fractions"],
                fit_result["intensity_fractions"],
            ), start=1,
        ):
            lines.extend([
                f"A{idx}: {amplitude:.6g}",
                f"tau{idx}: {tau_ns:.6g} ns",
                f"Amplitude fraction: {100.0 * amp_frac:.2f}%",
                f"Intensity fraction: {100.0 * int_frac:.2f}%",
            ])
        lines.extend([
            "",
            f"IRF shift: {fit_result['shift_irf_ns']:.6g} ns",
            f"Decay background: {fit_result['background_decay_counts']:.6g}",
            f"IRF background: {fit_result['background_irf_counts']:.6g}",
            f"Amplitude-weighted lifetime: {fit_result['amplitude_weighted_lifetime_ns']:.6g} ns",
            f"Intensity-weighted lifetime: {fit_result['intensity_weighted_lifetime_ns']:.6g} ns",
        ])
        return "\n".join(lines)

    def _write_tcspc_fluofit_style_dat(self, fit_result, output_dat_path, title=None):
        """Write the integrated fit result using the standalone script's FluoFit-style layout."""
        time_ns = np.asarray(fit_result["time_ns"])
        decay = np.asarray(fit_result["decay"])
        irf = np.asarray(fit_result["irf"])
        model = np.asarray(fit_result["fit"])
        residuals = np.asarray(fit_result["residuals"])
        n_exp = int(fit_result["n_exponentials"])
        amplitudes = np.asarray(fit_result["amplitudes_counts"])
        taus_ns = np.asarray(fit_result["taus_ns"])
        amplitude_fractions = np.asarray(fit_result["amplitude_fractions"])
        intensity_fractions = np.asarray(fit_result["intensity_fractions"])
        fixed_params = np.asarray(fit_result.get("fixed_params", []), dtype=bool)

        decay_name = os.path.basename(self.tcspc_decay_path or "decay.sdt")
        irf_name = os.path.basename(self.tcspc_irf_path or "irf.sdt")
        if title is None:
            title = f"Python {n_exp}-exponential reconvolution fit"

        lines = [
            "PicoQuant FluoFit-style Python Export",
            f"Saved : {datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')}",
            "",
            title,
            "",
            f"Model: Exp. [Reconv.] ({n_exp} Exponential Components)",
            "Plotted Data Set #0 Decay:",
            f'"{decay_name}" (0)',
            "Plotted Data Set #0 IRF:",
            f'"{irf_name}" (0)',
            f"X^2(reduced): {fit_result['reduced_chi2']:.4f} ; Fitted Data Points: {fit_result['fitted_points']}",
            "",
            "Data Set #0",
            "     Decay                 IRF                   Model Decay           Residuals",
            "     t[ns]     Intens.     t[ns]     Intens.     t[ns]     Intens.     t[ns]       diff.",
        ]

        for t, y, h, m, r in zip(time_ns, decay, irf, model, residuals):
            lines.append(
                f"{t:12.6f}{y:12.6f}{t:12.6f}{h:12.6f}"
                f"{t:12.6f}{m:12.6f}{t:12.6f}{r:12.6f}"
            )

        lines.extend([
            "",
            "Best Fit Parameters:",
            "",
            "Data Set #0",
            "Parameter              Value          Conf. Lower     Conf. Upper     Conf. Estimation",
        ])
        for idx, (amplitude, tau_ns) in enumerate(zip(amplitudes, taus_ns), start=1):
            amp_fixed = fixed_params.size and fixed_params[2 * (idx - 1)]
            tau_fixed = fixed_params.size and fixed_params[2 * (idx - 1) + 1]
            lines.append(
                f"A{idx} [Cnts]        {amplitude:14.5f}              ---            ---     "
                f"{'Fixed' if amp_fixed else 'Fitting'}"
            )
            lines.append(
                f"t{idx} [ns]          {tau_ns:14.5f}              ---            ---     "
                f"{'Fixed' if tau_fixed else 'Fitting'}"
            )

        shift_idx = 2 * n_exp
        lines.append(
            f"Bkgr. Dec [Cnts] {fit_result['background_decay_counts']:14.5f}              ---            ---     "
            f"{'Fixed' if fixed_params.size and fixed_params[shift_idx + 1] else 'Fitting'}"
        )
        lines.append(
            f"Bkgr. IRF [Cnts] {fit_result['background_irf_counts']:14.5f}              ---            ---     "
            f"{'Fixed' if fixed_params.size and fixed_params[shift_idx + 2] else 'Fitting'}"
        )
        lines.append(
            f"Shift IRF [ns]   {fit_result['shift_irf_ns']:14.5f}              ---            ---     "
            f"{'Fixed' if fixed_params.size and fixed_params[shift_idx] else 'Fitting'}"
        )

        lines.extend([
            "",
            "Average Lifetime:",
            f"tAv.1={fit_result['intensity_weighted_lifetime_ns']:.5f} ns (intensity weighted)",
            f"tAv.2={fit_result['amplitude_weighted_lifetime_ns']:.5f} ns (amplitude weighted)",
            "",
            "Fractional Intensities of the Positive Decay Components:",
        ])
        for idx, (tau_ns, fraction) in enumerate(zip(taus_ns, intensity_fractions), start=1):
            lines.append(f"t{idx} ({tau_ns:.5f} ns) : {100.0 * fraction:.2f}%")
        lines.extend(["", "Fractional Amplitudes of the Positive Decay Components:"])
        for idx, (tau_ns, fraction) in enumerate(zip(taus_ns, amplitude_fractions), start=1):
            lines.append(f"t{idx} ({tau_ns:.5f} ns) : {100.0 * fraction:.2f}%")

        lines.extend(["", "Fitted Decay and Exponential Components:"])
        component_header = "     Time [ns]" + "".join(
            f"      t{idx} [kCounts]" for idx in range(1, n_exp + 1)
        ) + "      Sum [kCounts]"
        lines.append(component_header)
        component_time = np.linspace(0.0, 22.0, 1000)
        components = [
            amplitude * np.exp(-component_time / tau_ns)
            for amplitude, tau_ns in zip(amplitudes, taus_ns)
        ]
        component_sum = np.sum(components, axis=0)
        for row_idx, t in enumerate(component_time):
            row = f"{t:12.5f}"
            for component in components:
                row += f"{component[row_idx] / 1000.0:16.8f}"
            row += f"{component_sum[row_idx] / 1000.0:16.8f}"
            lines.append(row)

        lines.extend(["", "Confidence Intervals:"])
        for idx in range(1, n_exp + 1):
            lines.append(f"A{idx} : [--- ; ---] Cnts")
            lines.append(f"t{idx} : [--- ; ---] ns")
        lines.append("Shift IRF : [--- ; ---] ns")
        lines.append("")

        with open(output_dat_path, "w", encoding="latin-1", newline="") as handle:
            handle.write("\r\n".join(lines))
        return output_dat_path

    def _create_range_adjuster_right_controls(self, parent, active_tab):
        """Create tab-specific numeric range inputs and a Re-graph button."""
        range_frame = ttk.LabelFrame(parent, text=f"{active_tab} Range Options", padding="10")
        range_frame.pack(fill='x', pady=5)

        ttk.Label(
            range_frame,
            text="Enter numeric Min/Max ranges, then click Re-graph to redraw the active graph with those limits.",
            wraplength=220
        ).pack(fill='x', pady=(0, 8))

        if not hasattr(self, "right_range_vars"):
            self.right_range_vars = {}
        if not hasattr(self, "right_range_values"):
            self.right_range_values = {}

        if active_tab == "Time-Domain Traces":
            range_specs = [
                ("Wavelength", "nm", self._safe_axis_min(getattr(self, "wavelength", None), 400), self._safe_axis_max(getattr(self, "wavelength", None), 800)),
                ("Signal", "ΔA", self._safe_data_min(-1.0), self._safe_data_max(1.0)),
            ]
        elif active_tab == "Wavelength-Domain Traces":
            range_specs = [
                ("Time", "ps", self._safe_axis_min(getattr(self, "times", None), 0), self._safe_axis_max(getattr(self, "times", None), 1000)),
                ("Signal", "ΔA", self._safe_data_min(-1.0), self._safe_data_max(1.0)),
            ]
        elif active_tab == "2D Map View":
            range_specs = [
                ("Wavelength", "nm", self._safe_axis_min(getattr(self, "wavelength", None), 400), self._safe_axis_max(getattr(self, "wavelength", None), 800)),
                ("Time", "ps", self._safe_positive_axis_min(getattr(self, "times", None), 1), self._safe_axis_max(getattr(self, "times", None), 1000)),
                ("Color", "ΔA", self._safe_data_min(-1.0), self._safe_data_max(1.0)),
            ]
        else:
            range_specs = []

        self.right_range_vars[active_tab] = {}
        saved_tab_values = self.right_range_values.get(active_tab, {})

        for label, units, default_min, default_max in range_specs:
            saved_min, saved_max = saved_tab_values.get(label, (default_min, default_max))
            min_var, max_var = self._add_numeric_range_inputs(
                range_frame,
                label=label,
                units=units,
                min_value=saved_min,
                max_value=saved_max
            )
            self.right_range_vars[active_tab][label] = {
                "min": min_var,
                "max": max_var,
                "units": units,
            }

        ttk.Button(
            range_frame,
            text="Re-graph",
            command=lambda tab=active_tab: self._right_panel_regraph_requested(tab)
        ).pack(fill='x', pady=(10, 2))

    def _add_numeric_range_inputs(self, parent, label, units, min_value, max_value):
        """Add one min/max numeric input row for a graph range."""
        row = ttk.Frame(parent)
        row.pack(fill='x', pady=4)
        row.grid_columnconfigure(1, weight=1)
        row.grid_columnconfigure(3, weight=1)

        min_var = tk.StringVar(value=f"{float(min_value):.4g}")
        max_var = tk.StringVar(value=f"{float(max_value):.4g}")

        ttk.Label(row, text=f"{label} ({units})").grid(row=0, column=0, columnspan=4, sticky='w')
        ttk.Label(row, text="Min").grid(row=1, column=0, sticky='w', padx=(0, 3), pady=2)
        ttk.Entry(row, textvariable=min_var, width=9).grid(row=1, column=1, sticky='ew', padx=(0, 6), pady=2)
        ttk.Label(row, text="Max").grid(row=1, column=2, sticky='w', padx=(0, 3), pady=2)
        ttk.Entry(row, textvariable=max_var, width=9).grid(row=1, column=3, sticky='ew', pady=2)

        return min_var, max_var

    def _right_panel_regraph_requested(self, active_tab):
        """Validate range fields and redraw the active graph with those ranges."""
        parsed_ranges = self._parse_right_panel_ranges(active_tab)
        if parsed_ranges is None:
            return

        if active_tab == "Time-Domain Traces":
            self._regraph_time_domain_from_right_panel(parsed_ranges)
        elif active_tab == "Wavelength-Domain Traces":
            self._regraph_wavelength_domain_from_right_panel(parsed_ranges)
        elif active_tab == "2D Map View":
            self._regraph_2d_map_from_right_panel(parsed_ranges)

    def _parse_right_panel_ranges(self, active_tab):
        """Return validated right-panel ranges as {label: (min, max, units)}."""
        tab_ranges = getattr(self, "right_range_vars", {}).get(active_tab, {})
        parsed_ranges = {}

        for label, spec in tab_ranges.items():
            try:
                min_value = float(spec["min"].get().strip())
                max_value = float(spec["max"].get().strip())
            except ValueError:
                messagebox.showerror(
                    "Input Error",
                    f"Invalid {label} range. Enter numeric Min and Max values."
                )
                return None

            if min_value == max_value:
                messagebox.showerror(
                    "Input Error",
                    f"Invalid {label} range. Min and Max cannot be the same value."
                )
                return None

            if min_value > max_value:
                min_value, max_value = max_value, min_value
                spec["min"].set(f"{min_value:.4g}")
                spec["max"].set(f"{max_value:.4g}")

            parsed_ranges[label] = (min_value, max_value, spec["units"])

        if not hasattr(self, "right_range_values"):
            self.right_range_values = {}
        self.right_range_values[active_tab] = {
            label: (values[0], values[1])
            for label, values in parsed_ranges.items()
        }

        return parsed_ranges

    def _set_axis_limits_from_right_ranges(self, ax, parsed_ranges, x_label, y_label=None):
        """Apply parsed right-panel x/y limits to an axes object."""
        if x_label in parsed_ranges:
            xmin, xmax, _ = parsed_ranges[x_label]
            ax.set_xlim(xmin, xmax)
        if y_label is not None and y_label in parsed_ranges:
            ymin, ymax, _ = parsed_ranges[y_label]
            ax.set_ylim(ymin, ymax)

    def _regraph_time_domain_from_right_panel(self, parsed_ranges):
        """Redraw Time-Domain Traces using right-panel wavelength and signal ranges."""
        if not hasattr(self, 'data') or self.data is None:
            messagebox.showwarning("No Data", "Load a data file first.")
            return

        if "Wavelength" in parsed_ranges:
            wl_min, wl_max, _ = parsed_ranges["Wavelength"]
            # Reuse the existing tab-local range variable so overlays and the Plot button stay in sync.
            if hasattr(self, "time_trace_wl_range_var"):
                self.time_trace_wl_range_var.set(f"{wl_min:g}, {wl_max:g}")

        self.plot_time_slices()
        self._set_axis_limits_from_right_ranges(
            self.ax_time, parsed_ranges, x_label="Wavelength", y_label="Signal"
        )
        if "Signal" in parsed_ranges:
            signal_min, signal_max, _ = parsed_ranges["Signal"]
            self.time_trace_manual_signal_range = (signal_min, signal_max)
        self.fig_time.tight_layout()
        self.canvas_time.draw()
        self.status_label.config(text="Status: Time-domain traces re-graphed with right-panel ranges.")

    def _regraph_wavelength_domain_from_right_panel(self, parsed_ranges):
        """Redraw Wavelength-Domain Traces using right-panel time and signal ranges."""
        if not hasattr(self, 'data') or self.data is None:
            messagebox.showwarning("No Data", "Load a data file first.")
            return

        self.plot_wavelength_slices()
        self._set_axis_limits_from_right_ranges(
            self.ax_wl, parsed_ranges, x_label="Time", y_label="Signal"
        )
        self.fig_wl.tight_layout()
        self.canvas_wl.draw()
        self.status_label.config(text="Status: Wavelength-domain traces re-graphed with right-panel ranges.")

    def _regraph_2d_map_from_right_panel(self, parsed_ranges):
        """Redraw the 2D map using right-panel wavelength, time, and color ranges."""
        if not hasattr(self, 'data') or self.data is None:
            messagebox.showwarning("No Data", "Load a data file first.")
            return

        if "Color" in parsed_ranges:
            color_min, color_max, _ = parsed_ranges["Color"]
            self.map_vmin_var.set(f"{color_min:.4g}")
            self.map_vmax_var.set(f"{color_max:.4g}")

        if "Time" in parsed_ranges and parsed_ranges["Time"][0] <= 0:
            messagebox.showerror(
                "Input Error",
                "The 2D map uses a logarithmic time axis, so its minimum time must be greater than zero."
            )
            return

        self.update_2d_map()
        self._set_axis_limits_from_right_ranges(
            self.map_ax, parsed_ranges, x_label="Wavelength", y_label="Time"
        )
        self.map_fig.tight_layout()
        self.canvas_tab2.draw()
        self.status_label.config(text="Status: 2D map re-graphed with right-panel ranges.")

    def _safe_axis_min(self, axis, fallback):
        """Return a finite axis minimum or a fallback placeholder."""
        if axis is None:
            return fallback
        value = np.nanmin(axis)
        return float(value) if np.isfinite(value) else fallback

    def _safe_positive_axis_min(self, axis, fallback):
        """Return the smallest finite positive axis value or a fallback."""
        if axis is None:
            return fallback
        values = np.asarray(axis, dtype=float)
        positive_values = values[np.isfinite(values) & (values > 0)]
        return float(np.min(positive_values)) if positive_values.size else fallback

    def _safe_axis_max(self, axis, fallback):
        """Return a finite axis maximum or a fallback placeholder."""
        if axis is None:
            return fallback
        value = np.nanmax(axis)
        return float(value) if np.isfinite(value) else fallback

    def _safe_data_min(self, fallback):
        """Return a finite data minimum or a fallback placeholder."""
        data = getattr(self, "data", None)
        if data is None:
            return fallback
        value = np.nanmin(data)
        return float(value) if np.isfinite(value) else fallback

    def _safe_data_max(self, fallback):
        """Return a finite data maximum or a fallback placeholder."""
        data = getattr(self, "data", None)
        if data is None:
            return fallback
        value = np.nanmax(data)
        return float(value) if np.isfinite(value) else fallback


    def _create_status_bar(self):
        """Creates the bottom status bar"""
        status_bar_frame = ttk.Frame(self, relief='sunken', padding="2 5")
        status_bar_frame.grid(row=1, column=0, columnspan=3, sticky='ew')
        
        self.status_label = ttk.Label(
            status_bar_frame, 
            text="Status: Ready. Load data to begin."
        )
        self.status_label.pack(side='left')
        
        self.progress_bar = ttk.Progressbar(
            status_bar_frame, 
            orient='horizontal', 
            length=200, 
            mode='indeterminate'
        )
        # self.progress_bar.pack(side='right', padx=10)
        # Uncomment above to show. We'll leave it hidden for now.

    def placeholder_command(self):
        """A placeholder function for button clicks."""
        messagebox.showinfo("Placeholder", "This function is not yet implemented.")
        
        # Example of updating status bar
        self.status_label.config(text="Status: Action triggered.")
        
        # Example of using progress bar
        # self.progress_bar.pack(side='right', padx=10)
        # self.progress_bar.start(10)
        # self.after(2000, self.progress_bar.stop)
        # self.after(2000, lambda: self.progress_bar.pack_forget())
        # self.after(2000, lambda: self.status_label.config(text="Status: Ready."))
    def load_data(self):
        file_path = filedialog.askopenfilename(
            title="Select psTA Data File",
            filetypes=[("All Files", "*.*")]
        )

        if not file_path:
            return

        self.status_label.config(text="Status: Running analysis...")
        self.update_idletasks()

        try:
            self.run_matlab_analysis(file_path)
            self.status_label.config(text="Status: Analysis complete.")
        except Exception as e:
            messagebox.showerror("MATLAB Error", str(e))
            self.status_label.config(text="Status: Error during analysis.")

    def run_matlab_analysis(self, data_file):
        ext = Path(data_file).suffix.lower()

        if ext == ".json" or ext == "":
            self._run_original_matlab_pipeline(data_file)

        elif ext == ".txt" or ext == ".csv":
            self._run_new_txt_pipeline(data_file)

        else:
            raise ValueError(f"Unsupported file type: {ext}")
        

    def _sync_working_dat_to_matlab(self):
        """Push the current working MATLAB object into the MATLAB workspace."""
        if self.proc_dat is None:
            self.proc_dat = self.raw_dat
        if self.proc_dat is None:
            raise ValueError("No data loaded.")
        self.eng.workspace['dat'] = self.proc_dat


    def _run_original_matlab_pipeline(self, data_file):
        self.pipeline_type = "original"
        dat = self.eng.TAExperiment(data_file)

        # Keep untouched + working copies
        self.raw_dat = dat
        self.proc_dat = dat
        self.eng.workspace['dat'] = dat

        arr = np.array(self.eng.eval("dat.TAMean"))
        times = np.array(self.eng.eval("dat.times")).flatten()
        wavelengths = np.array(self.eng.eval("dat.wavelengths")).flatten()

        self.raw_data = arr.copy()
        self.proc_data = arr.copy()

        self.data = self.proc_data
        self.times = times
        self.wavelength = wavelengths

        self.map_vmin_var.set(f"{arr.min():.4g}")
        self.map_vmax_var.set(f"{arr.max():.4g}")

        if hasattr(self, "right_dynamic_frame"):
            self._update_right_panel_for_active_tab()

        self.update_2d_map()

    def _infer_energy_axis_from_txt(self, data_file, calibration_coefficients=None):
        """
        Build a Soliton probe axis matching the number of data channels.

        When A, E, and T0 calibration coefficients are supplied, use the
        calibration from the Soliton example: sqrt(E^2 - A/(pixel + T0)).
        Otherwise retain pixel/channel numbers rather than inventing an energy
        scale that could be scientifically misleading.
        """
        # Read only the first scan file to determine shape
        df = pd.read_csv(
            data_file,
            sep=None,          # auto-detect delimiter
            engine="python",
            header=None
        )

        # MATLAB readmatrix will see the first column as time and the rest as data
        n_pixels = df.shape[1] - 1
        if n_pixels <= 0:
            raise ValueError(f"Could not infer a valid energy axis from {data_file}")

        pixels = np.arange(1, n_pixels + 1, dtype=float)
        if calibration_coefficients is None:
            return pixels

        A, E, T0 = calibration_coefficients
        energy_squared = E ** 2 - A / (pixels + T0)
        if np.any(energy_squared <= 0):
            raise ValueError(
                "The Soliton calibration coefficients produce non-positive energy values."
            )
        return np.sqrt(energy_squared)

    def _get_soliton_import_options(self, data_file):
        """Collect the acquisition settings required by SolitonTAExperimentSelfContained."""
        previous = getattr(self, "soliton_import_options", {})
        n_scans = simpledialog.askinteger(
            "Soliton Import",
            "Number of scans to combine:",
            parent=self,
            initialvalue=previous.get("n_scans", 1),
            minvalue=1,
        )
        if n_scans is None:
            return None

        if n_scans > 1 and "scan1" not in Path(data_file).name:
            messagebox.showerror(
                "Soliton Import",
                "For multiple scans, select the scan1 file. The Soliton class derives the remaining scan names from it."
            )
            return None

        is_mirrored = messagebox.askyesno(
            "Soliton Import",
            "Is this a mirrored delay scan?\n\nChoose Yes for a scan containing forward and reverse delay sweeps.",
            parent=self,
        )
        background_upper_bound = simpledialog.askfloat(
            "Soliton Background",
            "Background upper time bound (ps):\nChoose a pre-signal delay, often about -1.5 ps.",
            parent=self,
            initialvalue=previous.get("background_upper_bound", -1.5),
        )
        if background_upper_bound is None:
            return None

        calibration_text = simpledialog.askstring(
            "Soliton Energy Calibration",
            "Optional energy-calibration coefficients A, E, T0 (comma or space separated).\n"
            "Leave blank to use pixel/channel number on the horizontal axis.",
            parent=self,
            initialvalue=previous.get("calibration_text", ""),
        )
        if calibration_text is None:
            return None

        calibration_text = calibration_text.strip()
        calibration = None
        if calibration_text:
            try:
                values = [float(value) for value in calibration_text.replace(",", " ").split()]
                if len(values) != 3:
                    raise ValueError
                calibration = tuple(values)
            except ValueError:
                messagebox.showerror(
                    "Soliton Energy Calibration",
                    "Enter exactly three numeric values: A, E, T0."
                )
                return None

        chirp_text = simpledialog.askstring(
            "Soliton Chirp Correction",
            "Enter five chirp coefficients a4, a3, a2, a1, a0 (comma or space separated).\n"
            "Leave all values as zero for no chirp shift.",
            parent=self,
            initialvalue=previous.get("chirp_text", "0, 0, 0, 0, 0"),
        )
        if chirp_text is None:
            return None

        try:
            chirp_coefficients = [
                float(value) for value in chirp_text.replace(",", " ").split()
            ]
            if len(chirp_coefficients) != 5 or not np.all(np.isfinite(chirp_coefficients)):
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Soliton Chirp Correction",
                "Enter exactly five finite values: a4, a3, a2, a1, a0."
            )
            return None

        options = {
            "n_scans": n_scans,
            "is_mirrored": is_mirrored,
            "background_upper_bound": background_upper_bound,
            "calibration": calibration,
            "calibration_text": calibration_text,
            "chirp_coefficients": chirp_coefficients,
            "chirp_text": chirp_text,
        }
        self.soliton_import_options = options
        return options
    

    def _run_new_txt_pipeline(self, data_file):
        self.pipeline_type = "soliton"
        options = self._get_soliton_import_options(data_file)
        if options is None:
            raise ValueError("Soliton import was cancelled.")

        p = Path(data_file)

        folder_path = p.parent.as_posix() + "/"
        file_name = p.name

        self.eng.addpath(folder_path, nargout=0)

        dat = self.eng.SolitonTAExperimentSelfContained(
            folder_path, file_name, options["n_scans"]
        )

        self.raw_dat = dat
        self.proc_dat = dat
        self.eng.workspace["dat"] = dat

        energy_axis_np = self._infer_energy_axis_from_txt(
            data_file, options["calibration"]
        )
        energy_axis = matlab.double([energy_axis_np.tolist()])  # MATLAB row vector

        is_mirrored = options["is_mirrored"]
        upper_bound = options["background_upper_bound"]
        chirp_coeffs = matlab.double([options["chirp_coefficients"]])

        self.eng.eval("dat = dat;", nargout=0)

        self.eng.mainSoliton_SC(
            dat,
            energy_axis,
            is_mirrored,
            upper_bound,
            "chirpCorrectionCoefficients",
            chirp_coeffs,
            nargout=0
        )

        self.proc_dat = dat
        self.eng.workspace["dat"] = dat

        self._pull_current_dat_arrays()
        self.update_2d_map()

    def _pull_current_dat_arrays(self):
        """
        Pull TA mean, time axis, and x-axis from the current MATLAB dat object.
        Handles both original Soliton.
        """
        if self.proc_dat is None:
            raise ValueError("No MATLAB dat object is loaded.")

        self.eng.workspace["dat"] = self.proc_dat

        if self.pipeline_type == "soliton":
            arr = np.array(self.eng.eval("dat.TAMeanSortedBackgroundSub_CC"))
            times = np.array(self.eng.eval("dat.timesSorted_CC")).flatten()
            xaxis = np.array(self.eng.eval("dat.energyAxis_CC")).flatten()

        else:
            arr = np.array(self.eng.eval("dat.TAMean"))
            times = np.array(self.eng.eval("dat.times")).flatten()
            xaxis = np.array(self.eng.eval("dat.wavelengths")).flatten()

        self.proc_data = arr.copy()
        self.data = self.proc_data
        self.times = times
        self.wavelength = xaxis

        self.map_vmin_var.set(f"{np.nanmin(arr):.4g}")
        self.map_vmax_var.set(f"{np.nanmax(arr):.4g}")

        if hasattr(self, "right_dynamic_frame"):
            self._update_right_panel_for_active_tab()

    def _parse_slice_values(self, text: str, label: str):
        """Parse a comma-separated string into a sorted list of floats.
        Returns None and shows an error dialog on bad input."""
        try:
            vals = [float(v.strip()) for v in text.split(",") if v.strip()]
            if not vals:
                raise ValueError("empty")
            return sorted(vals)
        except ValueError:
            messagebox.showerror(
                "Input Error",
                f"Invalid {label} — enter comma-separated numbers.\n"
                f'Example: "0.1, 1, 10, 100"'
            )
            return None
        
    def _autoscale_time_trace_yaxis(self, wl_range):
        """Autoscale from TA data and visible overlays inside the wavelength window."""
        if self.data is None or self.wavelength is None:
            return

        wl_min, wl_max = wl_range
        mask = (self.wavelength >= wl_min) & (self.wavelength <= wl_max)
        if not np.any(mask):
            return

        visible_values = [self.data[:, mask].ravel()]

        # Loaded spectra should be scaled independently by their selected scale
        # factor, but still participate in the visible graph range when toggled on.
        for overlay_key, overlay in getattr(self, "time_trace_overlays", {}).items():
            if not self._get_overlay_visible_var(overlay_key).get():
                continue
            df = overlay.get("data")
            scale = overlay.get("scale")
            if df is None or scale is None:
                continue
            overlay_wavelengths = df["Wavelength (nm)"].to_numpy(dtype=float)
            overlay_values = float(scale) * df["Absorbance (AU)"].to_numpy(dtype=float)
            overlay_mask = (
                np.isfinite(overlay_wavelengths)
                & np.isfinite(overlay_values)
                & (overlay_wavelengths >= wl_min)
                & (overlay_wavelengths <= wl_max)
            )
            if np.any(overlay_mask):
                visible_values.append(overlay_values[overlay_mask])

        visible_data = np.concatenate(visible_values)
        visible_data = visible_data[np.isfinite(visible_data)]
        if visible_data.size == 0:
            return

        ymin = np.nanmin(visible_data)
        ymax = np.nanmax(visible_data)

        if np.isfinite(ymin) and np.isfinite(ymax):
            if ymin == ymax:
                pad = abs(ymin) * 0.05 if ymin != 0 else 1e-6
            else:
                pad = 0.05 * (ymax - ymin)
            self.ax_time.set_ylim(ymin - pad, ymax + pad)

    def _get_overlay_visible_var(self, overlay_key):
        """Return the Tk BooleanVar associated with a time-trace overlay."""
        if overlay_key == "ground_state":
            return self.ground_state_visible_var
        if overlay_key == "fluorescence":
            return self.fluorescence_visible_var
        raise ValueError(f"Unknown overlay type: {overlay_key}")

    def _get_overlay_toggle(self, overlay_key):
        """Return the Checkbutton widget associated with a time-trace overlay."""
        if overlay_key == "ground_state":
            return self.ground_state_toggle
        if overlay_key == "fluorescence":
            return self.fluorescence_toggle
        raise ValueError(f"Unknown overlay type: {overlay_key}")


    def _remove_time_trace_overlay_artist(self, overlay_key):
        """Remove one overlay line from the time-domain axis and clear its handle."""
        overlay = self.time_trace_overlays.get(overlay_key)
        if overlay is None:
            return

        old_line = overlay.get("line")
        if old_line is not None:
            try:
                old_line.remove()
            # A prior axis clear detaches the artist, and Matplotlib then raises
            # NotImplementedError instead of removing it a second time.
            except (ValueError, NotImplementedError, AttributeError):
                pass
            overlay["line"] = None

    def _validate_time_domain_data(self):
        """Return True when the time-domain plot has consistent arrays to graph."""
        if getattr(self, "data", None) is None or getattr(self, "wavelength", None) is None or getattr(self, "times", None) is None:
            messagebox.showwarning("No Data", "Load a data file first.")
            return False

        if self.data.ndim != 2:
            messagebox.showerror("Data Error", "Time-domain plotting expects a 2D data array.")
            return False

        if len(self.times) != self.data.shape[0]:
            messagebox.showerror(
                "Data Error",
                "The number of time points does not match the number of data rows."
            )
            return False

        if len(self.wavelength) != self.data.shape[1]:
            messagebox.showerror(
                "Data Error",
                "The number of wavelength points does not match the number of data columns."
            )
            return False

        if not np.any(np.isfinite(self.times)) or not np.any(np.isfinite(self.wavelength)):
            messagebox.showerror("Data Error", "Time or wavelength axes contain no finite values.")
            return False

        return True

    def _current_time_domain_wavelength_range(self):
        """Parse the tab-local wavelength range field, if the user entered one."""
        if getattr(self, "time_trace_wl_range_var", None) is None:
            return None

        range_text = self.time_trace_wl_range_var.get().strip()
        if not range_text:
            return None

        return self._parse_range_values(range_text, "wavelength")

    def _time_domain_wavelength_mask(self, wl_range):
        """Return a valid boolean mask for the selected time-domain wavelength range."""
        finite_mask = np.isfinite(self.wavelength)

        if wl_range is None:
            mask = finite_mask
        else:
            wl_min, wl_max = wl_range
            mask = finite_mask & (self.wavelength >= wl_min) & (self.wavelength <= wl_max)

        if not np.any(mask):
            if wl_range is None:
                messagebox.showerror("Range Error", "No finite wavelengths are available to plot.")
            else:
                messagebox.showerror(
                    "Range Error",
                    f"No wavelength channels found in {wl_range[0]:g}-{wl_range[1]:g}."
                )
            return None

        return mask

    def _load_absorbance_csv(self, file_path):
        """
        Load a spectrum CSV with required columns:
            Wavelength (nm), Absorbance (AU)
        and optional column:
            Std.Dev.
        Returns a clean DataFrame sorted by wavelength.
        """
        # 1. Handle encodings safely
        df = None
        for encoding in ["utf-8", "utf-16", "utf-8-sig", "latin1"]:
            try:
                df = pd.read_csv(file_path, encoding=encoding) 
                break
            except UnicodeDecodeError:
                continue
        if df is None:
            messagebox.showerror("CSV Error", "Could not decode the CSV file.")
            return None

        # Clean whitespace from column headers
        df.columns = [str(c).strip() for c in df.columns]

        # 2. Extract Data Series based on headers or positions
        if "Wavelength (nm)" in df.columns and "Absorbance (AU)" in df.columns:
            wavelength_raw = df["Wavelength (nm)"]
            absorbance_raw = df["Absorbance (AU)"]
            std_dev_raw = df.get("Std.Dev.", None)  # Safe fetch if it exists
        else:
            # Fallback to positional parsing if headers don't match
            if df.shape[1] < 2:
                messagebox.showerror(
                    "CSV Error",
                    "CSV must contain at least two columns: wavelength and absorbance."
                )
                return None
            wavelength_raw = df.iloc[:, 0]
            absorbance_raw = df.iloc[:, 1]
            # Assume 3rd column is Std.Dev if it exists
            std_dev_raw = df.iloc[:, 2] if df.shape[1] >= 3 else None

        # 3. Force values to numeric data types (coerce strings/errors to NaN)
        wavelength_numeric = pd.to_numeric(wavelength_raw, errors="coerce")
        absorbance_numeric = pd.to_numeric(absorbance_raw, errors="coerce")

        # 4. Filter for rows where BOTH main metrics are valid numbers
        valid_mask = wavelength_numeric.notna() & absorbance_numeric.notna()
        
        if not valid_mask.any():
            messagebox.showerror(
                "CSV Error",
                "No valid numeric wavelength/absorbance rows found."
            )
            return None

        # 5. Build a fresh, standardized DataFrame using the valid rows
        clean_data = {
            "Wavelength (nm)": wavelength_numeric[valid_mask],
            "Absorbance (AU)": absorbance_numeric[valid_mask]
        }

        # Handle standard deviation if present
        if std_dev_raw is not None:
            std_dev_numeric = pd.to_numeric(std_dev_raw, errors="coerce")
            clean_data["Std.Dev."] = std_dev_numeric[valid_mask]

        clean_df = pd.DataFrame(clean_data)

        # 6. Safely sort by Wavelength
        return clean_df.sort_values(by="Wavelength (nm)").reset_index(drop=True)

    def _get_time_trace_wavelength_mask(self, overlay_wavelengths):
        """Apply the current time-domain wavelength range to an overlay spectrum."""
        range_text = getattr(self, "time_trace_wl_range_var", tk.StringVar(value="")).get().strip()
        if not range_text:
            return np.ones_like(overlay_wavelengths, dtype=bool), None

        wl_range = self._parse_range_values(range_text, "wavelength")
        if wl_range is None:
            return None, None

        wl_min, wl_max = wl_range
        mask = (overlay_wavelengths >= wl_min) & (overlay_wavelengths <= wl_max)
        return mask, wl_range

    def plot_scaled_absorbance_overlay(self, overlay_key):
        """
        Select a CSV, ask for a scale factor, then overlay
        scale factor * Absorbance (AU) on the Time-Domain Traces plot.
        """
        if overlay_key not in getattr(self, "time_trace_overlays", {}):
            messagebox.showerror("Overlay Error", f"Unknown overlay type: {overlay_key}")
            return

        overlay = self.time_trace_overlays[overlay_key]
        label = overlay["label"]

        file_path = filedialog.askopenfilename(
            title=f"Select {label} CSV",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if not file_path:
            return

        scale = simpledialog.askfloat(
            f"{label} Scaling Factor",
            f"Enter scaling factor for {label}:\n\nPlotted value = scaling factor * Absorbance (AU)",
            parent=self,
            initialvalue=1.0
        )
        if scale is None:
            return

        try:
            df = self._load_absorbance_csv(file_path)
        except Exception as e:
            messagebox.showerror("CSV Error", str(e))
            self.status_label.config(text=f"Status: Could not load {label} CSV.")
            return

        # Do not enable toggles or overwrite a good overlay with a failed import.
        if df is None or df.empty:
            self.status_label.config(text=f"Status: Could not load {label} CSV.")
            return

        overlay.update({
            "data": df,
            "scale": float(scale),
            "file_path": file_path,
        })

        # Set visibility before plotting so re-importing an overlay after it was hidden
        # does not create a hidden line while the checkbox says it is shown.
        self._get_overlay_visible_var(overlay_key).set(True)
        self._get_overlay_toggle(overlay_key).config(state="normal")

        plotted = self._plot_or_refresh_time_trace_overlay(overlay_key, redraw=False, warn_if_empty=True)
        if (
            plotted
            and getattr(self, "data", None) is not None
            and getattr(self, "wavelength", None) is not None
            and getattr(self, "time_trace_manual_signal_range", None) is None
        ):
            self._autoscale_time_trace_yaxis(self._time_trace_visible_wavelength_range())
        self._refresh_time_trace_legend()
        self.fig_time.tight_layout()
        self.canvas_time.draw()

        if plotted:
            self.status_label.config(
                text=f"Status: {label} plotted with scale factor {float(scale):g}."
            )
        else:
            self.status_label.config(
                text=f"Status: {label} loaded, but no points are inside the current wavelength range."
            )

    def _plot_or_refresh_time_trace_overlay(self, overlay_key, redraw=True, warn_if_empty=False):
        """Draw or redraw one stored absorbance overlay on the time-domain axis."""
        overlay = self.time_trace_overlays.get(overlay_key)
        if overlay is None:
            return False

        df = overlay.get("data")
        scale = overlay.get("scale")
        if df is None or scale is None:
            return False

        self._remove_time_trace_overlay_artist(overlay_key)

        required_cols = {"Wavelength (nm)", "Absorbance (AU)"}
        if not required_cols.issubset(set(df.columns)):
            messagebox.showerror(
                "CSV Error",
                f"{overlay['label']} data is missing required columns: Wavelength (nm), Absorbance (AU)."
            )
            return False

        wavelengths = df["Wavelength (nm)"].to_numpy(dtype=float)
        scaled_absorbance = float(scale) * df["Absorbance (AU)"].to_numpy(dtype=float)
        valid = np.isfinite(wavelengths) & np.isfinite(scaled_absorbance)
        wavelengths = wavelengths[valid]
        scaled_absorbance = scaled_absorbance[valid]

        if wavelengths.size == 0:
            if warn_if_empty:
                messagebox.showwarning(
                    "CSV Warning",
                    f"{overlay['label']} has no valid numeric wavelength/absorbance points."
                )
            return False

        mask, wl_range = self._get_time_trace_wavelength_mask(wavelengths)
        if mask is None:
            return False
        if not np.any(mask):
            if warn_if_empty:
                messagebox.showwarning(
                    "Range Warning",
                    f"{overlay['label']} has no wavelengths inside the selected visible range."
                )
            return False

        (line,) = self.ax_time.plot(
            wavelengths[mask],
            scaled_absorbance[mask],
            lw=2.0,
            ls="--",
            label=f"{overlay['label']} x {float(scale):g}"
        )
        overlay["line"] = line

        visible = self._get_overlay_visible_var(overlay_key).get()
        line.set_visible(visible)

        if redraw:
            self._refresh_time_trace_legend()
            self.fig_time.tight_layout()
            self.canvas_time.draw()

        return True

    def _refresh_time_trace_overlays(self):
        """Redraw every loaded overlay after the time-domain axis is cleared or replotted."""
        for overlay_key in self.time_trace_overlays:
            self._plot_or_refresh_time_trace_overlay(overlay_key, redraw=False, warn_if_empty=False)

    def toggle_time_trace_overlay(self, overlay_key):
        """Show/hide a plotted absorbance overlay without deleting it."""
        overlay = self.time_trace_overlays.get(overlay_key)
        if overlay is None:
            return

        line = overlay.get("line")

        # After Clear, Plot, or axis re-creation, the stored data can exist while the
        # Matplotlib line handle is gone. Recreate it so the checkbox always works.
        if line is None and overlay.get("data") is not None and overlay.get("scale") is not None:
            self._plot_or_refresh_time_trace_overlay(overlay_key, redraw=False, warn_if_empty=False)
            line = overlay.get("line")

        if line is None:
            return

        line.set_visible(self._get_overlay_visible_var(overlay_key).get())
        if (
            getattr(self, "data", None) is not None
            and getattr(self, "wavelength", None) is not None
            and getattr(self, "time_trace_manual_signal_range", None) is None
        ):
            self._autoscale_time_trace_yaxis(self._time_trace_visible_wavelength_range())
        self._refresh_time_trace_legend()
        self.fig_time.tight_layout()
        self.canvas_time.draw()

    def _refresh_time_trace_legend(self):
        """Show only visible time traces and overlays in the legend."""
        handles, labels = self.ax_time.get_legend_handles_labels()
        visible_pairs = [
            (handle, label)
            for handle, label in zip(handles, labels)
            if getattr(handle, "get_visible", lambda: True)()
        ]

        legend = self.ax_time.get_legend()
        if legend is not None:
            legend.remove()

        if visible_pairs:
            visible_handles, visible_labels = zip(*visible_pairs)
            self.ax_time.legend(visible_handles, visible_labels, fontsize=8, framealpha=0.7)

    def _time_trace_visible_wavelength_range(self):
        """Return the current time-trace wavelength window, or the full data range."""
        wl_range = self._current_time_domain_wavelength_range()
        if wl_range is not None:
            return wl_range
        return float(np.nanmin(self.wavelength)), float(np.nanmax(self.wavelength))

    def plot_time_slices(self):
        """Plot ΔA vs wavelength for each requested time value (Tab 1)."""
        if not self._validate_time_domain_data():
            return

        # Pressing Plot or applying a wavelength window returns to automatic
        # signal scaling. Re-graph applies an explicit Signal range afterwards.
        self.time_trace_manual_signal_range = None

        targets = self._parse_slice_values(self.time_slice_var.get(), "time")
        if targets is None:
            return

        wl_range = self._current_time_domain_wavelength_range()
        if wl_range is None and getattr(self, "time_trace_wl_range_var", None) is not None:
            # _current_time_domain_wavelength_range already displays the parse error.
            range_text = self.time_trace_wl_range_var.get().strip()
            if range_text:
                return

        plot_mask = self._time_domain_wavelength_mask(wl_range)
        if plot_mask is None:
            return

        self.ax_time.cla()
        self.ax_time.set_title("Time-Domain Traces")
        self.ax_time.set_xlabel("Wavelength (nm)")
        self.ax_time.set_ylabel("ΔA")
        self.ax_time.grid(True)
        self.ax_time.axhline(0, color='gray', lw=0.8, ls='--')

        cmap = plt.get_cmap("plasma")
        colors = [cmap(i / max(len(targets) - 1, 1)) for i in range(len(targets))]

        for t_req, color in zip(targets, colors):
            idx = int(np.nanargmin(np.abs(self.times - t_req)))
            t_actual = self.times[idx]
            self.ax_time.plot(
                self.wavelength[plot_mask],
                self.data[idx, :][plot_mask],
                color=color,
                lw=1.5,
                label=f"{t_actual:.3g} ps"
            )

        if wl_range is not None:
            self.ax_time.set_xlim(*wl_range)
            self._autoscale_time_trace_yaxis(wl_range)
        else:
            full_range = (float(np.nanmin(self.wavelength)), float(np.nanmax(self.wavelength)))
            self.ax_time.set_xlim(*full_range)
            self._autoscale_time_trace_yaxis(full_range)

        self._refresh_time_trace_overlays()
        self._refresh_time_trace_legend()
        self.fig_time.tight_layout()
        self.canvas_time.draw()
        self.status_label.config(text="Status: Time-domain traces plotted.")

    def plot_wavelength_slices(self):
        """Plot ΔA vs time for each requested wavelength value.

        If an averaging range is provided, each kinetic trace is averaged over
        wl_req +/- range_nm.
        Example:
            wavelengths = 500, 550, 600
            range = 10
            averages over 490-510, 540-560, 590-610 nm
        """
        if not hasattr(self, 'data') or self.data is None:
            messagebox.showwarning("No Data", "Load a data file first.")
            return

        targets = self._parse_slice_values(self.wl_slice_var.get(), "wavelength")
        if targets is None:
            return

        try:
            avg_range = float(self.wl_avg_range_var.get().strip())
            if avg_range < 0:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Input Error",
                'Invalid wavelength averaging range — enter a non-negative number.\n'
                'Example: "10" averages each trace over wavelength +/- 10 nm.'
            )
            return

        self.ax_wl.cla()
        self.ax_wl.set_title("Wavelength-Domain Traces")
        self.ax_wl.set_xlabel("Time (ps)")
        self.ax_wl.set_ylabel("ΔA")
        self.ax_wl.grid(True)
        self.ax_wl.axhline(0, color='gray', lw=0.8, ls='--')

        cmap = plt.get_cmap("viridis")
        colors = [cmap(i / max(len(targets) - 1, 1)) for i in range(len(targets))]

        for wl_req, color in zip(targets, colors):
            if avg_range == 0:
                # Original behavior: use nearest wavelength channel
                idx = int(np.argmin(np.abs(self.wavelength - wl_req)))
                wl_actual = self.wavelength[idx]
                trace = self.data[:, idx]
                label = f"{wl_actual:.4g} nm"

            else:
                wl_min = wl_req - avg_range
                wl_max = wl_req + avg_range
                mask = (self.wavelength >= wl_min) & (self.wavelength <= wl_max)

                if not np.any(mask):
                    messagebox.showwarning(
                        "Range Warning",
                        f"No wavelength channels found in {wl_min:g}-{wl_max:g} nm."
                    )
                    continue

                # Average ΔA over the selected wavelength columns at each time point
                trace = np.nanmean(self.data[:, mask], axis=1)

                actual_min = np.nanmin(self.wavelength[mask])
                actual_max = np.nanmax(self.wavelength[mask])
                label = f"{wl_req:.4g} nm avg ({actual_min:.4g}-{actual_max:.4g})"

            self.ax_wl.plot(
                self.times,
                trace,
                color=color,
                lw=1.5,
                label=label
            )

        self.ax_wl.legend(fontsize=8, framealpha=0.7)
        self.fig_wl.tight_layout()
        self.canvas_wl.draw()

        if avg_range == 0:
            self.status_label.config(text="Status: Wavelength-domain traces plotted.")
        else:
            self.status_label.config(
                text=f"Status: Wavelength-domain traces averaged over +/- {avg_range:g} nm."
            )

    def _parse_range_values(self, text: str, label: str):
        """Parse a comma-separated pair of floats into (low, high)."""
        try:
            vals = [float(v.strip()) for v in text.split(",") if v.strip()]
            if len(vals) != 2:
                raise ValueError("Expected exactly two numbers.")
            low, high = vals
            if low > high:
                low, high = high, low
            return low, high
        except ValueError:
            messagebox.showerror(
                "Input Error",
                f"Invalid {label} range — enter two comma-separated numbers.\n"
                f'Example: "450, 700"'
            )
            return None


    def apply_time_trace_wavelength_range(self):
        """Limit the visible wavelength range in the time-domain traces plot only."""
        if not self._validate_time_domain_data():
            return

        parsed = self._parse_range_values(self.time_trace_wl_range_var.get(), "wavelength")
        if parsed is None:
            return

        wl_min, wl_max = parsed
        mask = self._time_domain_wavelength_mask((wl_min, wl_max))
        if mask is None:
            return

        # Replot instead of only changing axis limits, so time traces, overlays, and
        # legends are all regenerated from the same range state.
        self.plot_time_slices()
        self.status_label.config(
            text=f"Status: Time-trace wavelength view set to {wl_min:g} - {wl_max:g} nm."
        )

    def reset_time_trace_wavelength_range(self):
        """Reset the time-domain traces x-axis to the full wavelength span."""
        if self.wavelength is None:
            return
        self.time_trace_wl_range_var.set(f"{np.nanmin(self.wavelength):.4g}, {np.nanmax(self.wavelength):.4g}")

        if self._validate_time_domain_data():
            self.plot_time_slices()
        else:
            self.ax_time.set_xlim(np.nanmin(self.wavelength), np.nanmax(self.wavelength))
            self.fig_time.tight_layout()
            self.canvas_time.draw()

        self.status_label.config(text="Status: Time-trace wavelength view reset.")

    def _clear_slice_axis(self, ax, canvas, xlabel, ylabel, title):
        """Reset a slice plot to its blank placeholder state."""
        ax.cla()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True)
        ax.text(
            0.5, 0.5, "Load data, then enter values above",
            ha='center', va='center', transform=ax.transAxes,
            fontsize=12, color='gray', alpha=0.5
        )

        if ax is getattr(self, "ax_time", None):
            for overlay in getattr(self, "time_trace_overlays", {}).values():
                overlay["line"] = None

        canvas.draw()


    def reset_map_clamp(self):
        """Reset min/max fields to the actual data range and redraw."""
        if self.data is not None:
            self.map_vmin_var.set(f"{self.data.min():.4g}")
            self.map_vmax_var.set(f"{self.data.max():.4g}")
            self.update_2d_map()

    def update_2d_map(self):
        time_min = self._safe_positive_axis_min(getattr(self, "times", None), None)
        if time_min is None:
            messagebox.showwarning(
                "Invalid Time Axis",
                "The 2D map requires at least one finite time value greater than zero for its logarithmic scale."
            )
            return

        finite_times = np.asarray(self.times, dtype=float)
        finite_times = finite_times[np.isfinite(finite_times)]
        time_max = float(np.max(finite_times))
        if time_max <= time_min:
            time_max = time_min * 10

        self.map_fig.clear()
        
        # Re-add the main axis
        self.map_ax = self.map_fig.add_subplot(111)
        # Parse vmin/vmax from entry fields; fall back to data range if invalid
        try:
            vmin = float(self.map_vmin_var.get())
        except (ValueError, AttributeError):
            vmin = self.data.min() if self.data is not None else None

        try:
            vmax = float(self.map_vmax_var.get())
        except (ValueError, AttributeError):
            vmax = self.data.max() if self.data is not None else None

        if not np.isfinite(vmin) or not np.isfinite(vmax):
            messagebox.showerror("Map Color Range", "Map color limits must be finite numbers.")
            return
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        if vmin == vmax:
            color_padding = abs(vmin) * 0.05 if vmin != 0 else 1e-12
            vmin -= color_padding
            vmax += color_padding

        # Keep zero at the neutral midpoint of Matplotlib's RdBu map, even
        # when the user selects a one-sided clamp range.
        norm_vmin = vmin if vmin < 0 else -max(abs(vmax), 1e-12)
        norm_vmax = vmax if vmax > 0 else max(abs(vmin), 1e-12)
        color_norm = TwoSlopeNorm(vmin=norm_vmin, vcenter=0.0, vmax=norm_vmax)

        # Clamp the data so values outside [vmin, vmax] are pinned to the boundary
        display_data = np.clip(self.data, vmin, vmax)
        
        # Plot data
        im = self.map_ax.imshow(
            display_data,
            aspect='auto',
            origin='lower',
            extent=[
                self.wavelength.min(),
                self.wavelength.max(),
                time_min,
                time_max
            ],
            cmap="RdBu",
            norm=color_norm,
        )

        self.map_ax.set_title("2D Map View (TA Mean)")
        self.map_ax.set_xlabel("Wavelength (nm)")
        self.map_ax.set_ylabel("Time (ps, log scale)")
        self.map_ax.set_yscale("log")
        self.map_ax.set_ylim(time_min, time_max)

        # Add colorbar if toggled
        if self.colorbar_toggle_var.get():
            self.map_colorbar = self.map_fig.colorbar(im, ax=self.map_ax)
        else:
            self.map_colorbar = None

        # Redraw the canvas
        self.canvas_tab2.draw()
            
    def background_correction(self):
        """
        Apply background correction to the current working copy.
        """
        if self.raw_dat is None:
            messagebox.showwarning("No Data", "Load data before running background correction.")
            return

        upper_val = simpledialog.askfloat(
            "Background Correction",
            "Enter upper bound for background subtraction:",
            parent=self
        )
        if upper_val is None:
            return

        self.status_label.config(
            text=f"Status: Running background correction (upper={upper_val})..."
        )
        self.update_idletasks()

        try:
            # Start from the current working copy
            self._sync_working_dat_to_matlab()

            # IMPORTANT: assign the result back to dat in MATLAB
            if self.pipeline_type == "soliton":
                self.eng.eval(
                    f"backgroundSubtractSoliton_SC(dat, {float(upper_val)});",
                    nargout=0
                )

                # Keep the current chirp-corrected field in sync after new background subtraction.
                self.eng.eval(
                    "chirpCorrectionSoliton_SC(dat, dat.chirpCorrectionCoefficients);",
                    nargout=0
                )

            else:
                self.eng.eval(
                    f"subtractTABackground(dat, {float(upper_val)});",
                    nargout=0
                )

            self.proc_dat = self.eng.workspace["dat"]

            self._pull_current_dat_arrays()

            self.update_2d_map()

            self.status_label.config(text="Status: Background correction complete.")
            messagebox.showinfo("Success", "Background subtraction finished successfully.")

        except Exception as e:
            messagebox.showerror("MATLAB Error", str(e))
            self.status_label.config(text=f"Status: Error during background correction. {e}")

    def dispersion_correction(self):
        """Show options for dispersion correction."""
        if getattr(self, "eng", None) is None:
            return messagebox.showerror("MATLAB Not Running", "MATLAB engine not available.")
        if self.data is None:
            return messagebox.showwarning("No Data", "Load data before running dispersion correction.")

        # Small top-level chooser window
        chooser = tk.Toplevel(self)
        chooser.title("Dispersion Correction Options")
        chooser.transient(self)
        chooser.grab_set()
        chooser.geometry("450x200")

        ttk.Label(chooser, text="Choose dispersion correction action:", padding=8).pack()

        if self.pipeline_type == "soliton":
            ttk.Label(
                chooser,
                text="Soliton data uses its five-coefficient chirp polynomial.",
                wraplength=400,
            ).pack(padx=20, pady=(0, 4))
            ttk.Button(
                chooser,
                text="Apply Soliton Chirp Correction",
                command=lambda: [chooser.destroy(), self._open_soliton_chirp_dialog()],
            ).pack(fill="x", padx=20, pady=4)
            chooser.bind("<Escape>", lambda e: chooser.destroy())
            return

        # Using lambda allows us to destroy the window and call the function in one line
        ttk.Button(chooser, text="Correct Dispersions (MATLAB popup)", 
                   command=lambda: [chooser.destroy(), self._call_matlab_correctdispersion()]).pack(fill='x', padx=20, pady=4)
        ttk.Button(chooser, text="Apply External Dispersion Correction", 
                   command=lambda: [chooser.destroy(), self._open_apply_external_dialog()]).pack(fill='x', padx=20, pady=4)

        chooser.bind("<Escape>", lambda e: chooser.destroy())

    def _call_matlab_correctdispersion(self):
        """Call dispersion correction on the current working copy."""
        if self.pipeline_type == "soliton":
            self._open_soliton_chirp_dialog()
            return

        self._update_status("Launching MATLAB dispersion-correction popup...")
        try:
            self._sync_working_dat_to_matlab()

            cmd = "dat.correctDispersion();"
            self.eng.eval(cmd, nargout=0)

            self.proc_dat = self.eng.workspace['dat']
            self._refresh_dat_from_matlab("Dispersion correction (MATLAB) finished.")
        except Exception as e:
            messagebox.showerror("MATLAB Error", f"Error while calling correctdispersion:\n{e}")
            self._update_status("Error calling correctdispersion.")

    def _open_apply_external_dialog(self):
        """Prompt user for coefficients using a custom Toplevel and apply them in MATLAB."""
        if self.pipeline_type == "soliton":
            self._open_soliton_chirp_dialog()
            return

        default_vec = [-7.44208608527694e-06, 0.00672243358275736, -0.951026821920219]

        dlg = tk.Toplevel(self)
        dlg.title("Apply External Dispersion Correction")
        dlg.transient(self)
        dlg.grab_set()
        dlg.resizable(False, False)
        dlg.geometry("520x160")

        ttk.Label(dlg, text="Enter three dispersion coefficients (comma or space separated):", 
                  wraplength=480, padding=8).pack(anchor="w")

        entry_var = tk.StringVar(value=", ".join(f"{v:.18g}" for v in default_vec))
        ttk.Entry(dlg, textvariable=entry_var, width=60).pack(padx=8, pady=6, fill='x')

        status_lbl = ttk.Label(dlg, text="", foreground="red")
        status_lbl.pack(anchor="w", padx=8)

        btn_frame = ttk.Frame(dlg)
        btn_frame.pack(fill='x', pady=10, padx=8)

        def on_ok(event=None):
            txt = entry_var.get().strip()
            try:
                parts = [p for p in txt.replace(",", " ").split() if p]
                if len(parts) != 3:
                    raise ValueError("Please enter exactly 3 numbers.")
                coeffs = [float(p) for p in parts]
            except Exception as ex:
                status_lbl.config(text=f"Invalid input: {ex}")
                return

            dlg.destroy()

            self._update_status("Applying external dispersion correction...")
            try:
                self._sync_working_dat_to_matlab()

                vec_str = f"[{' '.join(map(str, coeffs))}]"
                self.eng.eval(
                    f"applyExternalDispersionCorrection(dat, {vec_str});",
                    nargout=0
                )

                self.proc_dat = self.eng.workspace['dat']
                self._refresh_dat_from_matlab("External dispersion correction applied.")
            except Exception as e:
                messagebox.showerror("MATLAB Error", f"Error applying external correction:\n{e}")
                self._update_status("Error during external dispersion correction.")
        def on_cancel(event=None):
            dlg.destroy()

        ttk.Button(btn_frame, text="Cancel", command=on_cancel).pack(side='right', padx=6)
        ttk.Button(btn_frame, text="Apply", command=on_ok).pack(side='right')

        dlg.bind("<Return>", on_ok)
        dlg.bind("<Escape>", on_cancel)

    def _open_soliton_chirp_dialog(self):
        """Apply the five-coefficient Soliton chirp polynomial without MATLAB popups."""
        default_coeffs = [0.0] * 5
        try:
            self._sync_working_dat_to_matlab()
            loaded_coeffs = np.asarray(
                self.eng.eval("dat.chirpCorrectionCoefficients"), dtype=float
            ).reshape(-1)
            if loaded_coeffs.size == 5 and np.all(np.isfinite(loaded_coeffs)):
                default_coeffs = loaded_coeffs.tolist()
        except Exception:
            # A newly loaded object may not yet expose coefficients; zeros are the
            # documented no-correction default for this class.
            pass

        dlg = tk.Toplevel(self)
        dlg.title("Soliton Chirp Correction")
        dlg.transient(self)
        dlg.grab_set()
        dlg.resizable(False, False)
        dlg.geometry("650x190")

        ttk.Label(
            dlg,
            text=(
                "Enter five chirp coefficients a4, a3, a2, a1, a0.\n"
                "The Soliton correction is a4·pixel⁴ + a3·pixel³ + a2·pixel² + a1·pixel + a0 (ps).\n"
                "All zeros keeps the unshifted, background-subtracted data."
            ),
            wraplength=620,
            padding=8,
        ).pack(anchor="w")

        entry_var = tk.StringVar(value=", ".join(f"{value:.12g}" for value in default_coeffs))
        ttk.Entry(dlg, textvariable=entry_var, width=78).pack(padx=8, pady=6, fill="x")
        status_lbl = ttk.Label(dlg, text="", foreground="red")
        status_lbl.pack(anchor="w", padx=8)

        btn_frame = ttk.Frame(dlg)
        btn_frame.pack(fill="x", pady=10, padx=8)

        def on_apply(event=None):
            try:
                coeffs = [float(value) for value in entry_var.get().replace(",", " ").split()]
                if len(coeffs) != 5 or not np.all(np.isfinite(coeffs)):
                    raise ValueError("Enter exactly five finite numbers.")
            except ValueError as exc:
                status_lbl.config(text=f"Invalid coefficients: {exc}")
                return

            dlg.destroy()
            self._update_status("Applying Soliton chirp correction...")
            try:
                self._sync_working_dat_to_matlab()
                self.eng.workspace["soliton_chirp_coeffs"] = matlab.double([coeffs])
                self.eng.eval(
                    "chirpCorrectionSoliton_SC(dat, soliton_chirp_coeffs);",
                    nargout=0,
                )
                self.proc_dat = self.eng.workspace["dat"]
                self._refresh_dat_from_matlab("Soliton chirp correction applied.")
            except Exception as exc:
                messagebox.showerror("MATLAB Error", f"Could not apply Soliton chirp correction:\n{exc}")
                self._update_status("Error during Soliton chirp correction.")

        ttk.Button(btn_frame, text="Cancel", command=dlg.destroy).pack(side="right", padx=6)
        ttk.Button(btn_frame, text="Apply", command=on_apply).pack(side="right")
        dlg.bind("<Return>", on_apply)
        dlg.bind("<Escape>", lambda event: dlg.destroy())

    def _update_status(self, msg):
        """Helper to update status label and force UI refresh."""
        self.status_label.config(text=f"Status: {msg}")
        self.update_idletasks()

    def _refresh_dat_from_matlab(self, success_msg):
        """
        Refresh Python-side arrays from the current MATLAB dat object.
        Works for both old TAExperiment and new SolitonTAExperimentSelfContained.
        """
        try:
            self.proc_dat = self.eng.workspace["dat"]

            self._pull_current_dat_arrays()

            if hasattr(self, "update_2d_map"):
                self.update_2d_map()

            self._update_status(success_msg)
            messagebox.showinfo("Done", success_msg)

        except Exception as e:
            messagebox.showerror(
                "MATLAB Error",
                f"Could not refresh dat from MATLAB workspace:\n{e}"
            )
            self._update_status("Error refreshing MATLAB data.")
            
    def on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            try:
                if hasattr(self, "eng"):
                    self.eng.quit()
            except:
                pass
            self.destroy()


if __name__ == "__main__":
    # Set high-DPI awareness for Windows
    try:
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    except:
        pass # Fails on non-Windows systems
    app = PsTAAnalysisApp()
    app.mainloop()
