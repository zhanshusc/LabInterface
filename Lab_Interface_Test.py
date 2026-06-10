import tkinter as tk
from tkinter import ttk
from tkinter import messagebox  
from tkinter import filedialog
from tkinter import simpledialog
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import matlab.engine
import threading
from pathlib import Path

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
class PsTAAnalysisApp(tk.Tk):
    """
    A Python tkinter application replicating the psTA Analysis Suite layout
    described for MATLAB App Designer.
    """ 
    def __init__(self):
        super().__init__()

        # Main App Window Setup
        self.title("psTA Analysis Suite (Python)")
        self.geometry("1200x700")
        
        # Main Grid Layout (3 Columns)
        # Configure the root window's grid
        # Column 0: Left Panel (weight 1)
        self.grid_columnconfigure(0, weight=1, minsize=250) 
        # Column 1: Middle Panel (weight 3 - gets more space)
        self.grid_columnconfigure(1, weight=3, minsize=500)
        # Column 2: Right Panel (weight 1)
        self.grid_columnconfigure(2, weight=1, minsize=250)
        # Row 0: Main content
        self.grid_rowconfigure(0, weight=1)
        # Row 1: Status Bar
        self.grid_rowconfigure(1, weight=0)

        # Create the Three Major Panels
        # We use ttk.Frame for a modern look
        self.left_panel = ttk.Frame(self, padding="10")
        self.middle_panel = ttk.Frame(self, padding="10")
        self.right_panel = ttk.Frame(self, padding="10")

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
            data_io_frame, text ="Load TCSPC Data", command=self.load_data).pack(fill='x', pady=2)
        
        
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
        tab_control = ttk.Notebook(frame)
        
        # Tab Frames
        tab1 = ttk.Frame(tab_control, padding="10")
        tab2 = ttk.Frame(tab_control, padding="10")
        tab3 = ttk.Frame(tab_control, padding="10")
        tab4 = ttk.Frame(tab_control, padding="10")
        
        # Add tabs to the notebook
        tab_control.add(tab1, text='Time-Domain Traces')
        tab_control.add(tab4, text='Wavelength-Domain Traces')
        tab_control.add(tab2, text='2D Map View')
        tab_control.add(tab3, text='Fitting')
        
        # Make the notebook fill the middle panel
        tab_control.pack(expand=1, fill='both')

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

        ttk.Label(tab4_controls, text="Wavelengths (nm), comma-separated:").pack(side='left', padx=5)
        self.wl_slice_var = tk.StringVar(value="500, 550, 600, 650")
        ttk.Entry(tab4_controls, textvariable=self.wl_slice_var, width=35).pack(side='left', padx=5)
        ttk.Label(tab4_controls, text="Average range (+/- nm):").pack(side='left', padx=(10, 5))
        self.wl_avg_range_var = tk.StringVar(value="0")
        ttk.Entry(tab4_controls, textvariable=self.wl_avg_range_var, width=8).pack(side='left', padx=5)
        ttk.Button(
            tab4_controls, text="Plot", command=self.plot_wavelength_slices
        ).pack(side='left', padx=2)
        ttk.Button(
            tab4_controls, text="Clear",
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
        
        canvas_tab3_top = FigureCanvasTkAgg(before_fig, master=tab3_main_frame)
        canvas_tab3_top.draw()
        canvas_tab3_top.get_tk_widget().grid(row=0, column=0, sticky="nsew", pady=2)

        # Bottom Plot: TCSPC
        after_fig = Figure(dpi=100)
        after_ax = after_fig.add_subplot(111)
        after_ax.set_title("TCSPC")
        after_ax.set_xlabel("Time (ps/ns)")
        after_ax.set_ylabel("Emission Intensity")
        after_ax.grid(True)
        after_ax.text(0.5, 0.5, "TCSPC", horizontalalignment='center', verticalalignment='center', transform=after_ax.transAxes, fontsize=14, color='gray', alpha=0.5)
        
        canvas_tab3_bot = FigureCanvasTkAgg(after_fig, master=tab3_main_frame)
        canvas_tab3_bot.draw()
        canvas_tab3_bot.get_tk_widget().grid(row=1, column=0, sticky="nsew", pady=2)


    def _create_right_panel(self):
        """Populates the Right Panel (Analysis & Export)"""
        frame = self.right_panel

        # Section: Kinetic Fitting
        fit_frame = ttk.LabelFrame(frame, text="Kinetic Fitting", padding="10")
        fit_frame.pack(fill='x', pady=5)
        
        ttk.Label(fit_frame, text="Fit Mode:").pack()
        self.fit_mode_combo = ttk.Combobox(
            fit_frame, 
            values=["single exp", "bi exp", "tri exp"]
        )
        self.fit_mode_combo.current(0)
        self.fit_mode_combo.pack(fill='x', pady=2)

        # Fit Range
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
        
        # Results Box
        ttk.Label(fit_frame, text="Fit Results:").pack()
        self.results_text = tk.Text(fit_frame, height=8, width=30)
        self.results_text.insert(
            '1.0', "τ₁, τ₂, amplitudes, χ²...\n"
        )
        self.results_text.config(state='disabled', bg='#f0f0f0') # Read-only
        self.results_text.pack(fill='x', expand=True, pady=5)
        
        # Section: Export Options 
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

        self.update_2d_map()

    def _infer_energy_axis_from_txt(self, data_file, e_min=1.5, e_max=3.5):
        """
        Infer the number of probe-energy channels from the TXT file and build
        a matching energy axis.
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

        energy_axis = np.linspace(e_min, e_max, n_pixels)
        return energy_axis
    

    def _run_new_txt_pipeline(self, data_file):
        self.pipeline_type = "soliton"
        p = Path(data_file)

        folder_path = str(p.parent) + "/"
        file_name = p.name

        self.eng.addpath(folder_path, nargout=0)

        dat = self.eng.SolitonTAExperimentSelfContained(folder_path, file_name, 1)

        self.raw_dat = dat
        self.proc_dat = dat
        self.eng.workspace["dat"] = dat

        energy_axis_np = self._infer_energy_axis_from_txt(data_file, e_min=1.5, e_max=3.5)
        energy_axis = matlab.double([energy_axis_np.tolist()])  # MATLAB row vector

        is_mirrored = False
        upper_bound = -1.5

        # IMPORTANT:
        # In the current MATLAB class, all-zero chirp coefficients trigger getpts().
        # That will open an interactive MATLAB point-picking window and can make the
        # Python UI look frozen. Use a non-interactive path or edit MATLAB as below.
        chirp_coeffs = matlab.double([[0, 0, 0, 0, 0]])

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
        """Autoscale the time-trace y-axis using only data inside the wavelength window."""
        if self.data is None or self.wavelength is None:
            return

        wl_min, wl_max = wl_range
        mask = (self.wavelength >= wl_min) & (self.wavelength <= wl_max)
        if not np.any(mask):
            return

        visible_data = self.data[:, mask]
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
        print(df)
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

        overlay.update({
            "data": df,
            "scale": float(scale),
            "file_path": file_path,
        })

        self._plot_or_refresh_time_trace_overlay(overlay_key)
        self._get_overlay_visible_var(overlay_key).set(True)
        self._get_overlay_toggle(overlay_key).config(state="normal")

        self.status_label.config(
            text=f"Status: {label} plotted with scale factor {float(scale):g}."
        )

    def _plot_or_refresh_time_trace_overlay(self, overlay_key):
        """Draw or redraw one stored absorbance overlay on the time-domain axis."""
        overlay = self.time_trace_overlays[overlay_key]
        df = overlay.get("data")
        scale = overlay.get("scale")
        if df is None or scale is None:
            return

        # Remove the previous artist before redrawing, especially after axis clears/range changes.
        old_line = overlay.get("line")
        if old_line is not None:
            try:
                old_line.remove()
            except ValueError:
                pass
            overlay["line"] = None

        wavelengths = df["Wavelength (nm)"].to_numpy(dtype=float)
        scaled_absorbance = float(scale) * df["Absorbance (AU)"].to_numpy(dtype=float)

        mask, wl_range = self._get_time_trace_wavelength_mask(wavelengths)
        if mask is None:
            return
        if not np.any(mask):
            messagebox.showwarning(
                "Range Warning",
                f"{overlay['label']} has no wavelengths inside the selected visible range."
            )
            return

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

        self._refresh_time_trace_legend()
        self.fig_time.tight_layout()
        self.canvas_time.draw()

    def _refresh_time_trace_overlays(self):
        """Redraw every loaded overlay after the time-domain axis is cleared or replotted."""
        for overlay_key in self.time_trace_overlays:
            self._plot_or_refresh_time_trace_overlay(overlay_key)

    def toggle_time_trace_overlay(self, overlay_key):
        """Show/hide a plotted absorbance overlay without deleting it."""
        overlay = self.time_trace_overlays[overlay_key]
        line = overlay.get("line")
        if line is None:
            return

        line.set_visible(self._get_overlay_visible_var(overlay_key).get())
        self._refresh_time_trace_legend()
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

    def plot_time_slices(self):
        """Plot ΔA vs wavelength for each requested time value (Tab 1)."""
        if not hasattr(self, 'data') or self.data is None:
            messagebox.showwarning("No Data", "Load a data file first.")
            return

        targets = self._parse_slice_values(self.time_slice_var.get(), "time")
        if targets is None:
            return

        self.ax_time.cla()
        self.ax_time.set_title("Time-Domain Traces")
        self.ax_time.set_xlabel("Wavelength (nm)")
        self.ax_time.set_ylabel("ΔA")
        self.ax_time.grid(True)
        self.ax_time.axhline(0, color='gray', lw=0.8, ls='--')

        cmap = plt.get_cmap("plasma")
        colors = [cmap(i / max(len(targets) - 1, 1)) for i in range(len(targets))]

        wl_range = None
        if getattr(self, "time_trace_wl_range_var", None) is not None:
            range_text = self.time_trace_wl_range_var.get().strip()
            if range_text:
                wl_range = self._parse_range_values(range_text, "wavelength")

        if wl_range is not None:
            wl_min, wl_max = wl_range
            plot_mask = (self.wavelength >= wl_min) & (self.wavelength <= wl_max)
        else:
            plot_mask = slice(None)

        for t_req, color in zip(targets, colors):
            idx = int(np.argmin(np.abs(self.times - t_req)))
            t_actual = self.times[idx]
            self.ax_time.plot(
                self.wavelength[plot_mask], self.data[idx, :][plot_mask],
                color=color, lw=1.5, label=f"{t_actual:.3g} ps"
            )

        if wl_range is not None:
            self.ax_time.set_xlim(*wl_range)
            self._autoscale_time_trace_yaxis(wl_range)
        else:
            self.ax_time.set_xlim(self.wavelength.min(), self.wavelength.max())
            self._autoscale_time_trace_yaxis((self.wavelength.min(), self.wavelength.max()))

        self._refresh_time_trace_overlays()
        self._refresh_time_trace_legend()
        self.fig_time.tight_layout()
        self.canvas_time.draw()
 
 
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
        if not hasattr(self, 'data') or self.data is None:
            messagebox.showwarning("No Data", "Load a data file first.")
            return

        parsed = self._parse_range_values(self.time_trace_wl_range_var.get(), "wavelength")
        if parsed is None:
            return

        wl_min, wl_max = parsed
        self.ax_time.set_xlim(wl_min, wl_max)
        self._autoscale_time_trace_yaxis((wl_min, wl_max))
        self._refresh_time_trace_overlays()
        self._refresh_time_trace_legend()
        self.fig_time.tight_layout()
        self.canvas_time.draw()
        self.status_label.config(
            text=f"Status: Time-trace wavelength view set to {wl_min:g} - {wl_max:g} nm."
        )


    def reset_time_trace_wavelength_range(self):
        """Reset the time-domain traces x-axis to the full wavelength span."""
        if self.wavelength is None:
            return
        self.time_trace_wl_range_var.set(f"{self.wavelength.min():.4g}, {self.wavelength.max():.4g}")
        self.ax_time.set_xlim(self.wavelength.min(), self.wavelength.max())
        self._autoscale_time_trace_yaxis((self.wavelength.min(), self.wavelength.max()))
        self._refresh_time_trace_overlays()
        self._refresh_time_trace_legend()
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
                self.times.min(),
                self.times.max()
            ],
            vmin=vmin,
            vmax=vmax
        )

        self.map_ax.set_title("2D Map View (TA Mean)")
        self.map_ax.set_xlabel("Wavelength (nm)")
        self.map_ax.set_ylabel("Time")

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

        # Using lambda allows us to destroy the window and call the function in one line
        ttk.Button(chooser, text="Correct Dispersions (MATLAB popup)", 
                   command=lambda: [chooser.destroy(), self._call_matlab_correctdispersion()]).pack(fill='x', padx=20, pady=4)
        ttk.Button(chooser, text="Apply External Dispersion Correction", 
                   command=lambda: [chooser.destroy(), self._open_apply_external_dialog()]).pack(fill='x', padx=20, pady=4)

        chooser.bind("<Escape>", lambda e: chooser.destroy())

    def _call_matlab_correctdispersion(self):
        """Call dispersion correction on the current working copy."""
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