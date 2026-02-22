import tkinter as tk
from tkinter import ttk
from tkinter import messagebox  
from tkinter import filedialog
from tkinter import simpledialog
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import matlab.engine
import threading

# Using messagebox for placeholders

# Import Matplotlib libraries for embedding plots


class PsTAAnalysisApp(tk.Tk):
    """
    A Python tkinter application replicating the psTA Analysis Suite layout
    described for MATLAB App Designer.
    """
    def __init__(self):
        super().__init__()

        # --- Main App Window Setup ---
        self.title("psTA Analysis Suite (Python)")
        self.geometry("1200x700")
        
        # --- Main Grid Layout (3 Columns) ---
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

        # --- Create the Three Major Panels ---
        # We use ttk.Frame for a modern look
        self.left_panel = ttk.Frame(self, padding="10")
        self.middle_panel = ttk.Frame(self, padding="10")
        self.right_panel = ttk.Frame(self, padding="10")

        # Place panels on the grid
        self.left_panel.grid(row=0, column=0, sticky="nsew")
        self.middle_panel.grid(row=0, column=1, sticky="nsew")
        self.right_panel.grid(row=0, column=2, sticky="nsew")

        # --- Populate Each Panel ---
        self._create_left_panel()
        self._create_middle_panel()
        self._create_right_panel()

        # --- Bottom Status Bar ---
        self._create_status_bar()

        # --- Data Containers ---
        # Data is the container for the primary data - ie experiment being analzyed 
        # data 2 is a place holder for TCSPC for future use
        # 1 set of wavelength and times place holder for now, if need be duplica
        self.data = None
        self.data2 = None
        self.wavelength = None
        self.times = None

        # --- Matlab Environment Start ---
        self.eng = matlab.engine.start_matlab()

        # Add path to MATLAB script
        self.eng.addpath(r"./TAExperiment.m", nargout=0)
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.map_colorbar = None

    

    def _create_left_panel(self):
        """Populates the Left Panel (Data Control & Processing)"""
        frame = self.left_panel
        
        # --- Section: Data I/O ---
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

        # --- Section: Preprocessing ---
        preproc_frame = ttk.LabelFrame(frame, text="Preprocessing", padding="10")
        preproc_frame.pack(fill='x', pady=5)

        # Background Correction controls (user input + run)
        ttk.Button(
            preproc_frame, text="Background Correction", command=self.background_correction
        ).pack(fill='x', pady=2)
        
        ttk.Button(
            preproc_frame, text="Dispersion Correction", command=self.dispersion_correction
        ).pack(fill='x', pady=2)

        # --- Sub-section: Outlier Removal ---
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
        # --- End Outlier Sub-section ---

        ttk.Button(
            preproc_frame, text="Recompute Average", command=self.placeholder_command
        ).pack(fill='x', pady=5)
        

    def _create_middle_panel(self):
        """Populates the Middle Panel (Visualization & Diagnostics)"""
        frame = self.middle_panel
        
        # Create the Tab Group (Notebook)
        tab_control = ttk.Notebook(frame)
        
        # --- Create Tab Frames ---
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

        # --- Populate Tab 1: Time-Domain Traces ---
        # Add controls
        tab1_controls = ttk.Frame(tab1)
        tab1_controls.pack(fill='x', pady=5)
        
        ttk.Label(tab1_controls, text="Wavelength:").pack(side='left', padx=5)
        self.trace_wavelength_combo = ttk.Combobox(
            tab1_controls, 
            values=["(no data loaded)", "500", "520", "550"]
        )
        self.trace_wavelength_combo.current(0)
        self.trace_wavelength_combo.pack(side='left', padx=5)
        
        self.show_scans_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            tab1_controls, text="Show individual scans", variable=self.show_scans_var
        ).pack(side='left', padx=5)
        
        self.show_mean_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            tab1_controls, text="Show mean trace", variable=self.show_mean_var
        ).pack(side='left', padx=5)
        
        # Add Matplotlib Plot Canvas (UIAxes_Traces)
        trace_fig = Figure(figsize=(5, 4), dpi=100)
        trace_ax = trace_fig.add_subplot(111)
        trace_ax.set_title("Time-Domain Traces")
        trace_ax.set_xlabel("Time (ps/ns)")
        trace_ax.set_ylabel("ΔA")
        trace_ax.grid(True)
        trace_ax.text(0.5, 0.5, "UIAxes_Traces", horizontalalignment='center', verticalalignment='center', transform=trace_ax.transAxes, fontsize=16, color='gray', alpha=0.5)

        canvas_tab1 = FigureCanvasTkAgg(trace_fig, master=tab1)
        canvas_tab1.draw()
        canvas_tab1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


         # --- Populate Tab 1: Time-Domain Traces ---
        # Add controls
        tab4_controls = ttk.Frame(tab4)
        tab4_controls.pack(fill='x', pady=5)
        
        ttk.Label(tab4_controls, text="Wavelength:").pack(side='left', padx=5)
        self.trace_wavelength_combo = ttk.Combobox(
            tab4_controls, 
            values=["(no data loaded)", "500", "520", "550"]
        )
        self.trace_wavelength_combo.current(0)
        self.trace_wavelength_combo.pack(side='left', padx=5)
        
        self.show_scans_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            tab4_controls, text="Show individual scans", variable=self.show_scans_var
        ).pack(side='left', padx=5)
        
        self.show_mean_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            tab4_controls, text="Show mean trace", variable=self.show_mean_var
        ).pack(side='left', padx=5)
        
        # Add Matplotlib Plot Canvas (UIAxes_Traces)
        trace_fig = Figure(figsize=(5, 4), dpi=100)
        trace_ax = trace_fig.add_subplot(111)
        trace_ax.set_title("Wavelength-Domain Traces")
        trace_ax.set_xlabel("Wavelengths / nm")
        trace_ax.set_ylabel("ΔA")
        trace_ax.grid(True)
        trace_ax.text(0.5, 0.5, "UIAxes_Traces", horizontalalignment='center', verticalalignment='center', transform=trace_ax.transAxes, fontsize=16, color='gray', alpha=0.5)

        canvas_tab4 = FigureCanvasTkAgg(trace_fig, master=tab4)
        canvas_tab4.draw()
        canvas_tab4.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # --- Populate Tab 3: 2D Map View ---
        tab2_controls = ttk.Frame(tab2)
        tab2_controls.pack(fill='x', pady=5)

        ttk.Button(
            tab2_controls, text="Normalize", command=self.placeholder_command
        ).pack(side='left', padx=5)
        
        self.colorbar_toggle_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            tab2_controls, text="Show Colorbar", variable=self.colorbar_toggle_var
        ).pack(side='left', padx=5)

        # Add Matplotlib Plot Canvas (UIAxes_Map)
        self.map_fig = Figure(figsize=(5, 4), dpi=100)
        self.map_ax = self.map_fig.add_subplot(111)
        self.map_ax.set_title("2D Map View (dat.TAmean)")
        self.map_ax.set_xlabel("Wavelength (nm)" )
        self.map_ax.set_ylabel("Time (log scale)")
        self.map_ax.text(0.5, 0.5, "UIAxes_Map", horizontalalignment='center', verticalalignment='center', transform=self.map_ax.transAxes, fontsize=16, color='gray', alpha=0.5)

        self.canvas_tab2 = FigureCanvasTkAgg(self.map_fig, master=tab2)
        self.canvas_tab2.draw()
        self.canvas_tab2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # --- Populate Tab 4: Residuals / Cleaned Comparison ---
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

        # --- Section: Kinetic Fitting ---
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
        
        # --- Section: Export Options ---
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

        self.status_label.config(text="Status: Running MATLAB analysis...")
        self.update_idletasks()

        try:
            self.run_matlab_analysis(file_path)
            self.status_label.config(text="Status: Analysis complete.")
        except Exception as e:
            messagebox.showerror("MATLAB Error", str(e))
            self.status_label.config(text="Status: Error during analysis.")

    def run_matlab_analysis(self, data_file):

        # Run MATLAB function
        dat = self.eng.TAExperiment(data_file)
        self.eng.workspace['dat'] = dat 
        # Extract arrays
        arr = np.array(self.eng.eval("dat.TAMean"))
        times = np.array(self.eng.eval("dat.times")).flatten()
        wavelengths = np.array(self.eng.eval("dat.wavelengths")).flatten()

        # Store for later use
        self.data = arr
        self.times = times
        self.wavelength = wavelengths

        # Update the 2D map
        self.update_2d_map()

    def update_2d_map(self):
        self.map_fig.clear()
        
        # Re-add the main axis
        self.map_ax = self.map_fig.add_subplot(111)
        
        # Plot data
        im = self.map_ax.imshow(
            self.data,
            aspect='auto',
            origin='lower',
            extent=[
                self.wavelength.min(),
                self.wavelength.max(),
                self.times.min(),
                self.times.max()
            ]
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
        Ask user for upper bound after button press,
        then call subtractTABackground(dat, upper_bound)
        """

        if self.data is None:
            messagebox.showwarning("No Data", "Load data before running background correction.")
            return

        # Ask user for float input
        upper_val = simpledialog.askfloat(
            "Background Correction",
            "Enter upper bound for background subtraction:",
            parent=self
        )
        # User pressed cancel
        if upper_val is None:
            return
        self.status_label.config(text=f"Status: Running background correction (upper={upper_val})...")
        self.update_idletasks()
        try:
            # Get current MATLAB dat
            dat = self.eng.workspace['dat']
            # If MATLAB function RETURNS updated dat
            self.eng.subtractTABackground(dat, float(upper_val), nargout=0)
            corrected_dat = self.eng.workspace['dat']

            # Store back in workspace
            self.eng.workspace['dat'] = corrected_dat

            # Extract updated arrays
            arr = np.array(self.eng.eval("dat.TAMean"))
            times = np.array(self.eng.eval("dat.times")).flatten()
            wavelengths = np.array(self.eng.eval("dat.wavelengths")).flatten()

            # Update Python-side data
            self.data = arr
            self.times = times
            self.wavelength = wavelengths

            # Refresh 2D map
            self.update_2d_map()

            self.status_label.config(text="Status: Background correction complete.")
            messagebox.showinfo("Success", "Background subtraction finished successfully.")

        except Exception as e:
            messagebox.showerror("MATLAB Error", str(e))
            self.status_label.config(text=f"Status: Error during background correction.{e}")

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
        """Call dat.correctdispersion() or correctdispersion(dat) in MATLAB."""
        self._update_status("Launching MATLAB dispersion-correction popup...")
        try:
            # Let MATLAB handle the method vs function resolution natively
            cmd = "dat.correctDispersion();"
            self.eng.eval(cmd, nargout=0)
            
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
                # Parse numbers: accept commas and/or spaces
                parts = [p for p in txt.replace(",", " ").split() if p]
                if len(parts) != 3:
                    raise ValueError("Please enter exactly 3 numbers.")
                coeffs = [float(p) for p in parts]
            except Exception as ex:
                status_lbl.config(text=f"Invalid input: {ex}")
                return

            dlg.destroy()
            
            # Apply to MATLAB
            self._update_status("Applying external dispersion correction...")
            try:
                vec_str = f"[{' '.join(map(str, coeffs))}]"
                self.eng.eval(f"applyExternalDispersionCorrection(dat, {vec_str});", nargout=0)
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
        """Pull dat arrays from MATLAB workspace and update Python state."""
        try:
            self.data = np.array(self.eng.eval("dat.TAMean"))
            self.times = np.array(self.eng.eval("dat.times")).flatten()
            self.wavelength = np.array(self.eng.eval("dat.wavelengths")).flatten()

            # Refresh map if the method exists
            if hasattr(self, 'update_2d_map'):
                self.update_2d_map()
                
            self._update_status(success_msg)
            messagebox.showinfo("Done", success_msg)
            
        except Exception as e:
            messagebox.showerror("MATLAB Error", f"Could not refresh dat from MATLAB workspace: {e}")


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