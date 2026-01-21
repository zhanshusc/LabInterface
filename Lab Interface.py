import tkinter as tk
from tkinter import ttk
from tkinter import messagebox  
from tkinter import filedialog
# Using messagebox for placeholders

# Import Matplotlib libraries for embedding plots
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

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

        ttk.Button(
            preproc_frame, text="Background Correction", command=self.placeholder_command
        ).pack(fill='x', pady=2)
        
        ttk.Button(
            preproc_frame, text="Dispersion Correction", command=self.placeholder_command
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

        # --- Populate Tab 2: 2D Map View ---
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
        map_fig = Figure(figsize=(5, 4), dpi=100)
        map_ax = map_fig.add_subplot(111)
        map_ax.set_title("2D Map View (dat.TAmean)")
        map_ax.set_xlabel("Wavelength (nm)" )
        map_ax.set_ylabel("Time (log scale)")
        map_ax.text(0.5, 0.5, "UIAxes_Map", horizontalalignment='center', verticalalignment='center', transform=map_ax.transAxes, fontsize=16, color='gray', alpha=0.5)

        canvas_tab2 = FigureCanvasTkAgg(map_fig, master=tab2)
        canvas_tab2.draw()
        canvas_tab2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # --- Populate Tab 3: Residuals / Cleaned Comparison ---
        # This tab will have two plots, so we use a main frame
        tab3_main_frame = ttk.Frame(tab3)
        tab3_main_frame.pack(fill=tk.BOTH, expand=True)
        tab3_main_frame.grid_rowconfigure(0, weight=1)
        tab3_main_frame.grid_rowconfigure(1, weight=1)
        tab3_main_frame.grid_columnconfigure(0, weight=1)
        
        # Top Plot: Before
        before_fig = Figure(dpi=100)
        before_ax = before_fig.add_subplot(111)
        before_ax.set_title("psTA")
        before_ax.set_ylabel("ΔA")
        before_ax.grid(True)
        before_ax.text(0.5, 0.5, "psTA", horizontalalignment='center', verticalalignment='center', transform=before_ax.transAxes, fontsize=14, color='gray', alpha=0.5)
        
        canvas_tab3_top = FigureCanvasTkAgg(before_fig, master=tab3_main_frame)
        canvas_tab3_top.draw()
        canvas_tab3_top.get_tk_widget().grid(row=0, column=0, sticky="nsew", pady=2)

        # Bottom Plot: After
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
            """
            Opens a file dialog for the user to select a file.
            """
            # 1. Open the file explorer
            file_path = filedialog.askopenfilename(
                title="Select psTA Data File",
                # Define which files show up. Change "*.csv" to your specific extension if needed
                filetypes=[
                    ("All Files", "*.*")
                ]
            )
            # 2. Check if the user actually selected a file (didn't click Cancel)
            if file_path:
                print(f"File selected: {file_path}")
                
                # 3. Load the data (Example using standard file reading)
                self.data = file_path

if __name__ == "__main__":
    # Set high-DPI awareness for Windows
    try:
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    except:
        pass # Fails on non-Windows systems
    app = PsTAAnalysisApp()
    app.mainloop()