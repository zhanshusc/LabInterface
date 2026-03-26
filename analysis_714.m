clear; close all; clc;
%11/15/2023
%%

dat = TAExperiment('PZI-Au-Bnc2Cz-400ex_0p3od-150nj_thf');
data= TAExperiment('MeCHX bLANK');

subtractTABackground(dat,-0.9);
subtractTABackground(data,-.5);

%applyExternalDispersionCorrection(data1,[1.03417442372052e-08,-1.23050580408648e-05,0.00523735124076471,-0.507337168663455]);

dat.wavelengths=dat.wavelengths+3;
dat.times=dat.times-0.3;

dat.correctDispersion();
applyExternalDispersionCorrection (dat,[-7.44208608527694e-06	0.00672243358275736	-0.951026821920219]);
applyExternalDispersionCorrection (data,[-7.07742834499263e-06	0.00629237385777828	-0.986355049838935]);

peekTA(dats);
scans(data1);
data1.TAMean=movmean(data1.TAMean,5);
peekTA(dat);

TAscans_Mac=data1.TAMean;

dat.peekTimeTraces([520,550,490,430]);
peekSpectra(dat,[0.5,2,10,50,100,500,1000]);
ylim([-0.005 0.004]);
xlim([350 720]);

peekTA(dat);
colormap('jet');
clim([-0.01 0.003051]);
ylim([-1 1000]);

dat.peekTimeTraces([450,500,550,600,700]);
peekSpectra(dat,[-1,-0.6,-0.2,]);

xlim([415 720]);
ylim([-0.001 0.006]);

peekSpectra(dat,[1]);
peekSpectra(dat,[0.6]);

%% 
% Define wavelength range
wavelength_min = 415; wavelength_max = 720;

% Extract indices for the selected wavelength range
wavelength_indices = find(dat.wavelengths >= wavelength_min & dat.wavelengths <= wavelength_max);

% Create the concatenated data matrix with the selected wavelengths
data_concat = [0, dat.wavelengths(wavelength_indices)];
data_concat(2:size(dat.times, 1)+1, 1) = dat.times;
data_concat(2:end, 2:end) = dat.TAMean(:, wavelength_indices);

% Save the data to a file
save('MD24h-MeCHX150nj400nmex12112024.txt', 'data_concat', '-ascii');
%% 
% Enter the target wavelength
target_wavelength = 560; % Specify wavelength here

% Extract time axis
time = dat.times(:); % Ensure time is a column vector

% Find the index of the closest wavelength
[~, idx] = min(abs(dat.wavelengths - target_wavelength));

% Select indices within the time range of 1 to 900
time_mask = find((dat.times >= 1) & (dat.times <= 900));
time = dat.times(time_mask);
intensity = dat.TAMean(time_mask, idx); % Extract intensity for the selected time range


%% % Plot the data and the exponential decay fit
% --- Inputs ---
num_exp = 2; target_wavelength = 560;
time = time(:); y = intensity(:) / max(intensity);  % Normalize

% --- Fit Model ---
if num_exp == 1
    model = @(b, t) b(1)*exp(-t/b(2));  init = [1, 100];
else
    model = @(b, t) b(1)*exp(-t/b(2)) + b(3)*exp(-t/b(4));  init = [0.5, 50, 0.5, 1000];
end
opts = optimset('Display','off');
b = lsqcurvefit(model, init, time, y, [], [], opts);

% --- Lifetimes Label ---
if num_exp == 1
    label = sprintf('\\tau = %.1f ps', b(2));
else
    label = sprintf('\\tau_1 = %.1f ps, \\tau_2 = %.1f ps', b(2), b(4));
end

% --- Plot ---
figure;
plot(time, y, 'ko', 'MarkerSize', 4); hold on;
plot(time, model(b, time), 'r-', 'LineWidth', 2);
legend('Data', label);
xlabel('Time (ps)'); ylabel('Norm. Intensity');
title(sprintf('%d-Exp Fit at %d nm', num_exp, target_wavelength));
grid on; box on;


%% 
% Define the time points (in same units as dat.times)
t1 = 100;   % e.g., 1.5 ps
t2 = 800;   % e.g., 300 ps

% Define the normalization wavelength
norm_wavelength =  550;  % nm

% Find closest indices for the times and normalization wavelength
[~, idx1] = min(abs(dat.times - t1));
[~, idx2] = min(abs(dat.times - t2));
[~, norm_idx] = min(abs(dat.wavelengths - norm_wavelength));

% Extract the slices
slice1 = dat.TAMean(idx1, :);
slice2 = dat.TAMean(idx2, :);

% Normalize using the value at 550 nm
norm_slice1 = slice1 / slice1(norm_idx);
norm_slice2 = slice2 / slice2(norm_idx);

% Plot
figure;
plot(dat.wavelengths, norm_slice1, 'LineWidth', 2); hold on;
plot(dat.wavelengths, norm_slice2, '--', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel(sprintf('Normalized ΔA (at %.0f nm)', norm_wavelength));
legend(sprintf('%.2f ps', dat.times(idx1)), sprintf('%.2f ps', dat.times(idx2)));
title('Overlay of Two Normalized Time Slices');
xlim([420 720]);
ylim([0 1.1]);
grid on;
box on;
%% 
% Target time (in ps) and normalization wavelength (in nm)
target_time = 900; 
norm_wavelength = 700;

% Find nearest time index for 1000 ps
[~, idx_time_data] = min(abs(data.times - target_time));
[~, idx_time_dat]  = min(abs(dat.times - target_time));

% Find nearest wavelength index for 550 nm
[~, idx_wl_data] = min(abs(data.wavelengths - norm_wavelength));
[~, idx_wl_dat]  = min(abs(dat.wavelengths - norm_wavelength));

% Extract the time slices
slice_data = data.TAMean(idx_time_data, :);
slice_dat  = dat.TAMean(idx_time_dat, :);

% Normalize slices using the value at 550 nm
slice_data_norm = slice_data / slice_data(idx_wl_data);
slice_dat_norm  = slice_dat  / slice_dat(idx_wl_dat);

% Plotting
figure;
plot(data.wavelengths, slice_data_norm, 'b-', 'LineWidth', 2); hold on;
plot(dat.wavelengths,  slice_dat_norm,  'r--', 'LineWidth', 2);
xlim([420 720]);
ylim([0 1.1]);
xlabel('Wavelength (nm)');
ylabel('Normalized ΔA');
legend('data @ 1000 ps','dat @ 1000 ps');
title('Overlay of 1000 ps Slices from data and dat');
grid on; box on;
%% 
% Target time (in ps) and normalization wavelength (in nm)
target_time = 1000; 
norm_wavelength = 550;

% Smoothing window (number of points)
smooth_window = 5;

% Find nearest time index for 1000 ps
[~, idx_time_data] = min(abs(data.times - target_time));
[~, idx_time_dat]  = min(abs(dat.times - target_time));

% Find nearest wavelength index for 550 nm
[~, idx_wl_data] = min(abs(data.wavelengths - norm_wavelength));
[~, idx_wl_dat]  = min(abs(dat.wavelengths - norm_wavelength));

% Extract the time slices
slice_data = data.TAMean(idx_time_data, :);
slice_dat  = dat.TAMean(idx_time_dat, :);

% Apply smoothing (moving average)
slice_data_smooth = movmean(slice_data, smooth_window);
slice_dat_smooth  = movmean(slice_dat,  smooth_window);

% Normalize using the value at 550 nm (after smoothing)
slice_data_norm = slice_data_smooth / slice_data_smooth(idx_wl_data);
slice_dat_norm  = slice_dat_smooth  / slice_dat_smooth(idx_wl_dat);

% Plotting
figure;
plot(data.wavelengths, slice_data_norm, 'b-', 'LineWidth', 2); hold on;
plot(dat.wavelengths,  slice_dat_norm,  'r--', 'LineWidth', 2);
xlim([420 720]);
ylim([0 1.1]);
xlabel('Wavelength (nm)');
ylabel('Normalized ΔA');
legend('data @ 1000 ps','dat @ 1000 ps');
title('Overlay of Smoothed 1000 ps Slices from data and dat');
grid on; box on;
%% 
% Load data
sads = readmatrix('MAC-Au-BnBCz-THF-07152024 425 to 720.txt_amplitudes.dat');

% Define which columns to include as Y
y_cols = [2,3,4,5];

% Extract X and create mask for 410–720 nm range
x = sads(:, 1);
x_mask = (x >= 415) & (x <= 720);
x_filtered = x(x_mask);

% Plotting
figure;
hold on;

% Loop through each specified Y column and plot
for i = 1:length(y_cols)
    y = sads(:, y_cols(i));
    y_filtered = y(x_mask);
    plot(x_filtered, y_filtered, 'LineWidth', 2);
end

% Formatting
xlabel('Wavelength (nm)');
ylabel('Amplitude');
legend(arrayfun(@(c) sprintf('Col %d', c), y_cols, 'UniformOutput', false));
title('Selected Y Columns vs Wavelength');
box on;
set(gca, 'XGrid', 'off', 'YGrid', 'off'); % Disable grid lines

%% 
% Load data
sads = readmatrix('MAC-aU-bNCZ THF 410-720 NM 07152024_amplitudes.dat');

% Define which columns to include as individual Y plots
y_cols = [2,3,4];

% Define two columns to be summed and plotted together
sum_cols = [4, 5];  % Change these to the columns you want to sum

% Extract X and apply 410–720 nm mask
x = sads(:, 1);
x_mask = (x >= 410) & (x <= 720);
x_filtered = x(x_mask);

% Plotting
figure;
hold on;

% Plot individual y_cols
for i = 1:length(y_cols)
    y = sads(:, y_cols(i));
    y_filtered = y(x_mask);
    plot(x_filtered, y_filtered, 'LineWidth', 2);
end

% Plot sum of specified columns
y_sum = sads(:, sum_cols(1)) + sads(:, sum_cols(2));
y_sum_filtered = y_sum(x_mask);
plot(x_filtered, y_sum_filtered, 'k--', 'LineWidth', 2); % Black dashed line for sum

% Formatting
xlabel('Wavelength (nm)');
ylabel('Amplitude');
legend_labels = [arrayfun(@(c) sprintf('Col %d', c), y_cols, 'UniformOutput', false), ...
                 {sprintf('Sum of Col %d + Col %d', sum_cols(1), sum_cols(2))}];
legend(legend_labels);
title('Selected Y Columns vs Wavelength');
box on;
set(gca, 'XGrid', 'off', 'YGrid', 'off'); % Disable grid lines

%% 
% --- User-defined wavelengths (nm) ---
wavelengths_to_plot = [460, 700];  % Modify these as needed

% --- Extract axes and data from 'dat' ---
TA = dat.TAMean;              % [time x wavelength]
Time = dat.times;             % Time axis (ps)
ProbeAxis = dat.wavelengths;  % Wavelength axis (nm)

% --- Find indices of the nearest wavelengths ---
[~, w_idx1] = min(abs(ProbeAxis - wavelengths_to_plot(1)));
[~, w_idx2] = min(abs(ProbeAxis - wavelengths_to_plot(2)));

% --- Extract decay traces ---
trace1 = TA(:, w_idx1);
trace2 = TA(:, w_idx2);

% --- Apply smoothing (moving average or sgolay) ---
% You can change the window size (e.g., 11) and method ('moving', 'sgolay', etc.)
smooth_trace1 = smooth(trace1, 11, 'moving');
smooth_trace2 = smooth(trace2, 11, 'moving');

% --- Plotting ---
figure('Color', 'w'); hold on;
plot(Time, smooth_trace1, 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('%.0f nm (smoothed)', ProbeAxis(w_idx1)));
plot(Time, smooth_trace2, 'b--', 'LineWidth', 2, ...
    'DisplayName', sprintf('%.0f nm (smoothed)', ProbeAxis(w_idx2)));

legend('show', 'Location', 'northeast', 'FontSize', 12);
xlabel('Time (ps)', 'FontSize', 14);
ylabel('ΔA (smoothed)', 'FontSize', 14);
title('Smoothed psTA Traces at Two Wavelengths', 'FontSize', 16);
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');
grid on;

%% 

% Target wavelength
target_wavelength = 530; % nm

% Time window (ps)
time_window = 0.5;  % ±0.5 ps

% Extract time axis
time = dat.times(:); % Ensure column vector

% Find index of closest wavelength
[~, idx] = min(abs(dat.wavelengths - target_wavelength));

% Find indices within ±time_window
time_mask = find(time >= -time_window & time <= time_window);

% Extract data and convert to mOD
time_selected = time(time_mask);
intensity_selected = dat.TAMean(time_mask, idx) * 1000; % Convert to mOD

% --- Find the highest signal point ---
[~, max_idx] = max(intensity_selected);
t_peak = time_selected(max_idx);

% --- Shift time axis only (keep positive intensity) ---
time_centered = time_selected - t_peak;

% --- Gaussian Fit ---
gaussEqn = 'a*exp(-((x-b)^2)/(2*c^2)) + d';
startPoints = [max(intensity_selected), 0, 0.1, min(intensity_selected)];
fit_gauss = fit(time_centered, intensity_selected, gaussEqn, 'Start', startPoints);

% Generate smooth curve for plotting
x_fit = linspace(min(time_centered), max(time_centered), 400);
y_fit = feval(fit_gauss, x_fit);

% --- Calculate FWHM in femtoseconds with error ---
sigma = fit_gauss.c;
FWHM_ps = 2 * sqrt(2 * log(2)) * sigma;
FWHM_fs = FWHM_ps * 1000;  % ps → fs

% Compute error on sigma from confidence interval
ci = confint(fit_gauss, 0.95); % 95% confidence
sigma_err = (ci(2,3) - ci(1,3)) / 2; % error on c
FWHM_err_fs = 2 * sqrt(2 * log(2)) * sigma_err * 1000;

% --- Plot Data and Fit ---
figure('Color','w');
hold on;
p1 = plot(time_centered, intensity_selected, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5, ...
    'DisplayName', 'Data');
p2 = plot(x_fit, y_fit, 'r-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('FWHM = %.1f ± %.1f fs', FWHM_fs, FWHM_err_fs));

% Axis labels
xlabel('Time (ps)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('-ΔAbs. (mOD)', 'FontSize', 16, 'FontWeight', 'bold');

% Fix x-axis limits
xlim([-0.3 0.3]);  % Always limit between -0.3 ps and +0.3 ps

% Legend
legend([p1, p2], 'FontSize', 14, 'Location', 'best', 'Box', 'off');

% Grid and axes formatting
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'on', ...
         'TickDir', 'out', 'TickLength', [0.015, 0.015]);

% Export quality
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4]); % for saving as high-res figure

hold off;

% --- Display fit parameters and FWHM ---
disp('Gaussian Fit Parameters (after peak centering):');
disp(fit_gauss);
fprintf('Calculated FWHM = %.2f ± %.2f fs\n', FWHM_fs, FWHM_err_fs);
fprintf('Data shifted by %.4f ps so that max signal is at 0 ps.\n', t_peak);
%% time points in different scans. give the scan no and time points 
% ---- inputs you choose ----
scanIdx    = [1 4];                 % scans to plot
timePoints = [2 980];      % times (same units as dat.times)

tt = dat.times(:);
wl = dat.wavelengths(:).';

if numel(scanIdx) > 1
    tiledlayout(numel(scanIdx),1,'TileSpacing','compact');
end

for k = 1:numel(scanIdx)
    s = scanIdx(k);
    if iscell(dat.scans), S = dat.scans{1,s}; else, S = dat.scans(s); end
    M = S.TAMean;  % time x wavelength

    [~, idxT] = arrayfun(@(t) min(abs(tt - t)), timePoints(:));
    tUsed = tt(idxT);

    if numel(scanIdx) > 1, nexttile; end
    hold on
    try, colororder(turbo(numel(idxT))); catch, colororder(parula(numel(idxT))); end
    for j = 1:numel(idxT)
        plot(wl, M(idxT(j),:), 'LineWidth', 1.6);
    end
    xlabel('Wavelength'); ylabel('\DeltaA'); title(sprintf('Scan %d', s));
    legend(compose('t = %g', tUsed), 'Location','best'); box on
end
%% 
% Read two CSV files (each should have columns: wavelength, absorbance)-
% plot for UV vis before and after
data1 = readmatrix('Zr-Z-TolueneAcquisition 1 2025-09-10 16_02_10 »» Detector1 »» D1 345_350-750.trace.txt');
data2 = readmatrix('Zr-Z-TolueneAcquisition 1 2025-09-10 16_02_10 »» Detector1 »» D1 340_350-750.trace.txt');

% Extract columns
wl1 = data1(:,1); A1 = data1(:,2);
wl2 = data2(:,1); A2 = data2(:,2);

% Plot overlay
figure; hold on
plot(wl1, A1, 'b-', 'LineWidth', 1.6);
plot(wl2, A2, 'r-', 'LineWidth', 1.6);

xlabel('Wavelength (nm)');
ylabel('Absorbance');
legend('Spectrum 1','Spectrum 2');
title('Overlay of UV–Vis Spectra');
box on; grid on

