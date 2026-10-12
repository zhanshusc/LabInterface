# Soliton processing inputs

These inputs come from `MatlabItems/Example_of_using_SolitonTAExperimentSelfContained_class.m`.
Use the example values only for that acquisition; other experiments need their
own calibration, background cutoff, and chirp fit.

## Acquisition and processing

| Input | Example | Interface input |
| --- | --- | --- |
| First scan and data folder | `scan1_unref_DOD.txt`, in the supplied Pyrazine mirrored-scan folder | Select this file with Load TA Data; the folder is inferred. |
| Number of scans | `1` | Number of scan files to combine. Multiple scans must use matching names with `scan1`, `scan2`, etc. and matching matrix shapes. |
| Mirrored acquisition | `1` (yes) | Yes for forward and reverse delay sweeps in each file. The row count must be even. |
| Background upper bound | `-1.85` ps | A pre-signal delay within the measured range. The average includes all sorted delays through the nearest measured delay to this cutoff. |
| Energy calibration A, E, T0 | `276390, 12.4698, 1320.4` | Optional for pixel plots; required for calibrated energy/wavelength plots. |
| Chirp coefficients | `a4, a3, a2, a1, a0` | Five experiment-specific values, highest power first. Import with zeros to leave the data unshifted, then fit using Dispersion Correction if needed. |

Text files contain a whitespace-separated numeric matrix: delay in ps in the
first column, then TA signal for each detector channel. The example uses 1024
channels; the interface infers the count from the file.

Energy calibration is `sqrt(E^2 - A / (pixel + T0))`, in eV. Pixels start at 1.
Invalid and nonpositive energies are hidden from physical-axis plots, while
original detector numbers are retained for chirp fitting.

The chirp delay is `a4*pixel^4 + a3*pixel^3 + a2*pixel^2 + a1*pixel + a0`, in ps.
It is subtracted from each channel's time axis. Linear interpolation maps the
shifted signal back onto the common delay grid; outside-range values become zero.
Interactive fits of orders 1–4 are padded to five coefficients with leading zeros.

The MATLAB example contains no numerical nonzero chirp coefficients. It copies
the coefficients from its first processing call. In the included MATLAB and
Python implementations, the default is all zeros and no automatic fit occurs.
Supply a calibration fit or select time-zero points in the interface.

## Background subtraction

Processing order: load → unmirror if needed → sort and remove duplicate delays
→ average scans/sweeps → subtract background → correct chirp.

The background is a separate average for each channel. It is subtracted from
every time point in the sorted, unshifted scan mean. Changing or repeating the
cutoff always starts from this mean, so subtraction does not accumulate.
The background window averages to zero before chirp correction. After shifting
and interpolation, the same fixed window need not average to zero.

The standard deviation uses population normalization, matching MATLAB
`std(data, 1)`. The interface reports the actual cutoff and number of delays used.

## Optional spectra and metadata

The example also imports unirradiated and irradiated absorbance CSVs, with
wavelength, absorbance, and an optional variance column, plus notes about dilution
and path length. Fluorescence spectra and overlay scaling are optional. These
spectra are not required for background subtraction or chirp correction.

The following metadata are assigned separately and do not drive the core
processing. The current TA import dialog does not collect these fields.

| Metadata | Example value | Units |
| --- | --- | --- |
| Pump wavelength | 267 | nm |
| Pump bandwidth | 1000 | nm FWHM, as written; confirm for your experiment |
| Pump energy | 75 | nJ |
| Optical density | 0.46 | OD |
| Jet path length | 125 | µm |
| Pump/probe angle | 15 | degrees; explicitly a placeholder in the example |
| Pump spot size | 107 | µm FWHM at sample |
| Probe spot size | 1000 | µm FWHM at sample |
| Pump polarization | 35 | degrees in the lab frame |
| Probe polarization | 90 | degrees in the lab frame |

Use measured values for your experiment. Horizontal polarization is 0° and
vertical is 90°. The example's pump-energy assignment uses `pryTA` instead of
`pyrTA`; correct that object name if reusing its MATLAB assignments.
