"""Transient-absorption experiment model, ported from ``MatlabItems/TAExperiment.m``.

The class loads Burritos JSON acquisitions and Specter 7.16 text exports.  Its
public arrays use the same orientation as the MATLAB class: rows are time
delays and columns are spectrometer pixels.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np


DataType = Literal["TA", "Pump On", "Pump Off"]


@dataclass
class CalibrationPoint:
    wavelength: float
    pixel: float
    uncertainty: float


@dataclass
class TAScan:
    """Per-scan statistics used to make the experiment-wide weighted mean."""

    TAMean: np.ndarray
    TAVariance: np.ndarray
    TANShots: np.ndarray
    pumpOnMean: np.ndarray
    pumpOnVariance: np.ndarray
    pumpOnNShots: np.ndarray
    pumpOffMean: np.ndarray
    pumpOffVariance: np.ndarray
    pumpOffNShots: np.ndarray


class TAExperiment:
    """Load, inspect, and process one transient-absorption experiment.

    Parameters
    ----------
    file_name:
        Burritos JSON acquisition (default) or Specter 7.16 text export.
    experiment:
        ``"Burritos"`` or ``"Specter"``.  The spelling is case-insensitive.
    """

    def __init__(self, file_name: str | Path, experiment: str = "Burritos") -> None:
        self.fileName = str(file_name)
        path = Path(file_name)
        if not path.is_file():
            raise FileNotFoundError(f"File not found: {path}")

        self.json: dict[str, Any] | None = None
        self.CalibrationPoints: list[CalibrationPoint] = []
        self.analytes: list[Any] = []
        self.solvent = ""
        self.pumpWavelength = self.pumpBandwidth = self.pumpEnergy = 0.0
        self.opticalDensity = self.jetPathLength = self.pumpAngle = 0.0
        self.pumpSpotSize = self.probeSpotSize = 0.0
        self.dispersionFitCoefficients = np.array([], dtype=float)
        self.dispersionFit = np.array([], dtype=float)
        self.scans: list[TAScan] = []

        kind = experiment.lower()
        if kind == "burritos":
            self.loadFile()
        elif kind == "specter":
            self.loadFile716()
        else:
            raise ValueError("experiment must be 'Burritos' or 'Specter'")

    @staticmethod
    def _array(value: Any, *, shape: tuple[int, int] | None = None) -> np.ndarray:
        array = np.asarray(value, dtype=float)
        if shape is not None:
            if array.size != shape[0] * shape[1]:
                raise ValueError(f"Expected {shape[0]} by {shape[1]} values; got {array.size}.")
            array = array.reshape(shape)
        return array

    @staticmethod
    def _field(source: dict[str, Any], *names: str, default: Any = None) -> Any:
        """Return the first available spelling used by acquisition versions."""
        for name in names:
            if name in source:
                return source[name]
        return default

    def _empty_scan(self, n_times: int, n_pixels: int) -> TAScan:
        zeros = lambda: np.zeros((n_times, n_pixels), dtype=float)
        return TAScan(*(zeros() for _ in range(9)))

    def _scan_from_json(self, source: dict[str, Any]) -> TAScan:
        shape = (self.nTimes, self.nPixels)
        scan = self._empty_scan(*shape)
        # These names cover the MATLAB TAScan.populate714 convention and the
        # snake_case fields emitted by newer acquisition software.
        fields = {
            "TAMean": ("TAMean", "ta_mean", "taMean"),
            "TAVariance": ("TAVariance", "ta_variance", "taVariance"),
            "TANShots": ("TANShots", "ta_n_shots", "taNShots", "n_shots"),
            "pumpOnMean": ("pumpOnMean", "pump_on_mean"),
            "pumpOnVariance": ("pumpOnVariance", "pump_on_variance"),
            "pumpOnNShots": ("pumpOnNShots", "pump_on_n_shots"),
            "pumpOffMean": ("pumpOffMean", "pump_off_mean"),
            "pumpOffVariance": ("pumpOffVariance", "pump_off_variance"),
            "pumpOffNShots": ("pumpOffNShots", "pump_off_n_shots"),
        }
        for attribute, spellings in fields.items():
            value = self._field(source, *spellings)
            if value is not None:
                setattr(scan, attribute, self._array(value, shape=shape))

        # Burritos 7.14 files store a scan as one record per time delay:
        # scans -> scan -> spectrum -> transient_absorption/pump_on/pump_off.
        # The original MATLAB implementation delegates this layout to
        # TAScan.populate714; handle it directly here because the Python port
        # has no separate TAScan.m dependency.
        records = source.get("scan")
        if records is not None:
            if not isinstance(records, list) or len(records) != self.nTimes:
                raise ValueError(
                    "Burritos scan records must be a list with one entry for each time delay."
                )
            record_fields = {
                "transient_absorption": ("TAMean", "TAVariance", "TANShots"),
                "pump_on": ("pumpOnMean", "pumpOnVariance", "pumpOnNShots"),
                "pump_off": ("pumpOffMean", "pumpOffVariance", "pumpOffNShots"),
            }
            for row, record in enumerate(records):
                try:
                    spectrum = record["spectrum"]
                except (KeyError, TypeError) as error:
                    raise ValueError(f"Burritos scan record {row} has no spectrum data.") from error
                for source_name, attributes in record_fields.items():
                    try:
                        values = spectrum[source_name]
                        mean = self._array(values["mean"]).reshape(-1)
                        variance = self._array(values["variance"]).reshape(-1)
                        shots = self._array(values["num_shots"]).reshape(-1)
                    except (KeyError, TypeError, ValueError) as error:
                        raise ValueError(
                            f"Burritos scan record {row} has invalid {source_name} statistics."
                        ) from error
                    if any(array.size != self.nPixels for array in (mean, variance, shots)):
                        raise ValueError(
                            f"Burritos scan record {row} has a channel count that does not match num_pixels."
                        )
                    setattr(scan, attributes[0], np.asarray(getattr(scan, attributes[0])))
                    getattr(scan, attributes[0])[row] = mean
                    getattr(scan, attributes[1])[row] = variance
                    getattr(scan, attributes[2])[row] = shots
        return scan

    def loadFile(self) -> None:
        """Load a Burritos JSON acquisition."""
        with Path(self.fileName).open(encoding="utf-8") as handle:
            data = json.load(handle)
        self.json = data

        spectrometer = data["Spectrometer"]
        transient = data["TransientAbsorption"]
        self.nPixels = int(self._field(spectrometer, "num_pixels", "n_pixels"))
        self.buildPixelsVector()
        calibration = spectrometer["wavelength_calibration"]
        self.Calibration = np.array([
            float(calibration["slope"]), float(calibration["intercept"])
        ])
        self.buildWavelengthVector()

        for point in spectrometer.get("calibration_points") or []:
            self.CalibrationPoints.append(CalibrationPoint(
                float(point["wavelength"]), float(point["pixel"]),
                float(point.get("uncertainty", 0.0)),
            ))

        self.pumpWavelength = float(transient.get("pump_wavelength", 0.0))
        self.pumpBandwidth = float(transient.get("pump_bandwidth", 0.0))
        self.pumpEnergy = float(transient.get("pump_energy", 0.0))
        self.solvent = str(transient.get("solvent", ""))
        for attr, key in (("opticalDensity", "optical_density"),
                          ("jetPathLength", "jet_path_length"),
                          ("pumpAngle", "pump_angle"),
                          ("pumpSpotSize", "pump_spot_size"),
                          ("probeSpotSize", "probe_spot_size")):
            setattr(self, attr, float(transient.get(key, 0.0) or 0.0))

        self.times = self._array(transient["time_delays"]).reshape(-1)
        self.nTimes = self.times.size
        scan_data = transient.get("scans", [])
        self.scans = [self._scan_from_json(scan) for scan in scan_data]
        self.nScans = len(self.scans)
        self.combine_scans()

    def loadFile716(self) -> None:
        """Load a Specter 7.16 text export using the layout of the MATLAB loader."""
        lines = Path(self.fileName).read_text(encoding="utf-8", errors="replace").splitlines()
        self.nPixels = 256
        self.buildPixelsVector()
        entries: list[tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []

        def values(index: int) -> np.ndarray:
            result = np.fromstring(lines[index], sep=" ", dtype=float)
            if result.size != self.nPixels:
                raise ValueError(f"Expected {self.nPixels} values near line {index + 1}.")
            return result

        for index, line in enumerate(lines):
            if line.strip() == "# Spectrograph:  Slope and intercept":
                self.Calibration = np.fromstring(lines[index + 1], sep=" ", dtype=float)
                if self.Calibration.size != 2:
                    raise ValueError("Spectrograph calibration must contain slope and intercept.")
                self.buildWavelengthVector()
            elif line.strip() == "# Pump-Probe Delay":
                entries.append((float(lines[index + 1]), values(index + 4), values(index + 7),
                                values(index + 13), 1e-3 * values(index + 19),
                                (1e-3 * values(index + 22)) ** 2))

        if not entries:
            raise ValueError("No '# Pump-Probe Delay' records found in Specter file.")
        if not hasattr(self, "Calibration"):
            raise ValueError("No spectrograph calibration found in Specter file.")

        self.times = np.unique([entry[0] for entry in entries])
        self.nTimes = self.times.size
        self.nScans = int(np.ceil(len(entries) / self.nTimes))
        self.scans = [self._empty_scan(self.nTimes, self.nPixels) for _ in range(self.nScans)]
        for entry_index, (time, nshots, pump_off, pump_on, ta_mean, ta_var) in enumerate(entries):
            scan = self.scans[entry_index // self.nTimes]
            row = int(np.argmin(np.abs(self.times - time)))
            scan.TANShots[row] = nshots
            scan.pumpOffMean[row] = pump_off
            scan.pumpOnMean[row] = pump_on
            scan.TAMean[row] = ta_mean
            scan.TAVariance[row] = ta_var
            # 7.16 exports do not provide separate pump on/off variances.
            scan.pumpOnNShots[row] = nshots
            scan.pumpOffNShots[row] = nshots
        self.combine_scans()

    def buildPixelsVector(self) -> None:
        self.pixels = np.arange(1, self.nPixels + 1, dtype=float)

    def buildWavelengthVector(self) -> None:
        self.wavelengths = self.pixels * self.Calibration[0] + self.Calibration[1]

    def updateCalibration(self, calibrationPolynomial: Sequence[float]) -> None:
        calibration = np.asarray(calibrationPolynomial, dtype=float).reshape(-1)
        if calibration.size != 2:
            raise ValueError("Calibration must contain slope and intercept.")
        self.Calibration = calibration
        self.buildWavelengthVector()

    @staticmethod
    def _weighted_combine(means: np.ndarray, variances: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        valid = np.isfinite(means) & np.isfinite(variances)
        weights = np.where(valid, weights, 0.0)
        total = weights.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.nansum(weights * means, axis=0) / total
            variance = np.nansum(weights * variances, axis=0) / total
        return mean, np.nan_to_num(variance), total

    def combine_scans(self) -> None:
        """Recompute weighted experiment-wide means from all scans."""
        shape = (self.nTimes, self.nPixels)
        if not self.scans:
            self.TAMean = np.full(shape, np.nan)
            self.TAVariance = self.TANShots = np.zeros(shape)
            self.pumpOnMean = self.pumpOnVariance = self.pumpOnNShots = np.zeros(shape)
            self.pumpOffMean = self.pumpOffVariance = self.pumpOffNShots = np.zeros(shape)
            return
        def stack(name: str) -> np.ndarray:
            return np.stack([getattr(scan, name) for scan in self.scans])
        self.TAMean, self.TAVariance, self.TANShots = self._weighted_combine(
            stack("TAMean"), stack("TAVariance"), stack("TANShots"))
        self.pumpOnMean, self.pumpOnVariance, self.pumpOnNShots = self._weighted_combine(
            stack("pumpOnMean"), stack("pumpOnVariance"), stack("pumpOnNShots"))
        self.pumpOffMean, self.pumpOffVariance, self.pumpOffNShots = self._weighted_combine(
            stack("pumpOffMean"), stack("pumpOffVariance"), stack("pumpOffNShots"))

    def _nearest(self, values: np.ndarray, value: float) -> int:
        return int(np.nanargmin(np.abs(values - value)))

    def getTimeTrace(self, value: float | Sequence[float], Units: str = "Wavelengths", dataType: DataType = "TA") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        mean, variance, shots = self._data_arrays(dataType)
        values = np.asarray(value, dtype=float).reshape(-1)
        axis = self.wavelengths if Units.lower() == "wavelengths" else self.pixels
        indices = [self._nearest(axis, item) for item in values]
        if len(indices) == 1:
            index = indices[0]
            return mean[:, index], variance[:, index], shots[:, index]
        if len(indices) != 2:
            raise ValueError("value must be one coordinate or a two-coordinate range.")
        low, high = sorted(indices)
        selected_shots = shots[:, low:high + 1]
        total = selected_shots.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return (np.nansum(selected_shots * mean[:, low:high + 1], axis=1) / total,
                    np.nansum(selected_shots * variance[:, low:high + 1], axis=1) / total,
                    total)

    def getWavelengthTrace(self, time: float | Sequence[float], dataType: DataType = "TA") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        mean, variance, shots = self._data_arrays(dataType)
        values = np.asarray(time, dtype=float).reshape(-1)
        indices = [self._nearest(self.times, item) for item in values]
        if len(indices) == 1:
            index = indices[0]
            return mean[index], variance[index], shots[index]
        if len(indices) != 2:
            raise ValueError("time must be one delay or a two-delay range.")
        low, high = sorted(indices)
        selected_shots = shots[low:high + 1]
        total = selected_shots.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return (np.nansum(selected_shots * mean[low:high + 1], axis=0) / total,
                    np.nansum(selected_shots * variance[low:high + 1], axis=0) / total,
                    total)

    def _data_arrays(self, data_type: DataType) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        data = {"TA": (self.TAMean, self.TAVariance, self.TANShots),
                "Pump On": (self.pumpOnMean, self.pumpOnVariance, self.pumpOnNShots),
                "Pump Off": (self.pumpOffMean, self.pumpOffVariance, self.pumpOffNShots)}
        try:
            return data[data_type]
        except KeyError as error:
            raise ValueError("dataType must be 'TA', 'Pump On', or 'Pump Off'.") from error

    def subtractTABackground(self, upperLimit: float) -> None:
        """Subtract each scan's mean signal at delays up to ``upperLimit``."""
        index = self._nearest(self.times, upperLimit)
        for scan in self.scans:
            background = np.nanmean(scan.TAMean[:index + 1], axis=0)
            scan.TAMean = scan.TAMean - background
        self.combine_scans()

    def buildDispersionFit(self) -> None:
        if self.dispersionFitCoefficients.size == 0:
            raise ValueError("No dispersion fit coefficients have been supplied.")
        self.dispersionFit = np.polynomial.polynomial.polyval(self.pixels, self.dispersionFitCoefficients[::-1])

    def applyExternalDispersionCorrection(self, dispersionFitCoefficients: Sequence[float]) -> None:
        coefficients = np.asarray(dispersionFitCoefficients, dtype=float).reshape(-1)
        if coefficients.size == 0 or not np.all(np.isfinite(coefficients)):
            raise ValueError("Dispersion coefficients must be a non-empty sequence of finite values.")
        self.dispersionFitCoefficients = coefficients
        self.buildDispersionFit()
        self.dispersionCorrectData()

    def correctDispersion(self, coefficients: Sequence[float] | None = None) -> None:
        """Apply a supplied dispersion polynomial.

        MATLAB opens an App Designer dialog to collect its coefficients.  The
        Python model deliberately has no GUI dependency, so callers supply the
        coefficients (highest power first), as they do for
        :meth:`applyExternalDispersionCorrection`.
        """
        if coefficients is None:
            raise ValueError(
                "Python does not include MATLAB's interactive dispersion dialog; "
                "supply polynomial coefficients instead."
            )
        self.applyExternalDispersionCorrection(coefficients)

    @staticmethod
    def fitDispersionCoefficients(pixelPoints: Sequence[float], timePoints: Sequence[float],
                                  polynomialOrder: int = 2) -> np.ndarray:
        """Fit a dispersion polynomial through user-selected map points.

        Coefficients are returned highest-power first, matching ``numpy.polyfit``
        and the ordering used by the original MATLAB ``TAExperiment`` class.
        """
        x = np.asarray(pixelPoints, dtype=float).reshape(-1)
        y = np.asarray(timePoints, dtype=float).reshape(-1)
        try:
            degree = int(polynomialOrder)
        except (TypeError, ValueError) as error:
            raise ValueError("Polynomial order must be an integer.") from error
        if degree < 1:
            raise ValueError("Polynomial order must be at least one.")
        if x.size != y.size:
            raise ValueError("Pixel and time point arrays must have the same length.")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise ValueError("Selected dispersion points must all be finite.")
        if x.size < degree + 1 or np.unique(x).size < degree + 1:
            raise ValueError(
                f"An order-{degree} fit requires at least {degree + 1} distinct pixel positions."
            )
        return np.polyfit(x, y, degree)

    def correctDispersionFromPoints(self, pixelPoints: Sequence[float],
                                    timePoints: Sequence[float],
                                    polynomialOrder: int = 2) -> np.ndarray:
        """Fit clicked time-zero points and apply the resulting correction."""
        coefficients = self.fitDispersionCoefficients(
            pixelPoints, timePoints, polynomialOrder
        )
        self.applyExternalDispersionCorrection(coefficients)
        return coefficients.copy()

    def dispersionCorrectData(self) -> None:
        """Shift each TA pixel trace using MATLAB's linear extrapolation behavior."""
        if self.dispersionFit.shape != self.pixels.shape:
            raise ValueError("Build a dispersion fit before correcting data.")
        corrected = np.empty_like(self.TAMean)
        for column, shift in enumerate(self.dispersionFit):
            shifted_times = self.times - shift
            corrected[:, column] = self._interpolate_extrapolate(
                self.times, shifted_times, self.TAMean[:, column])
        self.TAMean = corrected

    @staticmethod
    def _interpolate_extrapolate(query: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        result = np.interp(query, x, y)
        if x.size > 1:
            left = query < x[0]
            right = query > x[-1]
            result[left] = y[0] + (query[left] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
            result[right] = y[-1] + (query[right] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
        return result

    def pumpPhotonEnergy(self) -> float:
        if self.pumpWavelength <= 0:
            raise ValueError("Pump wavelength value not valid.")
        return (6.636e-34 * 2.998e8) / (self.pumpWavelength * 1e-9)

    def numPumpPhotons(self) -> float:
        if self.pumpEnergy <= 0:
            raise ValueError("Pump energy value not valid.")
        return (self.pumpEnergy * 1e-9) / self.pumpPhotonEnergy()

    def numExcitedMolecules(self) -> float:
        if self.opticalDensity <= 0:
            raise ValueError("Optical density value not valid.")
        photons = self.numPumpPhotons()
        return photons - photons / 10 ** self.opticalDensity

    def pumpArea(self) -> float:
        if self.pumpSpotSize <= 0:
            raise ValueError("Pump spot size value not valid.")
        return np.pi * ((self.pumpSpotSize * 1e-6) / 2) ** 2

    def pumpCylinderVolume(self) -> float:
        if self.jetPathLength <= 0:
            raise ValueError("Jet path length value not valid.")
        return self.pumpArea() * self.jetPathLength * 1e-6 * 1000

    def concExcitedMolecules(self) -> float:
        return (self.numExcitedMolecules() / 6.02e23) / self.pumpCylinderVolume()

    def makeTitle(self) -> str:
        if self.analytes:
            names = []
            for analyte in self.analytes:
                name = getattr(analyte, "name", "Unknown")
                concentration = getattr(analyte, "concentration", np.nan)
                names.append(f"{name} ({concentration:.2f} \\u03bcM)")
            analyte_text = "  ".join(names)
        else:
            analyte_text = "Analyte Unknown"
        return f"{analyte_text}/ {self.solvent}"

    def makeSubtitle(self) -> str:
        return (f"\\u03bb_exc = {self.pumpWavelength:.0f} nm  "
                f"\\u0394\\u03bb = {self.pumpBandwidth:.0f} nm  E = {self.pumpEnergy:.1f} nJ")

    def peekTA(self, Units: str = "Wavelengths", TimeLimits: Sequence[float] | None = None):
        """Return a Matplotlib ``(figure, axes, colorbar)`` TA-map preview."""
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots()
        use_wavelengths = Units.lower() == "wavelengths"
        xaxis = self.wavelengths if use_wavelengths else self.pixels
        image = axes.pcolormesh(xaxis, self.times, self.TAMean, shading="auto", cmap="RdBu_r")
        finite = self.TAMean[np.isfinite(self.TAMean)]
        if finite.size:
            limit = min(abs(finite.min()), abs(finite.max()))
            if limit:
                image.set_clim(-limit, limit)
        axes.set_xlabel("Wavelength / nm" if use_wavelengths else "Pixels")
        axes.set_ylabel("Time / ps")
        axes.set_title(self.makeTitle())
        if TimeLimits is not None:
            axes.set_ylim(TimeLimits)
        colorbar = figure.colorbar(image, ax=axes)
        colorbar.set_label("Transient Absorption / OD")
        return figure, axes, colorbar

    def peekSpectra(self, times: Sequence[float], ShowConfidenceIntervals: bool = False,
                    Units: str = "Wavelengths", dataType: DataType = "TA"):
        """Return a Matplotlib spectral-slice preview."""
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots()
        xaxis = self.wavelengths if Units.lower() == "wavelengths" else self.pixels
        for time in np.asarray(times, dtype=float).reshape(-1):
            mean, variance, _ = self.getWavelengthTrace(time, dataType)
            line = axes.scatter(xaxis, mean, s=10, label=f"{time:.2f} / ps")
            if ShowConfidenceIntervals:
                uncertainty = 2 * np.sqrt(variance)
                axes.fill_between(xaxis, mean - uncertainty, mean + uncertainty,
                                  color=line.get_facecolor()[0], alpha=0.25)
        axes.set_xlabel("Wavelength / nm" if Units.lower() == "wavelengths" else "Pixels")
        axes.set_ylabel("Transient Absorption / OD" if dataType == "TA" else "Intensity")
        axes.set_title(self.makeTitle())
        axes.legend()
        return figure, axes

    def peekTimeTraces(self, xvalues: Sequence[float], ShowConfidenceIntervals: bool = False,
                       Units: str = "Wavelengths", dataType: DataType = "TA"):
        """Return a Matplotlib kinetic-trace preview."""
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots()
        unit = "nm" if Units.lower() == "wavelengths" else "pixels"
        for value in np.asarray(xvalues, dtype=float).reshape(-1):
            mean, variance, _ = self.getTimeTrace(value, Units, dataType)
            line = axes.scatter(self.times, mean, s=10, label=f"{value:.2f} / {unit}")
            if ShowConfidenceIntervals:
                uncertainty = 2 * np.sqrt(variance)
                axes.fill_between(self.times, mean - uncertainty, mean + uncertainty,
                                  color=line.get_facecolor()[0], alpha=0.25)
        axes.set_xlabel("Time / ps")
        axes.set_ylabel("Transient Absorption / OD" if dataType == "TA" else "Intensity")
        axes.set_title(self.makeTitle())
        axes.legend()
        return figure, axes


__all__ = ["TAExperiment", "TAScan", "CalibrationPoint"]
