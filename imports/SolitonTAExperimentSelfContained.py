"""Python port of ``SolitonTAExperimentSelfContained.m``.

The public field names intentionally match the MATLAB class, which lets code
written for the Lab Interface use the processed TA matrix, delay axis, and
energy axis without a MATLAB Engine session.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np


class SolitonTAExperimentSelfContained:
    """Load, combine, background-subtract, and chirp-correct Soliton TA scans."""

    def __init__(self, folderPath: str | Path, fileName: str, nScans: int) -> None:
        self.folderPath = Path(folderPath)
        self.fileName = fileName
        self.nScans = int(nScans)
        if self.nScans < 1:
            raise ValueError("nScans must be at least one.")

        self.energyAxis = np.array([], dtype=float)
        self.timesRaw = self.timesSorted = self.timesSorted_CC = np.array([], dtype=float)
        self.TACubeRaw = self.TACubeSorted = np.empty((0, 0, 0), dtype=float)
        self.TAMeanRaw = self.TAMeanSorted = np.empty((0, 0), dtype=float)
        self.TAMeanSortedBackgroundSub = self.TAMeanSortedBackgroundSub_CC = np.empty((0, 0), dtype=float)
        self.energyAxis_CC = np.array([], dtype=float)
        self.background = self.backgroundStdev = np.array([], dtype=float)
        self.backgroundSubtractionUpperBound = np.nan
        self.backgroundSubtractionUpperBoundInxed = -1
        self.chirpCorrectionCoefficients = np.zeros(5, dtype=float)
        self.nTimes = 0

        # Experiment metadata retained from the MATLAB class.
        self.pumpWavelength = self.pumpBandwidth = self.pumpEnergy = 0.0
        self.opticalDensity = self.jetPathLength = self.pumpAngle = 0.0
        self.pumpSpotSize = self.probeSpotSize = 0.0
        self.pumpPolarization = self.probePolarization = 0.0
        self.analytes: list[object] = []
        self.solvent = ""
        self.mainGraphTitle = ""
        self.UVVis_unirradiated: dict[str, object] = {}
        self.UVVis_irradiated: dict[str, object] = {}
        self.pumpSpectrum: dict[str, object] = {}
        self.fluorescenceSpectrum: dict[str, object] = {}

    def _first_scan_path(self) -> Path:
        path = self.folderPath / self.fileName
        if not path.is_file():
            raise FileNotFoundError(f"Soliton scan file not found: {path}")
        return path

    def _scan_path(self, scan_number: int) -> Path:
        """Mirror MATLAB's scan1-to-scanN filename substitution."""
        name = self.fileName.replace("scan1", f"scan{scan_number}")
        if name == self.fileName and scan_number != 1:
            raise ValueError("fileName must contain 'scan1' when nScans is greater than one.")
        return self.folderPath / name

    @staticmethod
    def _read_matrix(path: Path) -> np.ndarray:
        data = np.genfromtxt(path, dtype=float)
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(f"Expected a numeric time-plus-pixels matrix in {path}.")
        # MATLAB replaces both NaN and Inf with zero before any processing.
        return np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

    def loadSolitonScans_SC(self) -> None:
        """Load ``scan1`` through ``scanN`` and make the raw scan cube."""
        first = self._read_matrix(self._first_scan_path())
        cube = np.empty((*first.shape, self.nScans), dtype=float)
        cube[:, :, 0] = first
        for index in range(1, self.nScans):
            path = self._scan_path(index + 1)
            if not path.is_file():
                raise FileNotFoundError(f"Soliton scan file not found: {path}")
            matrix = self._read_matrix(path)
            if matrix.shape != first.shape:
                raise ValueError(f"{path.name} has shape {matrix.shape}; expected {first.shape}.")
            cube[:, :, index] = matrix

        self.timesRaw = cube[:, 0, :].mean(axis=1)
        self.TACubeRaw = cube[:, 1:, :]
        self.TAMeanRaw = self.TACubeRaw.mean(axis=2)
        self.nTimes = self.timesRaw.size

    @staticmethod
    def _sort_unique(times: np.ndarray, cube: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Implement MATLAB ``[B,I] = unique(times)`` for rows of a TA cube."""
        sorted_times, first_indices = np.unique(times, return_index=True)
        return sorted_times, cube[first_indices, :, :]

    def unmirrorScans_SC(self) -> None:
        """Split forward/reverse delay sweeps and average them as independent scans."""
        rows, pixels, scans = self.TACubeRaw.shape
        if rows % 2:
            raise ValueError("Mirrored Soliton scans must have an even number of time rows.")
        half = rows // 2
        unmirrored = np.empty((half, pixels, scans * 2), dtype=float)
        for index in range(scans):
            unmirrored[:, :, 2 * index] = self.TACubeRaw[:half, :, index]
            unmirrored[:, :, 2 * index + 1] = self.TACubeRaw[half:, :, index][::-1]
        self.timesSorted, self.TACubeSorted = self._sort_unique(self.timesRaw[:half], unmirrored)
        self.TAMeanSorted = self.TACubeSorted.mean(axis=2)

    def processSolitonTimes_SC(self) -> None:
        """Sort non-mirrored scans and drop duplicate delay rows."""
        self.timesSorted, self.TACubeSorted = self._sort_unique(self.timesRaw, self.TACubeRaw)
        self.TAMeanSorted = self.TACubeSorted.mean(axis=2)

    def mainSoliton_SC(self, energyAxis: Sequence[float], isMirrored: bool,
                       upperBoundForBS: float,
                       chirpCorrectionCoefficients: Sequence[float] | None = None) -> None:
        """Run the complete MATLAB ``mainSoliton_SC`` processing sequence."""
        self.loadSolitonScans_SC()
        if isMirrored:
            self.unmirrorScans_SC()
        else:
            self.processSolitonTimes_SC()
        axis = np.asarray(energyAxis, dtype=float).reshape(-1)
        if axis.size != self.TAMeanSorted.shape[1]:
            raise ValueError(
                f"energyAxis has {axis.size} values, but the scan has "
                f"{self.TAMeanSorted.shape[1]} data channels."
            )
        # Keep the processing axis non-negative.  Invalid calibration channels
        # remain represented as zero so matrix/channel alignment is preserved;
        # the GUI masks those channels from physical-axis plots.
        valid_energy = np.isfinite(axis) & (axis > 0)
        if not np.any(valid_energy):
            raise ValueError("energyAxis must contain at least one positive finite value.")
        self.energyAxis = np.where(valid_energy, axis, 0.0)
        self.backgroundSubtractSoliton_SC(upperBoundForBS)
        coefficients = (np.zeros(5) if chirpCorrectionCoefficients is None
                        else chirpCorrectionCoefficients)
        self.chirpCorrectionSoliton_SC(coefficients)

    def findNearestIV_SC(self, myArray: Sequence[float], myValue: float) -> tuple[int, float]:
        """Return the zero-based nearest index and its actual value."""
        array = np.asarray(myArray, dtype=float).reshape(-1)
        if not array.size:
            raise ValueError("Cannot search an empty axis.")
        index = int(np.nanargmin(np.abs(array - myValue)))
        return index, float(array[index])

    def backgroundSubtractSoliton_SC(self, upperBoundTime: float) -> None:
        index, actual_time = self.findNearestIV_SC(self.timesSorted, upperBoundTime)
        # MATLAB's 1:UBT_idx includes the nearest row.
        source = self.TAMeanSorted[:index + 1]
        self.background = source.mean(axis=0)
        self.backgroundStdev = source.std(axis=0, ddof=1) if source.shape[0] > 1 else np.zeros(source.shape[1])
        self.TAMeanSortedBackgroundSub = self.TAMeanSorted - self.background
        self.backgroundSubtractionUpperBound = actual_time
        self.backgroundSubtractionUpperBoundInxed = index

    def chirpCorrectionSoliton_SC(self, CCcoeffs: Sequence[float], showFigure: bool = False) -> None:
        """Apply the five-term, channel-index chirp polynomial.

        As in MATLAB, values shifted outside the measured delay range are zero.
        """
        coefficients = np.asarray(CCcoeffs, dtype=float).reshape(-1)
        if coefficients.size != 5 or not np.all(np.isfinite(coefficients)):
            raise ValueError("CCcoeffs must contain exactly five finite values: a4, a3, a2, a1, a0.")
        channels = np.arange(1, self.energyAxis.size + 1, dtype=float)
        shifts = np.polyval(coefficients, channels)
        corrected = np.empty_like(self.TAMeanSortedBackgroundSub)
        for column, shift in enumerate(shifts):
            shifted_times = self.timesSorted - shift
            corrected[:, column] = np.interp(
                self.timesSorted, shifted_times,
                self.TAMeanSortedBackgroundSub[:, column], left=0.0, right=0.0,
            )
        self.TAMeanSortedBackgroundSub_CC = corrected
        self.timesSorted_CC = self.timesSorted.copy()
        self.energyAxis_CC = self.energyAxis.copy()
        self.chirpCorrectionCoefficients = coefficients
        if showFigure:
            self.surfplotSorted_BS_CC_SC()

    def getBinnedTimeTraceForFitting_SC(self, wavelengthbounds: Sequence[float]) -> np.ndarray:
        bounds = np.asarray(wavelengthbounds, dtype=float).reshape(-1)
        if bounds.size != 2:
            raise ValueError("wavelengthbounds must contain lower and upper energy bounds.")
        first, _ = self.findNearestIV_SC(self.energyAxis, bounds[0])
        second, _ = self.findNearestIV_SC(self.energyAxis, bounds[1])
        low, high = sorted((first, second))
        return self.TAMeanSortedBackgroundSub_CC[:, low:high + 1].mean(axis=1)

    def importAbsSpectrum(self, filePath: str | Path, isIrradiated: bool,
                          notes: str = "", containsVariance: bool = True) -> None:
        data = np.genfromtxt(filePath, delimiter=",", dtype=float)
        if data.ndim != 2 or data.shape[1] < 2:
            # Whitespace-delimited exports are also accepted.
            data = np.genfromtxt(filePath, dtype=float)
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError("Absorbance file must contain wavelength and absorbance columns.")
        spectrum: dict[str, object] = {
            "filePath": str(filePath), "notes": notes, "wavelengths": data[:, 0],
            "eV": 1240.0 / data[:, 0], "absorbance": data[:, 1],
        }
        if containsVariance and data.shape[1] >= 3:
            spectrum["variance"] = data[:, 2]
        if isIrradiated:
            self.UVVis_irradiated = spectrum
        else:
            self.UVVis_unirradiated = spectrum

    def _import_intensity_spectrum(self, filePath: str | Path) -> dict[str, object]:
        data = np.genfromtxt(filePath, delimiter=",", dtype=float)
        if data.ndim != 2 or data.shape[1] < 2:
            data = np.genfromtxt(filePath, dtype=float)
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError("Spectrum file must contain wavelength and intensity columns.")
        return {"filePath": str(filePath), "wavelength": data[:, 0],
                "eV": 1240.0 / data[:, 0], "intensity": data[:, 1]}

    def importPumpSpectrum(self, filePath: str | Path) -> None:
        self.pumpSpectrum = self._import_intensity_spectrum(filePath)

    def importFluorescenceSpectrum(self, filePath: str | Path) -> None:
        self.fluorescenceSpectrum = self._import_intensity_spectrum(filePath)

    def saveMyTAForGlotaran_SC(self, toPath: str | Path | None = None,
                                asFileName: str = "TAMeans_BS_for_PyGloTarAn",
                                saveChirpCorrected: bool = False) -> np.ndarray:
        """Export the Glotaran matrix and return the same matrix in memory."""
        destination = Path(toPath) if toPath is not None else self.folderPath
        data = self.TAMeanSortedBackgroundSub_CC.T if saveChirpCorrected else self.TAMeanSortedBackgroundSub.T
        name = f"{asFileName}{'_CC' if saveChirpCorrected else ''}"
        matrix = np.empty((data.shape[0] + 1, data.shape[1] + 1), dtype=float)
        matrix[0, 0] = np.nan
        matrix[0, 1:] = self.timesSorted
        matrix[1:, 0] = self.energyAxis
        matrix[1:, 1:] = data
        np.savetxt(destination / name, matrix, delimiter=",", fmt="%.18g")
        return matrix

    def makeMainGraphTitle_SC(self, Sample: str, Solvent: str, Concentration: float,
                              pumpEnergy: float | None = None, pumpWavelength: float | None = None) -> None:
        energy = self.pumpEnergy if pumpEnergy is None else pumpEnergy
        wavelength = self.pumpWavelength if pumpWavelength is None else pumpWavelength
        self.mainGraphTitle = f"{Sample} {Concentration} mM, in {Solvent} excited by {energy} nJ {wavelength}"

    def surfplotRaw(self):
        return self._surface_plot(self.energyAxis, self.timesRaw, self.TAMeanRaw, "Raw Soliton TA")

    def surfplotRawSingleScan(self, scanNum: int):
        return self._surface_plot(self.energyAxis, self.timesRaw, self.TACubeRaw[:, :, scanNum - 1], "Raw Soliton TA")

    def surfplotSorted(self):
        return self._surface_plot(self.energyAxis, self.timesSorted, self.TAMeanSorted, "Sorted Soliton TA")

    def surfplotSorted_BS_CC_SC(self):
        return self._surface_plot(self.energyAxis_CC, self.timesSorted_CC,
                                  self.TAMeanSortedBackgroundSub_CC, "Background-subtracted, chirp-corrected Soliton TA")

    @staticmethod
    def _surface_plot(xaxis: np.ndarray, yaxis: np.ndarray, data: np.ndarray, title: str):
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots()
        image = axes.pcolormesh(xaxis, yaxis, data, shading="auto", cmap="RdBu_r", vmin=-0.005, vmax=0.005)
        axes.set(xlabel="Probe Energy / eV", ylabel="Time Delay / ps", title=title)
        figure.colorbar(image, ax=axes)
        return figure, axes


__all__ = ["SolitonTAExperimentSelfContained"]
