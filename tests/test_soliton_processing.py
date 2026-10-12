"""Numerical regressions for both source and packaged Soliton processors."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


class SolitonProcessingTests(unittest.TestCase):
    def models(self):
        for folder in (ROOT / 'imports', ROOT / 'LabInterfaceApp/imports'):
            spec = importlib.util.spec_from_file_location('soliton_model', folder / 'SolitonTAExperimentSelfContained.py')
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            yield module.SolitonTAExperimentSelfContained

    def test_background_and_reapplication(self):
        times = np.array([-3., -2., -1., 0., 1., 2.])
        data = np.array([[2., 4.], [4., 8.], [9., 11.], [12., 14.], [16., 19.], [20., 25.]])
        for model in self.models():
            with self.subTest(model=model), tempfile.TemporaryDirectory() as folder:
                np.savetxt(Path(folder) / 'scan1.txt', np.column_stack((times, data)))
                dat = model(folder, 'scan1.txt', 1)
                dat.mainSoliton_SC([1., 2.], False, -1.85)
                np.testing.assert_allclose(dat.background, [3., 6.])
                np.testing.assert_allclose(dat.backgroundStdev, [1., 2.])
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub_CC, data - [3., 6.])
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub[:2].mean(axis=0), 0.)
                np.testing.assert_array_equal(dat.TAMeanSorted, data)
                dat.chirpCorrectionSoliton_SC([0., 0., 0., 0., 1.])
                dat.backgroundSubtractSoliton_SC(-1.)
                expected = data - data[:3].mean(axis=0)
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub, expected)
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub_CC[:-1], expected[1:])
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub_CC[-1], 0.)
                dat.backgroundSubtractSoliton_SC(-1.)
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub, expected)
                for bad_bound in (np.nan, np.inf, -100., 100.):
                    with self.assertRaises(ValueError):
                        dat.backgroundSubtractSoliton_SC(bad_bound)

    def test_mirrored_multiscan_averaging(self):
        for model in self.models():
            with self.subTest(model=model), tempfile.TemporaryDirectory() as folder:
                times = [-2., -1., 0., 0., -1., -2.]
                np.savetxt(Path(folder) / 'scan1.txt', np.column_stack((times, [2., 4., 8., 10., 6., 4.])))
                np.savetxt(Path(folder) / 'scan2.txt', np.column_stack((times, [4., 6., 10., 12., 8., 6.])))
                dat = model(folder, 'scan1.txt', 2)
                dat.mainSoliton_SC([1.], True, -1.)
                np.testing.assert_allclose(dat.timesSorted, [-2., -1., 0.])
                np.testing.assert_allclose(dat.TAMeanSorted[:, 0], [4., 6., 10.])
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub_CC[:, 0], [-1., 1., 5.])

    def test_repository_example(self):
        scan = next((ROOT / 'TA_andAbsFilesForTommy').rglob('scan1_unref_DOD.txt'))
        for model in self.models():
            with self.subTest(model=model):
                raw = np.genfromtxt(scan)
                pixels = np.arange(1, raw.shape[1], dtype=float)
                squared_energy = 12.4698**2 - 2.7639e5 / (pixels + 1.3204e3)
                energy = np.sqrt(np.maximum(squared_energy, 0.))
                dat = model(scan.parent, scan.name, 1)
                dat.mainSoliton_SC(energy, True, -1.85)
                count = dat.backgroundSubtractionUpperBoundInxed + 1
                np.testing.assert_allclose(dat.background, dat.TAMeanSorted[:count].mean(axis=0))
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub_CC, dat.TAMeanSorted - dat.background)
                np.testing.assert_allclose(dat.TAMeanSortedBackgroundSub[:count].mean(axis=0), 0., atol=1e-14)


if __name__ == '__main__':
    unittest.main()
