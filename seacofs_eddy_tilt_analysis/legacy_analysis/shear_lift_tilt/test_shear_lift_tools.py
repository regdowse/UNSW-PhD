import sys
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for folder in (ROOT / "tilt_mechanisms", Path(__file__).resolve().parent):
    sys.path.insert(0, str(folder))

from shear_lift_tools import add_lagged_tilt_change, add_shear_relative_tilt


class ShearLiftGeometryTests(unittest.TestCase):
    def test_eastward_shear_projection(self):
        frame = pd.DataFrame({
            "Eddy": [1], "Day": [0], "Cyc": ["AE"],
            "TiltDir": [0.0], "TiltDis": [10.0],
            "se": [2.0], "sn": [0.0],
        })
        out = add_shear_relative_tilt(frame, "se", "sn", prefix="x")
        self.assertAlmostEqual(out.loc[0, "x_tilt_parallel_km"], 0.0)
        self.assertAlmostEqual(out.loc[0, "x_tilt_left_km"], 10.0)
        self.assertAlmostEqual(out.loc[0, "x_offset_deg"], 90.0)

    def test_northward_shear_projection(self):
        frame = pd.DataFrame({
            "Eddy": [1], "Day": [0], "Cyc": ["CE"],
            "TiltDir": [270.0], "TiltDis": [4.0],
            "se": [0.0], "sn": [3.0],
        })
        out = add_shear_relative_tilt(frame, "se", "sn", prefix="x")
        self.assertAlmostEqual(out.loc[0, "x_tilt_parallel_km"], 0.0, places=12)
        self.assertAlmostEqual(out.loc[0, "x_tilt_left_km"], 4.0)

    def test_exact_lag_and_current_basis(self):
        frame = pd.DataFrame({
            "Eddy": [1, 1, 1], "Day": [0, 1, 3], "Cyc": ["AE"] * 3,
            "TiltDir": [90.0, 90.0, 0.0], "TiltDis": [1.0, 3.0, 3.0],
            "se": [1.0] * 3, "sn": [0.0] * 3,
        })
        relative = add_shear_relative_tilt(frame, "se", "sn", prefix="x")
        out = add_lagged_tilt_change(relative, lag_days=1, prefix="x")
        self.assertAlmostEqual(out.loc[0, "x_1d_dtilt_parallel_km"], 2.0)
        self.assertTrue(np.isnan(out.loc[1, "x_1d_dtilt_parallel_km"]))


if __name__ == "__main__":
    unittest.main()
