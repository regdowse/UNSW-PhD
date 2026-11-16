"""Synthetic checks for scientific conventions, sampling and weighting."""
import unittest

import numpy as np
import pandas as pd

import ellipse_tilt_tools as et


class GeometryTests(unittest.TestCase):
    def test_axes_follow_pipeline_bearing_and_eigenvector_sign(self):
        for alpha in [0, 20, 45, 90, 179, 230]:
            a = np.deg2rad(alpha)
            major = np.array([np.cos(a), np.sin(a)])
            minor = np.array([-np.sin(a), np.cos(a)])
            q = np.outer(major, major) / 9 + np.outer(minor, minor)
            frame = pd.DataFrame([[q[0, 0], q[0, 1], q[1, 1]]], columns=et.QCOLS)
            out = et.ellipse_geometry(frame)
            expected = (np.degrees(np.arctan2(major[0], major[1])) + 20) % 180
            self.assertAlmostEqual(out.AxisRatio.iloc[0], 3)
            self.assertAlmostEqual(float(et.axial_difference(out.MajorBearing, expected)[0]), 0)
            self.assertAlmostEqual(float(et.axial_difference(expected + 180, expected)), 0)

    def test_invalid_and_circular_shapes(self):
        frame = pd.DataFrame([[1, 0, 1], [-1, 0, 1], [1, 2, 1], [np.nan, 0, 1]], columns=et.QCOLS)
        out = et.ellipse_geometry(frame)
        self.assertTrue(out.MajorBearing.isna().all())
        self.assertEqual(out.Q_valid.tolist(), [True, False, False, False])
        self.assertEqual(out.AxisRatio.iloc[0], 1)

    def test_depth_interpolation_and_no_extrapolation(self):
        p = pd.DataFrame({'Eddy': [1]*3, 'Day': [1]*3, 'Depth': [0, 100, 200],
                          'q11': [.25, .5, 1], 'q12': [0]*3, 'q22': [1]*3})
        out = et.sample_depth_geometry(p, [50, 100, 250])
        self.assertEqual(out.ShapeDepth.tolist(), [50, 100])
        self.assertAlmostEqual(out.q11.iloc[0], .375)
        self.assertEqual(out.DepthMethod.tolist(), ['interpolate', 'exact'])
        self.assertTrue(et.sample_depth_geometry(p, [50], max_gap_m=40).empty)
        self.assertTrue(et.sample_depth_geometry(p, [50], method='nearest', nearest_tolerance_m=20).empty)
        p.loc[1, 'q11'] = -1
        self.assertTrue(et.sample_depth_geometry(p, [50, 100, 150]).empty)
        with self.assertRaises(ValueError):
            et.sample_depth_geometry(pd.concat([p, p.iloc[[0]]]))

    def test_matched_keys_are_day_specific(self):
        g = pd.DataFrame({'Eddy':[1,1,1,1], 'Day':[1,1,2,3], 'ShapeDepth':[0,50,0,50]})
        matched = et.matched_depth_sample(g, [0,50])
        self.assertEqual(matched.Day.tolist(), [1,1])

    def test_rotation_does_not_bridge_excluded_days(self):
        g = pd.DataFrame({'ShapeDepth':[0]*4, 'Eddy':[1]*4, 'Day':[1,2,4,5],
                          'MajorBearing':[179,1,20,30], 'TiltDir':[179,1,20,30],
                          'AxisRatio':[2]*4, 'TiltDis':[10]*4})
        out = et.rotation_pairs(g)
        self.assertEqual(out.Day.tolist(), [2,5])
        np.testing.assert_allclose(out.AxisTurn, [2,10])
        np.testing.assert_allclose(out.TurnAgreement, 1)

    def test_within_slope_and_track_duplication_invariance(self):
        frames=[]
        for eddy in range(12):
            ar=np.linspace(1.2, 2.2, 8)
            frames.append(pd.DataFrame({'Eddy':eddy, 'Day':np.arange(8), 'Cyc':'AE',
                          'ShapeDepth':0, 'AxisRatio':ar, 'TiltDis':3*ar + 2*eddy}))
        g=pd.concat(frames, ignore_index=True)
        a=et.magnitude_summary(g, n_boot=30)
        b=et.magnitude_summary(pd.concat([g, g[g.Eddy.eq(0)]], ignore_index=True), n_boot=30)
        slope=a.loc[a.metric.str.startswith('Within'), 'estimate'].iloc[0]
        self.assertAlmostEqual(slope, 3)
        np.testing.assert_allclose(a[['estimate','low','high']], b[['estimate','low','high']])

    def test_empty_summaries_and_surface_only(self):
        surface=pd.DataFrame({'Eddy':[1], 'Day':[1], 'Cyc':['AE'], 'TiltDis':[10],
                              'TiltDir':[110], 'q11':[.25], 'q12':[0], 'q22':[1]})
        sampled=et.sample_depth_geometry(pd.DataFrame(columns=et.KEY+['Depth']+et.QCOLS))
        all_shapes=et.build_analysis_table(surface, sampled)
        self.assertAlmostEqual(all_shapes.AlignmentDeg.iloc[0], 0)
        empty=all_shapes.iloc[:0]
        self.assertTrue(et.magnitude_summary(empty).empty)
        self.assertTrue(et.alignment_summary(empty).empty)


if __name__ == '__main__':
    unittest.main()
