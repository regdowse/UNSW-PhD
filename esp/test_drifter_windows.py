import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import drifter_windows as w
import platform_weighting as pw


class WindowTests(unittest.TestCase):
    def test_common_period_uses_total_duration(self):
        d=pd.DataFrame({'time':pd.to_datetime(['2023-10-01 06:00','2023-10-10 06:00'])})
        s={'daily':pd.date_range('2023-10-01 06:00',periods=10)}
        t=w.common_times(d,s,[12,48,96])
        self.assertEqual(t[0],pd.Timestamp('2023-10-03 06:00'))
        self.assertEqual(t[-1],pd.Timestamp('2023-10-08 06:00'))
        with self.assertRaises(ValueError):w.common_times(d,s,[0])

    def test_no_restart_and_full_denominator(self):
        times=pd.date_range('2023-10-03 06:00',periods=3)
        data=pd.DataFrame(dict(x=[0.,1000.,0.],y=[0.,0.,1000.],u=0.,v=0.,
                               platform='Drifter',drifter_id='1'))
        with patch.object(w.dw,'snapshot',return_value=data),patch.object(w,'estimate',side_effect=pw.FitError('test_failure')) as fitter:
            result=w.run(None,None,None,[48],times)
        self.assertEqual(fitter.call_count,3)  # one first-day attempt per configuration
        self.assertEqual((result.status=='not_attempted_after_failure').sum(),6)
        c=w.counts(result)
        self.assertTrue((c.days_with_result==0).all())
        self.assertTrue((c.comparison_days==3).all())

    def test_matched_dates_only(self):
        rows=[]
        for hours in [24,48]:
            for day in range(3):
                rows.append(dict(configuration='Drifter only',window_hours=hours,time=day,
                    valid=hours==48 or day==0,xc_km=3 if hours==24 else 0,yc_km=4 if hours==24 else 0,
                    alpha=1.2 if hours==24 else 1.,inner_vorticity=-2e-5 if hours==24 else -1e-5))
        d=w.differences(pd.DataFrame(rows)).set_index('window_hours')
        self.assertEqual(d.loc[24,'matched_days'],1)
        self.assertEqual(d.loc[24,'median_centre_difference_km'],5)
        self.assertAlmostEqual(d.loc[24,'median_abs_vorticity_difference_1e5'],1)
        self.assertEqual(d.loc[48,'max_centre_difference_km'],0)

    def test_unavailable_reference_remains_missing(self):
        r=pd.DataFrame([dict(configuration='Drifter only',window_hours=h,time=1,valid=h==24,
                            xc_km=0,yc_km=0,alpha=1,inner_vorticity=-1e-5) for h in [24,48]])
        d=w.differences(r)
        self.assertTrue((d.matched_days==0).all())
        self.assertTrue(d.median_centre_difference_km.isna().all())

    def test_rank_deficiency_rejected(self):
        d=pd.DataFrame(dict(x=[0.,1.,2.],y=0.,u=0.,v=0.,platform='Drifter'))
        with self.assertRaises(pw.FitError):w.estimate(None,d,False,.1)

if __name__=='__main__':unittest.main()
