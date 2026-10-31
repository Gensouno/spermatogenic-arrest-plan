import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from check_competitive_set import adjust_family,validate_results
from run_independent_contexts import COHORTS


class CompetitiveFamilyGuards(unittest.TestCase):
    def fixture(self):
        rows=[]
        for i,c in enumerate(COHORTS):
            for policy in ['estimated','fixed_0p01']:
                rho=.5 if policy=='estimated' else .01
                rows.append(dict(cohort=c,policy=policy,PValue=.01*(i+1) if policy=='estimated' else .001*(i+1),
                    fixed_priorities=36,NGenes=30,universe_genes=9000,background_genes=8970,
                    correlation_raw_or_assumed=rho,variance_inflation_used=1+29*rho,
                    camera_df=4 if policy=='estimated' else 8998,Direction='Down',reverse_Direction='Up',
                    priority_mean_logFC=-1.,background_mean_logFC=0.))
        d=pd.DataFrame(rows);d['reverse_PValue']=d.PValue;d['BH'],d['BY']=adjust_family(d)
        return d

    def test_adjustment_known_joint_family(self):
        d=self.fixture();bh,by=adjust_family(d)
        self.assertTrue(np.allclose(bh[[1,3,5,7]],.008))
        self.assertTrue(np.allclose(by,bh*sum(1/i for i in range(1,9))))
        validate_results(d)

    def test_missing_sensitivity_cannot_shrink_family(self):
        with self.assertRaisesRegex(ValueError,'eight-test'):adjust_family(self.fixture().iloc[:-1])

    def test_separate_policy_correction_rejected(self):
        d=self.fixture();d.loc[d.policy.eq('fixed_0p01'),'BH']=.004
        with self.assertRaisesRegex(ValueError,'Joint-family'):validate_results(d)

    def test_reversed_contrast_cannot_keep_direction(self):
        d=self.fixture();d.loc[0,'reverse_Direction']='Down'
        with self.assertRaisesRegex(ValueError,'reversal'):validate_results(d)

    def test_estimated_policy_cannot_inherit_large_df(self):
        d=self.fixture();d.loc[0,'camera_df']=8998
        with self.assertRaisesRegex(ValueError,'degrees'):validate_results(d)

    def test_different_measured_genes_between_policies_rejected(self):
        d=self.fixture();d.loc[0,'NGenes']=29;d.loc[0,'background_genes']=8971;d.loc[0,'variance_inflation_used']=15
        with self.assertRaisesRegex(ValueError,'Policy-dependent'):validate_results(d)

    def test_nonfinite_probabilities_rejected(self):
        d=self.fixture();d.loc[0,'PValue']=np.nan
        with self.assertRaisesRegex(ValueError,'Invalid P'):adjust_family(d)


if __name__=='__main__':unittest.main()
