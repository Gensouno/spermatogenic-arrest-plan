import sys
from pathlib import Path
import unittest
from unittest.mock import patch
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import run_stage_specificity as module
from audit_stage_sources import ena_records,sra_records
from run_stage_resolution import descriptor_audit


class StageSpecificityGuards(unittest.TestCase):
    def fixture(self):
        target=pd.Series(dict(control_mean_logCPM=2.,Chen_peak_group='Meiotic',Chen_group_dominance_share=.6))
        background=pd.DataFrame(dict(symbol=['b','a','different_stage','different_dominance','far'],control_mean_logCPM=[2.2,2.2,2.1,2.1,4.],Chen_peak_group=['Meiotic','Meiotic','Postmeiotic','Meiotic','Meiotic'],Chen_group_dominance_share=[.6,.6,.6,.95,.6],logFC=[1.,-1.,8.,-8.,0.]))
        return target,background

    def test_neighbors_do_not_select_on_effects(self):
        target,bg=self.fixture();before=module.neighbors(target,bg,'Chen',1,True)
        bg.logFC=[100,-200,300,-400,500]
        self.assertEqual(before,module.neighbors(target,bg,'Chen',1,True))

    def test_stage_dominance_and_abundance_constraints(self):
        target,bg=self.fixture();ids,_=module.neighbors(target,bg,'Chen',1,True)
        self.assertEqual(ids,['a','b'])
        loose,_=module.neighbors(target,bg,'Chen',1,False)
        self.assertIn('different_stage',loose);self.assertNotIn('far',loose)

    def test_mathematical_boundary_is_not_lost_to_float_error(self):
        target,bg=self.fixture();bg=bg.iloc[:1].copy();bg['Chen_group_dominance_share']=.8;bg['control_mean_logCPM']=3.
        self.assertEqual(module.neighbors(target,bg,'Chen',1,True)[0],['b'])

    def test_neighborhood_size_does_not_change_target_weight(self):
        keys=dict(cohort='toy',reference='Chen',abundance_caliper=1.,comparison='abundance_only')
        records=[]
        for target,n,effect in [('x',5,-1.),('y',20,-3.)]:
            records.extend(dict(**keys,target=target,background=f'{target}{i}',target_logFC=0.,background_logFC=effect,abundance_gap=0.,dominance_gap=0.) for i in range(n))
        summary=pd.DataFrame([dict(**keys,matched_priorities=2)])
        _,result=module.summarize(pd.DataFrame(),pd.DataFrame(records),summary)
        self.assertEqual(result.background_mean_logFC.iloc[0],-2.)

    def test_empty_match_is_reported_without_invented_zero_effect(self):
        keys=dict(cohort='toy',reference='Chen',abundance_caliper=1.,comparison='abundance_only')
        empty=pd.DataFrame(columns=list(keys)+['target','target_logFC','background_logFC','abundance_gap','dominance_gap'])
        _,result=module.summarize(pd.DataFrame(),empty,pd.DataFrame([dict(**keys,matched_priorities=0)]))
        self.assertTrue(np.isnan(result.mean_logFC_difference.iloc[0]))

    def test_ena_partial_payload_and_privacy_are_explicit(self):
        data=b'<SAMPLE_SET><SAMPLE accession="SAMN1"><TITLE>a</TITLE><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>age</TAG><VALUE>p14</VALUE></SAMPLE_ATTRIBUTE><SAMPLE_ATTRIBUTE><TAG>biomaterial_provider</TAG><VALUE>private</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE></SAMPLE_SET>'
        result=ena_records(data);self.assertEqual(set(result),{'SAMN1'});self.assertNotIn('biomaterial_provider',result['SAMN1']['attributes'])

    def test_sra_ambiguous_biosample_link_fails(self):
        data=b'<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><SAMPLE><IDENTIFIERS><EXTERNAL_ID namespace="BioSample">SAMN1</EXTERNAL_ID><EXTERNAL_ID namespace="BioSample">SAMN2</EXTERNAL_ID></IDENTIFIERS></SAMPLE><EXPERIMENT accession="SRX1"/></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>'
        with self.assertRaises(ValueError):sra_records(data)

    def test_single_stage_peak_does_not_relabel_group_mean_peak(self):
        stage=pd.DataFrame([dict(symbol='x',Chen_measured=True,Chen_peak_group='Postmeiotic',Chen_peak_stage='MI',Chen_peak_group_mean_CPM=10.,Chen_group_dominance_share=.7,Green_measured=True,Green_peak_group='Postmeiotic',Green_peak_stage='GC9',Green_peak_group_mean_CPM=10.,Green_group_dominance_share=.7)])
        original=stage.copy();design=dict(Chen_groups={'Meiotic':['MI'],'Postmeiotic':['RS2']},Green_groups={'Postmeiotic':['GC9']})
        audit=descriptor_audit(stage,{'x'},design)
        self.assertEqual(audit.query("reference == 'Chen' and subset == 'fixed_priorities'").peak_stage_in_other_mean_peak_group.iloc[0],1)
        pd.testing.assert_frame_equal(stage,original)


if __name__=='__main__':unittest.main()
