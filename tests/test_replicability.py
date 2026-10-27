"""Replication-family arithmetic, age intersections and stage-context safeguards."""
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts')]
from replicability_statistics import intersection_union,partial_conjunction,adjust_joint_family,directional_tail_from_two_sided
from analyze_replicability import make_unit_evidence,PRIMARY,UNITS,loo_summaries
from annotate_replicability import summarize_context
from unittest.mock import patch
from check_replicability_context import unique_effects,transfer_summary


class ReplicabilityTests(unittest.TestCase):
    def test_one_small_study_cannot_establish_three_units(self):
        self.assertEqual(partial_conjunction([1e-20,.4,.7,.8],3),1)
        np.testing.assert_allclose(partial_conjunction([[.01,.02,.03,.9],[.1,.2,.3,.4]],3),[.06,.6])
        self.assertAlmostEqual(partial_conjunction([.01,.02,.03,.9],1),.04)
        self.assertAlmostEqual(partial_conjunction([.01,.02,.03,.9],4),.9)

    def test_partial_conjunction_axis_and_invalid_values(self):
        np.testing.assert_allclose(partial_conjunction([[.01,.1],[.02,.2]],2,axis=0),[.02,.2])
        for p,r in [([.1],0),([.1],2),([np.nan],1),([-.1],1),([1.1],1),([],1)]:
            with self.assertRaises(ValueError):partial_conjunction(p,r)

    def test_intersection_is_worst_age(self):
        np.testing.assert_allclose(intersection_union([[.001,.8],[.9,.002]]),[.8,.9])
        with self.assertRaises(ValueError):intersection_union([])

    def test_joint_family_by_known_fixture_and_no_rowwise_adjustment(self):
        p=np.array([[.001,.02],[.1,.8]])
        np.testing.assert_allclose(adjust_joint_family(p),np.array([[1/120,1/12],[5/18,1]]))
        self.assertFalse(np.allclose(adjust_joint_family(p),np.stack([adjust_joint_family(r) for r in p])))
        self.assertTrue((adjust_joint_family(p)>=adjust_joint_family(p,'bh')).all())

    def test_tails_preserve_direction_and_zero(self):
        p=np.array([.02,.04,1]);t=np.array([2,-2,0])
        np.testing.assert_allclose(directional_tail_from_two_sided(p,t,'up'),[.01,.98,.5])
        np.testing.assert_allclose(directional_tail_from_two_sided(p,t,'down'),[.99,.02,.5])
        with self.assertRaises(ValueError):directional_tail_from_two_sided([.5],[0],'up')

    def test_rfx2_is_one_unit_and_opposite_ages_veto_common_direction(self):
        rows=[]
        for g in PRIMARY:
            negative=g=='GSE68284'
            rows.append({'gse':g,'gene':'A','logFC':-1 if negative else 1,
                         'P_up':.99 if negative else .01,'P_down':.01 if negative else .99})
        model,p,pos,neg=make_unit_evidence(pd.DataFrame(rows),['A'])
        self.assertEqual(len(UNITS),4);self.assertEqual(p['up'].shape,(1,4))
        self.assertEqual(model.loc['A','Rfx2_direction'],'mixed_or_zero')
        self.assertEqual(model.loc['A','Rfx2_P_up'],.99)
        self.assertEqual(model.loc['A','Rfx2_P_down'],.99)
        self.assertEqual(model.loc['A','Rfx2_logFC_min'],-1)
        self.assertEqual(model.loc['A','Rfx2_logFC_max'],1)

    def test_context_equal_stage_means_not_sums(self):
        profile=pd.DataFrame([[9,8,8,8]],index=['A'],columns=['a','b','c','d'])
        d=summarize_context(profile,['A'],{'one':['a'],'three':['b','c','d']},'x')
        self.assertEqual(d.loc['A','x_peak_group'],'one')
        self.assertAlmostEqual(d.loc['A','x_group_dominance_share'],9/17)
        self.assertFalse(d.loc['A','x_high_dominance_descriptive'])

    def test_context_missing_zero_and_tie_not_arbitrary_argmax(self):
        profile=pd.DataFrame([[0,0],[20,20],[20,0]],index=['zero','tie','unique'],columns=['a','b'])
        d=summarize_context(profile,['zero','tie','unique','absent'],{'first':['a'],'second':['b']},'x')
        self.assertEqual(d.loc['zero','x_peak_group'],'unexpressed_in_source')
        self.assertEqual(d.loc['tie','x_peak_group'],'tied_groups')
        self.assertEqual(d.loc['absent','x_peak_group'],'not_measured')
        self.assertTrue(d.loc['unique','x_high_dominance_descriptive'])
        self.assertFalse(d.loc['tie','x_high_dominance_descriptive'])
        with self.assertRaises(ValueError):summarize_context(profile,['zero'],{'x':['a','b'],'y':['b']},'x')

    def test_missing_loo_variant_does_not_pass_stability(self):
        # Forty variants exist overall, but A is re-filtered out of one age variant.
        sizes=dict(zip(PRIMARY,[6,8,9,11,6]));rows=[]
        for g,n in sizes.items():
            for i in range(n):
                rows.append({'gse':g,'variant':f'drop_{i}','gene':'B','logFC':-1})
                if not (g=='GSE68284' and i==0):rows.append({'gse':g,'variant':f'drop_{i}','gene':'A','logFC':-1})
        table=pd.DataFrame(rows)
        def reader(path,*args,**kwargs):
            if str(path).endswith('effects.csv.gz'):return table
            g=Path(path).name.split('_')[0]
            return pd.DataFrame({'matrix_sample':range(sizes[g])})
        with patch('analyze_replicability.pd.read_csv',side_effect=reader):
            cohort,units,counts=loo_summaries(ROOT,['A','B'])
        r=units[(units.gene=='A')&(units.unit=='Rfx2')].iloc[0]
        self.assertFalse(r.all_variants_available);self.assertFalse(r.stable_down_1p5fold)
        self.assertEqual(counts.loc['A','LOO_stable_down_1p5fold_units'],3)
        self.assertEqual(counts.loc['B','LOO_stable_down_1p5fold_units'],4)

    def test_external_duplicate_symbols_all_excluded(self):
        table=pd.DataFrame({'symbol':['A','A','B',None],'logFC':[9,-1,.8,-2]})
        result=unique_effects(table,'symbol')
        self.assertEqual(list(result.index),['B'])
        self.assertEqual(result['B'],.8)

    def test_external_missing_never_zero_and_requires_full_fit_and_every_omission(self):
        effect=pd.Series({'A':-1,'B':1,'C':-.2})
        omissions={'drop1':pd.Series({'A':-1,'B':-1,'C':-1}),
                   'drop2':pd.Series({'B':-1,'C':-1})}
        d=transfer_summary(effect,omissions,['A','B','C','D'],[-1,-1,-1,-1]).set_index('gene')
        self.assertTrue(d.loc['A','directional_1p5fold'])
        self.assertFalse(d.loc['A','baseline_and_all_LOO_1p5fold'])
        self.assertFalse(d.loc['B','baseline_and_all_LOO_sign_match'])
        self.assertFalse(d.loc['C','baseline_and_all_LOO_1p5fold'])
        self.assertTrue(np.isnan(d.loc['D','external_logFC']))
        self.assertFalse(d.loc['D','tested_exact_unique_symbol'])


if __name__=='__main__':unittest.main()
