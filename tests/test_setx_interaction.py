"""Setx interaction direction, design aliasing, missingness and multiplicity."""
from pathlib import Path
import subprocess
import sys
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_setx_interaction import CELLS, audit_design, parse_member, summarize_candidates


class InteractionTests(unittest.TestCase):
    def test_filename_tokens_are_not_inferred_animals(self):
        r=parse_member('GSM7034763_d16_7844_WT_GES14_00144_AGTTCC_GRCm38.trimmed.genes.results.gz')
        self.assertEqual((r['source_code'],r['group'],r['age']),('7844','Control','P16'))
        self.assertNotIn('animal',r)
        for name in ['../GSM1_x','GSM1_d16_1_WT_GES14_1_AGTTCC_mm10_piRNA.trimmed.genes.results.gz']:
            with self.assertRaises(ValueError):parse_member(name)

    def design(self):
        return pd.DataFrame([dict(gsm=f'GSM{i*4+j}',age=c.split('_')[0],group=c.split('_')[1],
            source_code=str(i*4+j),barcode_token=f'index{i}_{j//2}') for i,c in enumerate(CELLS) for j in range(4)])

    def test_nested_barcode_effects_are_not_identifiable(self):
        d=self.design();a=audit_design(d)
        self.assertTrue(a['barcodes_nested_in_cells'])
        self.assertEqual((a['cell_design_rank'],a['cell_plus_barcode_rank'],a['cell_plus_barcode_columns']),(10,20,30))

    def test_repeated_code_retained_without_relabeling(self):
        d=self.design();d.loc[[24,28],'source_code']='same'
        original=d.copy(deep=True);a=audit_design(d)
        self.assertEqual(len(a['repeated_source_code_rows']),2)
        pd.testing.assert_frame_equal(d,original)

    def test_duplicate_or_unbalanced_design_fails(self):
        d=self.design();d.loc[0,'gsm']=d.loc[1,'gsm']
        with self.assertRaises(ValueError):audit_design(d)
        d=self.design();d.loc[0,'age']='P12'
        with self.assertRaises(ValueError):audit_design(d)

    def test_missing_candidate_variant_is_not_stability(self):
        variants=['main','pseudo_0p1','pseudo_1','exclude_ambiguous_pair']+[f'omit_GSM{i}' for i in range(40)]
        d=pd.DataFrame([dict(variant=v,symbol=s,delta_P18=-1.,BY=.01,BH=.005) for v in variants for s in ['A','B']])
        d=d[~((d.symbol=='B')&(d.variant=='omit_GSM0'))]
        candidates=pd.DataFrame({'gene':['A','B','C'],'descriptive_direction':['down']*3})
        result=summarize_candidates(d,candidates,{'A','B','C'}).set_index('symbol')
        self.assertTrue(result.loc['A','all_variants_down_1p5fold'])
        self.assertFalse(result.loc['B','all_variants_down_1p5fold'])
        self.assertFalse(result.loc['C','all_variants_measured'])

    def test_conditional_fdr_over_full_gene_family(self):
        # Candidate-only BH would report .01; the complete four-gene family gives .04.
        from scipy.stats import false_discovery_control
        p=np.array([.01,.7,.8,.9])
        self.assertAlmostEqual(false_discovery_control(p,method='bh')[0],.04)
        self.assertGreater(false_discovery_control(p,method='by')[0],.04)

    def test_R_design_effect_sign_filter_and_nested_F(self):
        code=r'''
source('code/analyze_setx_interaction.R')
meta<-data.table(gsm=paste0('GSM',1:40),age=rep(rep(ages,each=2),each=4),group=rep(rep(c('Control','KO'),each=4),5))
x<-matrix(0,3,40,dimnames=list(c('constant_offset','late_change','late_only'),meta$gsm))
# Add equal developmental shifts to both groups; these must cancel in interactions.
for(i in 1:40){
 a<-match(meta$age[i],ages);ko<-meta$group[i]=='KO'
 x[1,i]<-2^(2+a+.4*ko)-.5
 x[2,i]<-2^(2+a+ifelse(a==5,-1*ko,0))-.5
 x[3,i]<-if(a==5)2 else 0
}
f<-fit_interaction(x,meta,.5,FALSE)$table
stopifnot(max(abs(as.numeric(f[gene=='constant_offset',paste0('delta_',ages[-1]),with=FALSE])))<1e-10)
stopifnot(abs(f[gene=='late_change']$delta_P18+1)<1e-10,'late_only' %in% f$gene)
stopifnot(identical(colnames(make_design(meta)),cells),qr(make_contrasts()$interaction)$rank==4)
# Enough non-degenerate synthetic genes for robust empirical Bayes and nested-F audit.
set.seed(17);z<-matrix(rexp(200*40),200,40,dimnames=list(paste0('g',1:200),meta$gsm))
z[,meta$age=='P18'&meta$group=='KO']<-z[,meta$age=='P18'&meta$group=='KO']*4
full<-fit_interaction(z,meta)
stopifnot(full$nested_F_max_relative_error<1e-8,full$residual_df==30)
stopifnot(max(abs(full$table$BY-p.adjust(full$table$P,'BY')))<1e-14)
keep<-!meta$gsm %in% c('GSM26','GSM30');small<-fit_interaction(z[,keep],meta[keep])
stopifnot(small$samples==38,small$rank==10,small$residual_df==28)
bad<-try(fit_interaction(z[,40:1],meta),silent=TRUE);stopifnot(inherits(bad,'try-error'))
'''
        result=subprocess.run(['bash','envs/rscript.sh','-e',code],cwd=ROOT,text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)


if __name__=='__main__':unittest.main()
