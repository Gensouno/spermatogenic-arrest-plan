"""Missing-aware omission retention and controlled length-offset diagnostics."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from run_setx_stability import candidate_stability,check_inputs


class StabilityTests(unittest.TestCase):
    def example(self):
        required=['full','exclude_pair','omit_GSM1','omit_GSM2']
        d=pd.DataFrame([dict(variant=v,symbol=s,delta_P18=-1.,BY=.01) for v in required for s in ['A','B']])
        return required,d

    def test_missing_variant_or_symbol_is_not_retention(self):
        required,d=self.example();d=d[~((d.variant=='omit_GSM1')&(d.symbol=='B'))]
        out=candidate_stability(d,['A','B','C'],{'A','B','C'},required).set_index('symbol')
        self.assertTrue(out.loc['A','all_negative_1p5fold'])
        self.assertFalse(out.loc['B','all_variants_measured'])
        self.assertFalse(out.loc['B','conditional_BY05_all_variants_descriptive'])
        self.assertFalse(out.loc['C','all_negative'])
        with self.assertRaises(ValueError):candidate_stability(d[d.variant!='omit_GSM2'],['A','B'],{'A'},required)

    def test_effect_retention_does_not_imply_adjusted_threshold_retention(self):
        required,d=self.example();d.loc[(d.symbol=='A')&(d.variant=='omit_GSM1'),'BY']=.2
        out=candidate_stability(d,['A'],{'A'},required).iloc[0]
        self.assertTrue(out.all_negative_1p5fold)
        self.assertFalse(out.conditional_BY05_all_variants_descriptive)

    def test_zero_sign_or_insufficient_magnitude_fails(self):
        required,d=self.example();d.loc[(d.symbol=='A')&(d.variant=='omit_GSM1'),'delta_P18']=0
        d.loc[(d.symbol=='B')&(d.variant=='omit_GSM1'),'delta_P18']=-.1
        out=candidate_stability(d,['A','B'],{'A','B'},required).set_index('symbol')
        self.assertFalse(out.loc['A','all_negative']);self.assertTrue(out.loc['B','all_negative'])
        self.assertFalse(out.loc['B','all_negative_1p5fold'])

    def test_duplicate_candidate_fit_is_rejected(self):
        required,d=self.example()
        with self.assertRaises(ValueError):candidate_stability(pd.concat([d,d.iloc[:1]]),['A','B'],{'A'},required)

    def test_input_manifest_rejects_path_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'data').mkdir()
            (root/'data/setx_stability_inputs.json').write_text(json.dumps({'assets':[{'path':'../escape','bytes':1,'sha256':'0'*64}]}))
            with self.assertRaises(ValueError):check_inputs(root)

    def test_matched_count_models_hold_normalization_fixed(self):
        code=r'''
source('code/analyze_setx_stability.R')
set.seed(28)
meta<-data.table(gsm=paste0('GSM',1:40),age=rep(rep(ages,each=2),each=4),group=rep(rep(c('Control','KO'),each=4),5))
x<-matrix(rnbinom(250*40,mu=1000,size=40)+.25,250,40,dimnames=list(paste0('g',1:250),meta$gsm))
l<-matrix(200,250,40,dimnames=dimnames(x))
constant<-matched_count_pair(x,l,meta)
a<-constant$tables[['matched_count_length_offset']];b<-constant$tables[['matched_count_zero_offset_counterfactual']]
# Centering constant log-lengths leaves machine-roundoff offsets (~1e-15),
# propagated slightly by the likelihood solver; use a numerical, not exact, check.
stopifnot(max(abs(a$delta_P18-b$delta_P18))<1e-9,max(abs(a$P-b$P))<1e-9,constant$normalization_error<1e-12)
# Exactly eightfold count/length change in one age/genotype cell.
late<-meta$age=='P18'&meta$group=='KO';x[1,]<-1000.25;x[1,late]<-8*x[1,late];l[1,late]<-8*l[1,late]
changed<-matched_count_pair(x,l,meta)
a<-changed$tables[['matched_count_length_offset']];b<-changed$tables[['matched_count_zero_offset_counterfactual']]
stopifnot(abs((b[gene=='g1']$delta_P18-a[gene=='g1']$delta_P18)-3)<.1)
stopifnot(changed$normalization_error<1e-12)
stopifnot(max(abs(a$BY-p.adjust(a$P,'BY')))<1e-14)
bad<-try(matched_count_pair(x[,40:1],l[,40:1],meta),silent=TRUE);stopifnot(inherits(bad,'try-error'))
'''
        p=subprocess.run(['bash','envs/rscript.sh','-e',code],cwd=ROOT,text=True,capture_output=True)
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)


if __name__=='__main__':unittest.main()
