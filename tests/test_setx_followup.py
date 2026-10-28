"""Identity ambiguity and length-aware count-model scientific safeguards."""
from pathlib import Path
import subprocess
import sys
import unittest
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from audit_setx_sources import parse_soft,exact_relation,reconcile


class SourceTests(unittest.TestCase):
    def fixture(self):
        text='''^SERIES = GSE1
!Series_title = Test
^SAMPLE = GSM1
!Sample_title = wt testis, 16dpp, rep2
!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN1
!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX1
'''
        _,samples=parse_soft(text)
        ena=pd.DataFrame([dict(run_accession='SRR1',experiment_accession='SRX1',sample_accession='SAMN1',
            sample_alias='GSM1',secondary_sample_accession='SRS1',sample_title='wt testis, 16dpp, rep2',
            library_name='GSM1',library_layout='PAIRED',library_strategy='RNA-Seq',instrument_model='HiSeq',
            tax_id='10090',read_count='10',base_count='1000')])
        design=pd.DataFrame([dict(gsm='GSM1',age='P16',group='Control',source_code='7844',barcode_token='AGTTCC')])
        return samples,ena,design

    def test_accessions_do_not_establish_independent_animals(self):
        samples,ena,design=self.fixture();d=reconcile(samples,samples,ena,1,design)
        self.assertTrue(d.source_identity_reconciled.all())
        self.assertFalse(d.individual_animal_identity_confirmed.any())
        self.assertEqual(d.filename_source_code.iloc[0],'7844')

    def test_duplicate_run_and_count_mismatch_fail(self):
        samples,ena,design=self.fixture()
        with self.assertRaises(ValueError):reconcile(samples,samples,pd.concat([ena,ena]),2,design)
        with self.assertRaises(ValueError):reconcile(samples,samples,ena,2,design)

    def test_wrong_label_or_biosample_fails(self):
        samples,ena,design=self.fixture();ena.loc[0,'sample_title']='SetxKO testis, 16dpp, rep2'
        with self.assertRaises(ValueError):reconcile(samples,samples,ena,1,design)
        samples,ena,design=self.fixture();ena.loc[0,'sample_accession']='SAMN2'
        with self.assertRaises(ValueError):reconcile(samples,samples,ena,1,design)

    def test_duplicate_sample_and_ambiguous_relation_fail(self):
        with self.assertRaises(ValueError):parse_soft('^SAMPLE = GSM1\n^SAMPLE = GSM1\n')
        with self.assertRaises(ValueError):exact_relation({'!Sample_relation':['SRA: SRX1','SRA: SRX2']},'SRA',r'SRX\d+')

    def test_estimated_counts_offsets_and_known_interaction(self):
        code=r'''
source('code/analyze_setx_length.R')
set.seed(72)
meta<-data.table(gsm=paste0('GSM',1:40),age=rep(rep(ages,each=2),each=4),group=rep(rep(c('Control','KO'),each=4),5))
x<-matrix(rnbinom(400*40,mu=500,size=30)+.25,400,40,dimnames=list(paste0('g',1:400),meta$gsm))
l<-matrix(200,400,40,dimnames=dimnames(x))
late<-meta$age=='P18'&meta$group=='KO';x[1,late]<-x[1,late]/16
# Deliberately invalid expressed gene length must be excluded, never imputed.
l[2,1]<-0
f<-length_fit(x,l,meta)
stopifnot(!('g2' %in% f$table$gene),f$summary$excluded_nonpositive_length==1)
stopifnot(f$summary$noninteger_source_values>0,!f$summary$counts_rounded,f$summary$length_offset_max_error<1e-12)
stopifnot(f$summary$interaction_df==4,f$table[gene=='g1']$delta_P18 < -2)
stopifnot(max(abs(f$table$BY-p.adjust(f$table$P,'BY')))<1e-14)
bad<-try(length_fit(x[,40:1],l[,40:1],meta),silent=TRUE);stopifnot(inherits(bad,'try-error'))
'''
        p=subprocess.run(['bash','envs/rscript.sh','-e',code],cwd=ROOT,text=True,capture_output=True)
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)


if __name__=='__main__':unittest.main()
