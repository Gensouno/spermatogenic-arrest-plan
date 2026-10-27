"""New data scale, genome/piRNA identity and source-symbol safeguards."""
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts')]
from prepare_new_models import read_eef,select_genome_members,EEF_COLUMNS
from summarize_new_models import mapped_effects


class NewModelTests(unittest.TestCase):
    def test_annotation_columns_not_read_as_counts(self):
        d=pd.DataFrame({'gene_id':['E1'],'gene_name':['A'],**{c:[i] for i,c in enumerate(EEF_COLUMNS)},'gene_length':[17.3]})
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'source.tsv';d.to_csv(p,sep='\t',index=False)
            self.assertEqual(list(read_eef(p)[EEF_COLUMNS].iloc[0]),list(range(6)))
            d[EEF_COLUMNS[0]]=.2;d.to_csv(p,sep='\t',index=False)
            with self.assertRaises(ValueError):read_eef(p)

    def test_extra_source_library_cannot_be_silently_dropped(self):
        d=pd.DataFrame({'gene_id':['E1'],'gene_name':['A'],**{c:[1] for c in EEF_COLUMNS},'Ee_con4':[1]})
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'source.tsv';d.to_csv(p,sep='\t',index=False)
            with self.assertRaises(ValueError):read_eef(p)

    def test_piRNA_excluded_but_inventory_must_reconcile(self):
        genome=tarfile.TarInfo('GSM1_x_GRCm38.trimmed.genes.results.gz')
        pirna=tarfile.TarInfo('GSM1_x_mm10_piRNA.trimmed.genes.results.gz')
        self.assertEqual(select_genome_members([genome,pirna],{'GSM1'}),{'GSM1':genome})
        for members in [[genome],[genome,pirna,pirna],[genome,pirna,tarfile.TarInfo('../escape')]]:
            with self.assertRaises(ValueError):select_genome_members(members,{'GSM1'})

    def test_symbol_ambiguity_excluded_without_summing_effects(self):
        d=pd.DataFrame({'gene':['E1','E2','E3','E4'],'logFC':[1,-1,2,3]})
        mapping=pd.DataFrame({'gene':['E1','E2','E3','E4'],'symbol':['A','A','B','-']})
        result=mapped_effects(d,mapping)
        self.assertEqual(result.to_dict(),{'B':2})
        with self.assertRaises(ValueError):mapped_effects(pd.concat([d,d.iloc[:1]]),mapping)


if __name__=='__main__':unittest.main()
