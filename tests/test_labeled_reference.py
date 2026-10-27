import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import io
import contextlib
import json
import tempfile
import hashlib

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'));sys.path.insert(0,str(ROOT/'scripts'))
from prepare_labeled_reference import broad_label,parse_triplets
import fetch_labeled_ranges as fetch_ranges
import prepare_labeled_reference as prepare
from check_author_zygotene import join_author_labels


class LabeledReferenceTests(unittest.TestCase):
    def test_author_crosswalk_is_explicit(self):
        for value,expected in [('Spermatogonia','Spermatogonia'),('eP1','Meiotic'),('D','Meiotic'),
            ('MI','Meiotic'),('MII','Meiotic'),('S1','Postmeiotic'),('S11','Postmeiotic'),('Sertoli','Somatic'),('Outliers',None)]:
            self.assertEqual(broad_label(value),expected)
        for value in ['S0','S12','NewCell','invented']:
            with self.assertRaises(ValueError):broad_label(value)

    def test_triplet_block_preserves_indices(self):
        np.testing.assert_array_equal(parse_triplets(b'1 1 3\n4 1 7\n2 2 5\n',4,2),[[1,1,3],[4,1,7],[2,2,5]])

    def test_triplets_reject_fractional_negative_truncated_or_zero_values(self):
        for b in [b'1 1 1.5\n',b'1 1 -2\n',b'1 1\n',b'1 1 0\n',b'1 1 2x\n',b'1 1 2\n2 2\n',b'0 1 2\n',b'5 1 1\n',b'1 1 1 2\n2 3\n']:
            with self.assertRaises(ValueError):parse_triplets(b,4,2)

    def test_triplets_reject_duplicate_and_unsorted_entries(self):
        for b in [b'1 1 2\n1 1 3\n',b'2 1 2\n1 1 3\n',b'1 2 2\n1 1 3\n']:
            with self.assertRaises(ValueError):parse_triplets(b,4,2)

    def test_http_range_must_be_honored_and_etag_stable(self):
        class Response(io.BytesIO):
            status=206
            headers={'Content-Range':f'bytes 0-9/{fetch_ranges.SIZE}','ETag':'original'}
        remote=fetch_ranges.RemoteMatrix()
        with patch.object(fetch_ranges.urllib.request,'urlopen',return_value=Response(b'1234567890')):
            with remote.request(0,9):pass
        self.assertEqual(remote.etag,'original')
        for headers in [{'Content-Range':'wrong','ETag':'original'},
            {'Content-Range':f'bytes 0-9/{fetch_ranges.SIZE}','ETag':'changed'}]:
            response=Response(b'1234567890');response.headers=headers
            with patch.object(fetch_ranges.urllib.request,'urlopen',return_value=response):
                with self.assertRaises(ValueError):remote.request(0,9)

    def test_byte_boundary_is_exact_not_proportional_guess(self):
        data=b'header\n'+b''.join(f'{g} {c} {g+c}\n'.encode() for c in range(1,101) for g in range(1,51))
        class LocalMatrix(fetch_ranges.RemoteMatrix):
            def request(self,start,end):return io.BytesIO(data[start:end+1])
        with patch.object(fetch_ranges,'SIZE',len(data)):
            r=LocalMatrix();cumulative=np.arange(101)*50
            for c in [1,2,10,55,99]:
                result=r.boundary(c,cumulative,7)
                self.assertTrue(data[result:].startswith(f'1 {c} '.encode()))

    def test_small_end_to_end_aggregation_preserves_units_and_excludes_ambiguous_species(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);base=root/'data/raw/labeled_reference';base.mkdir(parents=True)
            genes=pd.DataFrame({'ID':['ENSMUSG00001','ENSMUSG00002','ENSMUSG00003','ENSMUSG00004','ENSG00001'],
                'Symbol':['GeneA','GeneB','Ambiguous','Ambiguous','HumanGene']})
            genes.to_csv(base/'genes.tsv',sep='\t',index=False)
            rows=[];sources=[];assets=[]
            for library,sample,label,age in [('do17815','B6','S9',67),('do26386','P5','Spermatogonia',5),('do18195','P15','eP2',15)]:
                start=len(rows)+1
                for cell in range(21):rows.append({'Sample':sample,'Library':library,'Barcode':str(cell),
                    'AnnotatedClusters':label,'total_counts':15,'total_features_by_counts':5})
                end=len(rows);name=library+'.triplets'
                (base/name).write_text(''.join(f'{g} {c} {g}\n' for c in range(start,end+1) for g in range(1,6)))
                assets.append({'path':str((base/name).relative_to(root)),
                    'range':{'library':library,'first_column':start,'last_column':end,'nonzeros':105}})
                sources.append({'Source Name':library,'Characteristics[age]':age,'Characteristics[strain]':'C57BL6J',
                    'Characteristics[genotype]':'wild type genotype','Comment[cell count]':21})
            cells=pd.DataFrame(rows,index=range(1,64));cells.to_csv(base/'cell_metadata.txt',sep=' ')
            pd.DataFrame(sources).to_csv(base/'E-MTAB-6946.sdrf.txt',sep='\t',index=False)
            (root/'data/labeled_reference_assets.json').write_text(json.dumps({'assets':assets,
                'source_matrix':{'header':['%%MatrixMarket matrix coordinate integer general','5 63 315']}}))
            with patch.object(prepare,'ROOT',root),patch.object(prepare,'WT_SAMPLES',{'B6','P5','P15'}),contextlib.redirect_stdout(io.StringIO()):
                prepare.main()
            result=pd.read_csv(root/'results/labeled_reference/type_profiles.csv.gz',index_col=0)
            self.assertEqual(result.index.tolist(),['GeneA','GeneB'])
            self.assertEqual(result.shape,(2,6))
            self.assertEqual(result.loc['GeneA','do17815:S9:RNA_sum'],21)
            self.assertAlmostEqual(result.loc['GeneA','do17815:S9:equal_cell_mean'],1e6/15)
            cells.loc[1,'total_counts']=16;cells.to_csv(base/'cell_metadata.txt',sep=' ')
            with patch.object(prepare,'ROOT',root),patch.object(prepare,'WT_SAMPLES',{'B6','P5','P15'}),contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'author per-cell totals'):prepare.main()

    def test_fine_annotation_join_preserves_original_matrix_columns_and_library(self):
        original=pd.DataFrame({'Library':['A','B','A'],'Barcode':['same','same','other']},index=[10,20,30])
        fine=pd.DataFrame({'Library':['B','A'],'Barcode':['same','same'],'P15Clusters':['Zygotene','Leptotene']})
        joined=join_author_labels(original,fine)
        self.assertEqual(joined.original_matrix_column.tolist(),[10,20,30])
        self.assertEqual(joined.P15Clusters.iloc[:2].tolist(),['Leptotene','Zygotene'])
        self.assertTrue(pd.isna(joined.P15Clusters.iloc[2]))
        with self.assertRaises(ValueError):join_author_labels(original,pd.concat([fine,fine.iloc[:1]]))

    def test_missing_frozen_window_restores_only_matching_content(self):
        data=b'1 1 2\n'
        class Remote:
            def request(self,start,end):return io.BytesIO(data)
        with tempfile.TemporaryDirectory() as directory:
            cache=Path(directory)
            row={'url':fetch_ranges.URL,'etag':'fixed','bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),
                'range':{'library':'do17815','start_byte':71,'end_byte_exclusive':71+len(data)}}
            with patch.object(fetch_ranges,'CACHE',cache),patch.dict(fetch_ranges.os.environ,{'HTTPS_PROXY':'http://proxy.invalid','https_proxy':'http://proxy.invalid'}),contextlib.redirect_stdout(io.StringIO()):
                fetch_ranges.restore_window(row,Remote())
                self.assertEqual((cache/'do17815.triplets').read_bytes(),data)
                row['range']['library']='do26386';row['sha256']='0'*64
                with self.assertRaisesRegex(ValueError,'content differs'):fetch_ranges.restore_window(row,Remote())
                self.assertFalse((cache/'do26386.triplets').exists())


if __name__=='__main__':unittest.main()
