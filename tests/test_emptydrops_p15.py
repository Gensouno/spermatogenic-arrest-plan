"""Failure-mode checks for source boundaries, count reconciliation and strata."""
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'code')]
from emptydrops_range_client import MatrixRanges
from fetch_emptydrops_p15 import verify_frozen
from prepare_emptydrops_p15 import aggregate,broad_label,mark_overlap,unique_mouse,verify_cells


class Response(io.BytesIO):
    def __init__(self,body,start,total,etag='"fixed"',status=206):
        super().__init__(body);self.status=status
        self.headers={'ETag':etag,'Content-Range':f'bytes {start}-{start+len(body)-1}/{total}'}


class EmptyDropsTests(unittest.TestCase):
    def test_fixed_labels(self):
        self.assertEqual(broad_label('Leptotene'),'Meiotic')
        self.assertEqual(broad_label('Sertoli_2'),'Somatic')
        self.assertIsNone(broad_label('Outliers'))
        with self.assertRaises(ValueError):broad_label('NewCell')

    def test_library_barcode_pair(self):
        old=pd.DataFrame({'Library':['a'],'Barcode':['same']})
        new=pd.DataFrame({'Library':['a','b'],'Barcode':['same','same']})
        np.testing.assert_array_equal(mark_overlap(new,old),[True,False])
        with self.assertRaises(ValueError):mark_overlap(new,pd.concat([old,old]))

    def test_ambiguous_species_filter(self):
        genes=pd.DataFrame({'ID':['ENSMUSG1','ENSMUSG2','ENSMUSG3','ENSG4'],'Symbol':['dup','dup','ok','human']})
        np.testing.assert_array_equal(unique_mouse(genes),[False,False,True,False])
        genes.loc[0,'ID']='bad'
        with self.assertRaises(ValueError):unique_mouse(genes)

    def test_count_reconciliation(self):
        cells=pd.DataFrame({'total_counts':[10,20],'total_features_by_counts':[1,2]})
        verify_cells(np.array([10,20]),np.array([1,2]),cells)
        with self.assertRaises(ValueError):verify_cells(np.array([10,21]),np.array([1,2]),cells)
        with self.assertRaises(ValueError):verify_cells(np.array([10,20]),np.array([2,2]),cells)

    def test_streaming_strata_and_means(self):
        cells=pd.DataFrame({'total_counts':[2,8],'total_features_by_counts':[1,2]})
        groups=pd.DataFrame({'cells':[2,1,1]})
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'counts';path.write_bytes(b'1 1 2\n1 2 2\n2 2 6\n')
            sums,means,_,_=aggregate(path,pd.DataFrame(index=range(2)),cells,1,2,3,groups,
                [np.array([0,0]),np.array([1,-1]),np.array([-1,2])])
            np.testing.assert_array_equal(sums[:,0],sums[:,1]+sums[:,2])
            np.testing.assert_allclose(means[:,0],[625000,375000])
            np.testing.assert_allclose(means[:,0],(means[:,1]+means[:,2])/2)
            with self.assertRaises(ValueError):aggregate(path,pd.DataFrame(index=range(2)),cells,1,2,4,groups,[np.array([0,0])])

    def test_proxy_required(self):
        remote=MatrixRanges('https://ftp.ebi.ac.uk/biostudies/a',100)
        with patch.dict(os.environ,{},clear=True),self.assertRaises(ValueError):remote.request(0,9)

    def test_range_etag_and_status(self):
        remote=MatrixRanges('https://ftp.ebi.ac.uk/biostudies/a',100,etag='"fixed"')
        env={k:'http://proxy' for k in ['HTTP_PROXY','http_proxy','HTTPS_PROXY','https_proxy']}
        for tag,status in [('"changed"',206),('W/"fixed"',206),('"fixed"',200)]:
            with patch.dict(os.environ,env),patch('urllib.request.urlopen',return_value=Response(b'123',0,100,tag,status)),self.assertRaises(ValueError):remote.request(0,2)
        with patch.dict(os.environ,env),patch('urllib.request.urlopen',return_value=Response(b'123',0,100)) as mock:
            remote.request(0,2).close()
            self.assertEqual(mock.call_args.args[0].get_header('If-match'),'"fixed"')

    def test_boundary_exact_column_transition(self):
        body=b'%%MatrixMarket matrix coordinate integer general\n2 5 10\n'+b''.join(f'{g} {c} {g+c}\n'.encode() for c in range(1,6) for g in [1,2])
        header=body.index(b'1 1 2')
        remote=MatrixRanges('https://ftp.ebi.ac.uk/biostudies/a',len(body))
        def request(start,end):return io.BytesIO(body[start:end+1])
        with patch.object(remote,'request',side_effect=request):
            self.assertEqual(remote.boundary(1,np.arange(0,11,2),header),header)
            self.assertEqual(remote.boundary(6,np.arange(0,11,2),header),len(body))
            self.assertEqual(remote.boundary(3,np.arange(0,11,2),header),body.index(b'1 3 4'))

    def test_verified_chunk_resume_and_corruption(self):
        body=b'1 1 5\n';remote=MatrixRanges('https://ftp.ebi.ac.uk/biostudies/a',len(body),etag='"fixed"')
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'window'
            with patch.object(remote,'request',return_value=io.BytesIO(body)) as mock:
                result=remote.download(path,0,len(body)-1);self.assertEqual(mock.call_count,1)
            with patch.object(remote,'request',side_effect=AssertionError('Network forbidden')):remote.download(path,0,len(body)-1)
            self.assertEqual(result['sha256'],hashlib.sha256(body).hexdigest())
            (path.with_name(path.name+'.parts')/'0000.bin').write_bytes(b'1 1 6\n')
            with self.assertRaises(ValueError):remote.download(path,0,len(body)-1)

    def test_frozen_cache_never_silently_replaced(self):
        import fetch_emptydrops_p15 as module
        with tempfile.TemporaryDirectory() as temp,patch.object(module,'ROOT',Path(temp)):
            row={'path':'input','bytes':3,'sha256':hashlib.sha256(b'yes').hexdigest()}
            with self.assertRaises(FileNotFoundError):verify_frozen(row)
            (Path(temp)/'input').write_bytes(b'bad')
            with self.assertRaises(ValueError):verify_frozen(row,True)

    def test_missing_frozen_cache_restores_and_verifies(self):
        import fetch_emptydrops_p15 as module
        env={k:'http://proxy' for k in ['HTTP_PROXY','http_proxy','HTTPS_PROXY','https_proxy']}
        with tempfile.TemporaryDirectory() as temp,patch.object(module,'ROOT',Path(temp)),patch.dict(os.environ,env):
            row={'path':'input','url':'https://ftp.ebi.ac.uk/biostudies/a','bytes':3,'sha256':hashlib.sha256(b'yes').hexdigest()}
            with patch('urllib.request.urlopen',return_value=io.BytesIO(b'yes')):verify_frozen(row,True)
            self.assertEqual((Path(temp)/'input').read_bytes(),b'yes')
            (Path(temp)/'input').unlink()
            with patch('urllib.request.urlopen',return_value=io.BytesIO(b'bad')),self.assertRaises(ValueError):verify_frozen(row,True)
            self.assertFalse((Path(temp)/'input').exists())

    def test_chunk_directory_symlink_refused(self):
        remote=MatrixRanges('https://ftp.ebi.ac.uk/biostudies/a',10,etag='"fixed"')
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'window';(Path(temp)/'target').mkdir()
            path.with_name('window.parts').symlink_to(Path(temp)/'target',target_is_directory=True)
            with self.assertRaises(ValueError):remote.download(path,0,9)


if __name__=='__main__':unittest.main()
