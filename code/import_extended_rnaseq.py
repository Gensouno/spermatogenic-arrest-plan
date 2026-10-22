#!/usr/bin/env python3
"""Read selected tar members in memory; keep count and TPM units distinct."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    sheet = pd.read_csv(ROOT / 'data/samples.tsv',sep='\t')
    dest = ROOT / 'data/interim/extended'
    dest.mkdir(parents=True,exist_ok=True)
    manifest = []
    for gse, suffix, unit in [('GSE101511','.count.txt.gz','transcript_feature_counts'),
                               ('GSE294209','_gene_tpm_estimates.txt.gz','gene_TPM')]:
        archive = ROOT / f'data/raw/{gse}/{gse}_RAW.tar'
        selected = sheet[(sheet.gse == gse) & (sheet.tissue_class == 'whole_testis')].copy()
        matrices, names = [], []
        with tarfile.open(archive) as tar:
            for member in tar.getmembers():
                if not member.isfile() or not member.name.endswith(suffix):
                    continue
                gsm = re.search(r'GSM\d+',member.name).group()
                if gsm not in selected.sample_id.tolist():
                    raise ValueError(f'Unexpected selected sample {gsm}')
                tab = pd.read_csv(io.BytesIO(gzip.decompress(tar.extractfile(member).read())),
                                  sep='\t',header=None,index_col=0)
                if tab.shape[1] != 1 or tab.index.has_duplicates:
                    raise ValueError('Ambiguous input rows')
                matrices.append(tab.iloc[:,0].rename(gsm))
                names.append(member.name)
        matrix = pd.concat(matrices,axis=1)
        if matrix.isna().any().any() or set(matrix.columns) != set(selected.sample_id):
            raise ValueError('Missing values or mismatched samples')
        matrix = matrix[selected.sample_id.tolist()]
        if not np.isfinite(matrix.to_numpy()).all() or matrix.to_numpy().min() < 0:
            raise ValueError('Invalid abundance')
        if unit.endswith('counts') and not np.allclose(matrix,np.round(matrix)):
            raise ValueError('Noninteger feature counts')
        selected['matrix_sample'] = selected.sample_id
        matrix.index.name = 'input_id'
        matrix.to_csv(dest / f'{gse}_raw.tsv.gz',sep='\t')
        selected.to_csv(ROOT / f'data/bulk/{gse}_samples.tsv',sep='\t',index=False)
        manifest.append({'gse':gse,'unit':unit,'rows':len(matrix),'samples':len(selected),
                         'selected_members':names,'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest()})
        print(gse, matrix.shape, unit, flush=True)
    (dest / 'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__ == '__main__':
    main()
