#!/usr/bin/env python3
"""Source-checked sample aliases and audit-only FPKM import."""
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


def lrh1_samples(source):
    if len(source) != 16:
        raise ValueError('Expected 16 official Lrh1 samples')
    rows = []
    for i in range(16):
        gsm = f'GSM{5455535+i}'
        record = source[gsm]
        p1 = i < 8
        age = 'P1' if p1 else 'P3'
        group = 'Control' if i % 8 < 4 else 'KO'
        driver = 'NanosCre' if p1 else 'DhhCre'
        alias = f'{age.lower()}_testis_{"noCre" if group == "Control" else driver}_rep{i%4+1}'
        genotype = ('Lrh1 flox/-' if p1 else 'Lrh1 flox/flox')
        if group == 'KO':
            genotype += '; ' + ('Nanos3-Cre' if p1 else 'Dhh-Cre')
        chars = record['Sample_characteristics_ch1']
        title = record['Sample_title'][0]
        match = re.search(r'\[([^\]]+)\]$', title)
        if (not match or match[1].replace('Nanos3Cre', 'NanosCre') != alias
                or 'genotype: ' + genotype not in chars
                or 'age: Postnatal Day ' + age[1:] not in chars
                or 'tissue: testis' not in chars
                or record.get('Sample_organism_ch1') != ['Mus musculus']):
            raise ValueError('Official sample evidence does not match prespecified alias: ' + gsm)
        rows.append({'gsm': gsm, 'cohort': f'Lrh1_{age}_{driver}', 'age': age, 'group': group,
                     'author_sample_name': alias, 'genotype': genotype,
                     'official_title': title, 'official_characteristics': ' | '.join(chars),
                     'alias_note': 'Only Nanos3Cre -> NanosCre spelling alias; no expression-based relabeling.'})
    return pd.DataFrame(rows)


def parse_fpkm(content):
    if not content.startswith(b'\xff\xfe'):
        raise ValueError('Expected source UTF-16 little-endian BOM')
    frame = pd.read_csv(io.StringIO(content.decode('utf-16')), sep='\t', index_col=0)
    if frame.shape[1] != 1 or not frame.columns[0].endswith('.FPKM'):
        raise ValueError('Expected a single explicitly labeled FPKM column')
    a = frame.iloc[:, 0]
    if (a.index.has_duplicates or not a.index.str.startswith('ENSMUST').all()
            or not np.isfinite(a).all() or (a < 0).any()):
        raise ValueError('Invalid transcript FPKM table')
    return a


def main():
    out = ROOT / 'results/followup'
    metadata = json.loads((out / 'metadata.json').read_text())
    lrh1_samples(metadata['GSE180204']['samples']).to_csv(out / 'Lrh1_source_samples.csv', index=False)
    rows = []
    with tarfile.open(ROOT / 'data/raw/followup/GSE181132/GSE181132_RAW.tar') as archive:
        members = [m for m in archive.getmembers() if m.isfile()]
        if len(members) != 4:
            raise ValueError('Expected four Ppp2ca processed profiles')
        for member in members:
            gsm = re.search(r'GSM\d+', member.name)[0]
            source = metadata['GSE181132']['samples'][gsm]
            listed = [v for k, values in source.items() if k.startswith('Sample_supplementary_file') for v in values]
            if not any(u.rsplit('/', 1)[-1] == member.name for u in listed):
                raise ValueError('Archive member differs from official sample attachment')
            vector = parse_fpkm(gzip.decompress(archive.extractfile(member).read()))
            expected = re.search(r'_([ZM]-\d+)\.txt', member.name)[1] + '.FPKM'
            if vector.name != expected:
                raise ValueError('Internal sample label differs from filename')
            a = vector.to_numpy(float)
            rows.append({'gsm': gsm, 'member': member.name, 'internal_column': vector.name,
                'features': len(vector), 'id_type': 'Ensembl transcript', 'unit': 'FPKM',
                'fraction_noninteger': float((np.abs(a - np.round(a)) > 1e-6).mean()),
                'official_characteristics': ' | '.join(source['Sample_characteristics_ch1']),
                'profile_sha256': hashlib.sha256(vector.rename('value').sort_index().to_csv().encode()).hexdigest(),
                'decision': 'Audit only; not admitted as raw gene counts'})
    frame = pd.DataFrame(rows)
    if frame.gsm.nunique() != 4 or frame.profile_sha256.nunique() != 4:
        raise ValueError('Missing/duplicated Ppp2ca profiles')
    frame.to_csv(out / 'Ppp2ca_unit_audit.csv', index=False)
    print(frame[['gsm', 'features', 'unit', 'fraction_noninteger']].to_string(index=False))


if __name__ == '__main__':
    main()
