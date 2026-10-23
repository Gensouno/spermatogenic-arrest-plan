#!/usr/bin/env python3
"""Import reviewed expansion counts separately from the frozen eight-cohort baseline."""
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
RAW = ROOT / 'data/raw/expansion_assets'
OUT = ROOT / 'data/interim/expansion'
RESULTS = ROOT / 'results/data_expansion'


def validate_counts(x):
    values = x.to_numpy(dtype=float)
    if (x.index.has_duplicates or x.columns.has_duplicates or not np.isfinite(values).all()
            or np.any(values < 0) or np.any(np.abs(values - np.round(values)) > 1e-6)
            or np.any(values.sum(axis=0) <= 0)):
        raise ValueError('Invalid or non-integer counts; do not round unverified normalized abundance')


def parse_single_table(content):
    text = content.decode()
    sep = '\t' if '\t' in text.splitlines()[0] else ','
    table = pd.read_csv(io.StringIO(text), sep=sep, index_col=0)
    if table.shape[1] != 1:
        raise ValueError('Expected one expression column per sample')
    return table.iloc[:, 0]


def featurecounts(content):
    table = pd.read_csv(io.BytesIO(content), sep='\t', comment='#', index_col=0)
    if list(table.columns[:5]) != ['Chr', 'Start', 'End', 'Strand', 'Length'] or table.shape[1] != 6:
        raise ValueError('Unexpected featureCounts annotation/count columns')
    return table.iloc[:, 5]


def source(gse):
    return json.loads((RESULTS / f'{gse}_metadata.json').read_text())


def sample(gse, gsm, group, cohort, role, age, expected_genotype):
    record = source(gse)['samples'][gsm]
    characteristics = record['characteristics_ch1']
    genotype = next((v.split(': ', 1)[1] for v in characteristics if v.lower().startswith('genotype: ')), '')
    if genotype != expected_genotype:
        raise ValueError(f'Unexpected genotype for {gsm}: {genotype}')
    return {'cohort': cohort, 'gse': gse, 'sample_id': gsm, 'group': group,
            'genotype': genotype, 'age': age, 'role': role, 'official_title': ' | '.join(record['title']),
            'official_source': ' | '.join(record['source_name_ch1']),
            'official_characteristics': ' | '.join(characteristics),
            'metadata_sha256': source(gse)['sha256']}


def save(cohort, x, meta, id_type):
    validate_counts(x)
    if list(x.columns) != [r['sample_id'] for r in meta]:
        raise ValueError('Count column and GSM order differ')
    x.index.name = 'input_id'
    x.to_csv(OUT / f'{cohort}_raw.tsv.gz', sep='\t')
    pd.DataFrame(meta).assign(id_type=id_type).to_csv(OUT / f'{cohort}_samples.tsv', sep='\t', index=False)
    print(cohort, x.shape, meta[0]['role'], flush=True)


def audit_nkapl_archive():
    rows = []
    with tarfile.open(RAW / 'GSE232415/GSE232415_RAW.tar') as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            gsm = re.search(r'GSM\d+', member.name).group()
            expected = re.search(r'_(.+?)\.counts?\.csv', member.name).group(1)
            vector = parse_single_table(gzip.decompress(archive.extractfile(member).read()))
            # Ignore the column label when detecting duplicated sample vectors.
            payload = vector.rename('value').sort_index().to_csv().encode()
            rows.append({'gsm': gsm, 'member': member.name, 'expected_internal_label': expected,
                'actual_internal_label': vector.name, 'header_matches_filename': vector.name == expected,
                'profile_sha256': hashlib.sha256(payload).hexdigest(), 'genes': len(vector)})
    audit = pd.DataFrame(rows)
    audit['duplicated_profile'] = audit.profile_sha256.duplicated(keep=False)
    audit.to_csv(RESULTS / 'GSE232415_file_identity_audit.csv', index=False)
    return audit


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT / 'data/expansion_assets.json').read_text())
    for row in manifest['assets']:
        path = ROOT / row['path']
        if path.stat().st_size != row['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError(f'Expansion input changed: {row["path"]}')
        metadata = ROOT / f'data/raw/expansion_metadata/{row["gse"]}_family.soft.gz'
        if hashlib.sha256(metadata.read_bytes()).hexdigest() != row['metadata_sha256'] or source(row['gse'])['sha256'] != row['metadata_sha256']:
            raise ValueError(f'Expansion source metadata changed: {row["gse"]}')
    # Explicit source-table count columns: never accidentally import FPKM.
    table = pd.read_csv(RAW / 'GSE259242/GSE259242_raw_counts.txt.gz', sep='\t')
    aliases = [f'D2-9d-Control-{i}_count' for i in range(1, 4)] + [f'D2-9d-cko-{i}_count' for i in range(1, 4)]
    x = table.set_index('id')[aliases].copy()
    meta = [sample('GSE259242', f'GSM{8111945+i}', 'Control' if i < 3 else 'KO',
                   'GSE259242_P9', 'new_arrest_comparison', 'P9',
                   'control' if i < 3 else 'conditional Dcaf2 knockout') for i in range(6)]
    x.columns = [r['sample_id'] for r in meta]
    save('GSE259242_P9', x, meta, 'ensembl_gene')
    with tarfile.open(RAW / 'GSE117706/GSE117706_RAW.tar') as archive:
        members = {re.search(r'GSM\d+', m.name).group(): m for m in archive.getmembers() if m.isfile()}
        for age, first in [('P2', 3307579), ('P7', 3307585)]:
            cohort = f'GSE117706_{age}'
            meta, vectors = [], []
            for i in range(6):
                gsm = f'GSM{first+i}'
                vector = featurecounts(gzip.decompress(archive.extractfile(members[gsm]).read()))
                expected = ('control' if i < 3 else 'KO') + f'_{age}_{i%3+1}'
                if expected not in vector.name or expected not in members[gsm].name:
                    raise ValueError(f'featureCounts source sample differs: {gsm}')
                meta.append(sample('GSE117706', gsm, 'Control' if i < 3 else 'KO', cohort,
                                   'HET_control_sensitivity', age, 'Pramef12+/-' if i < 3 else 'Pramef12-/-'))
                vectors.append(vector.rename(gsm))
            save(cohort, pd.concat(vectors, axis=1), meta, 'ensembl_gene')
    with tarfile.open(RAW / 'GSE235777/GSE235777_RAW.tar') as archive:
        vectors, meta, summaries = [], [], []
        for member in sorted(archive.getmembers(), key=lambda m: m.name):
            if not member.isfile():
                continue
            gsm = re.search(r'GSM\d+', member.name).group()
            index = int(gsm[3:]) - 7509320
            if index not in range(6):
                raise ValueError('Unexpected Msl3 sample')
            table = pd.read_csv(io.BytesIO(gzip.decompress(archive.extractfile(member).read())), sep='\t', header=None, index_col=0)
            if table.shape[1] != 3:
                raise ValueError('STAR gene count table must contain three strand modes')
            gene = table.index.str.startswith('ENSMUSG')
            summaries.append({'gsm': gsm, 'unstranded_gene_total': float(table.loc[gene].iloc[:, 0].sum()),
                'strand1_gene_total': float(table.loc[gene].iloc[:, 1].sum()),
                'strand2_gene_total': float(table.loc[gene].iloc[:, 2].sum())})
            vectors.append(table.loc[gene].iloc[:, 0].rename(gsm))
            meta.append(sample('GSE235777', gsm, 'Control' if index < 3 else 'KO', 'GSE235777_P45',
                'phenotype_negative_comparator', 'P45', 'MSL3 WT' if index < 3 else 'MSL3 cKO'))
        pd.DataFrame(summaries).to_csv(RESULTS / 'Msl3_STAR_strand_totals.csv', index=False)
        save('GSE235777_P45', pd.concat(vectors, axis=1), meta, 'ensembl_gene')
    audit = audit_nkapl_archive()
    print('NKAPL/SOX30 whole-testis files BLOCKED: internal-label mismatch and duplicated profiles;',
          int((~audit.header_matches_filename).sum()), 'mismatches,', int(audit.duplicated_profile.sum()), 'profiles in duplicate groups', flush=True)
    unit_checks = []
    for gse, reader in [('GSE278121', lambda p: pd.read_csv(p, sep='\t', index_col=0).dropna(axis=1, how='all')),
                        ('GSE290783', lambda p: pd.read_csv(p, sep=r'\s+', index_col=0))]:
        x = reader(RAW / f'{gse}/{gse}_raw_counts.txt.gz')
        a = x.to_numpy(float)
        unit_checks.append({'gse': gse, 'genes': len(x), 'samples': x.shape[1],
            'fraction_noninteger_values': float((np.abs(a - np.round(a)) > 1e-6).mean()),
            'decision': 'unit_pending; not rounded or analyzed as raw counts'})
    pd.DataFrame(unit_checks).to_csv(RESULTS / 'unresolved_count_units.csv', index=False)


if __name__ == '__main__':
    main()
