#!/usr/bin/env python3
"""New, untouched sorted WT profiles; fixed reference compatibility diagnostics."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tarfile

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from prepare_expansion import parse_single_table, source

ROOT = Path(__file__).resolve().parents[1]


def tar_profiles(gse, selected, unit):
    vectors, rows = [], []
    with tarfile.open(ROOT / f'data/raw/expansion_assets/{gse}/{gse}_RAW.tar') as archive:
        members = {re.search(r'GSM\d+', m.name).group(): m for m in archive.getmembers() if m.isfile()}
        for gsm, stage in selected.items():
            m = members[gsm]
            vector = parse_single_table(gzip.decompress(archive.extractfile(m).read()))
            record = source(gse)['samples'][gsm]
            if gse == 'GSE232415':
                expected = re.search(r'_(.+?)\.counts\.csv', m.name).group(1)
                if expected != vector.name or 'genotype: Wild-type' not in record['characteristics_ch1']:
                    raise ValueError('Unresolved sorted NKAPL sample identity')
                if not any(f'({stage})' in value for value in record['characteristics_ch1'] if value.startswith('cell type:')):
                    raise ValueError('Sorted NKAPL stage differs from source evidence')
            if gse == 'GSE158471' and 'genotype: Cnot4-WT' not in record['characteristics_ch1']:
                raise ValueError('Only WT source profiles are admitted')
            if gse == 'GSE158471' and not any(value.startswith('Stage: ' + stage) for value in record['characteristics_ch1']):
                raise ValueError('CNOT4 stage differs from source evidence')
            vector = vector[~vector.index.str.startswith('__')]
            if vector.index.has_duplicates or not np.isfinite(vector).all() or (vector < 0).any():
                raise ValueError('Invalid sorted expression')
            vectors.append(vector.rename(gsm))
            rows.append({'gse': gse, 'gsm': gsm, 'stage': stage, 'unit': unit,
                         'original_title': ' | '.join(record['title']),
                         'original_characteristics': ' | '.join(record['characteristics_ch1']),
                         'original_column': vector.name, 'tar_member': m.name})
    return pd.concat(vectors, axis=1), rows


def main():
    out = ROOT / 'results/data_expansion/sorted_validation'
    out.mkdir(parents=True, exist_ok=True)
    datasets = []
    selections = [('GSE158471', {'GSM4801338': 'LZ', 'GSM4801339': 'LZ', 'GSM4801340': 'PD', 'GSM4801341': 'PD'}, 'FPKM'),
                  ('GSE232415', {'GSM8618769': 'PS', 'GSM8618770': 'PS', 'GSM8618771': 'RS', 'GSM8618772': 'RS'}, 'counts')]
    for gse, selected, unit in selections:
        x, meta = tar_profiles(gse, selected, unit)
        datasets.append((gse, x, meta))
    x = pd.read_csv(ROOT / 'data/interim/expansion/H3mm15_sorted_symbol_counts.tsv.gz', sep='\t', index_col=0)
    stages = ['SG', 'preL', 'LZ', 'PD', 'Mll', 'RS']
    meta = []
    columns = []
    for j, stage in enumerate(stages):
        for rep in [1, 2]:
            gsm = f'GSM{8818498 + 2*j + rep-1}'
            record = source('GSE290629')['samples'][gsm]
            if 'genotype: WT' not in record['characteristics_ch1'] or f'rep {rep}' not in record['title'][0]:
                raise ValueError('H3mm15 WT source metadata mismatch')
            expected_type = {'SG': 'spermatogonia', 'preL': 'preleptotene', 'LZ': 'leptotene, zygotene',
                             'PD': 'pachytene, diplotene', 'Mll': 'meiosis II', 'RS': 'round spermatid'}[stage]
            if 'cell type: ' + expected_type not in record['characteristics_ch1']:
                raise ValueError('H3mm15 stage alias differs from source evidence')
            column = f'{stage}-{rep}P'
            columns.append(column)
            meta.append({'gse': 'GSE290629', 'gsm': gsm, 'stage': stage, 'unit': 'counts',
                         'original_title': ' | '.join(record['title']),
                         'original_characteristics': ' | '.join(record['characteristics_ch1']),
                         'original_column': column,
                         'mapping_note': 'Explicit source alias: stage + replicate + P; Mll corresponds to source meiosis II. No expression-based relabeling.'})
    x = x[columns].copy()
    x.columns = [r['gsm'] for r in meta]
    datasets.append(('GSE290629', x, meta))
    expected = {'SG': 'Spermatogonia', 'LZ': 'Meiotic', 'PD': 'Meiotic',
                'PS': 'Meiotic', 'Mll': 'Meiotic', 'RS': 'Postmeiotic'}
    # preL is a transition, not a forced binary SPG/meiotic truth label.
    weights, metadata = [], []
    for gse, x, meta in datasets:
        refs = {name: pd.read_csv(ROOT / f'results/stage_reference/pooled_{name}_reference.csv.gz', index_col=0)
                for name in ['original_4', 'centroid_4', 'centroid_13']}
        # Existing exported stage references cover a particular gene universe;
        # use exactly that pre-existing universe, do not refit reference genes.
        genes = refs['centroid_4'].index.intersection(x.index).sort_values()
        refs = {name: s.loc[genes].div(s.loc[genes].sum(), axis=1) * 1e6 for name, s in refs.items()}
        y = x.loc[genes].div(x.loc[genes].sum(), axis=1) * 1e6
        scale = 1 / np.sqrt(refs['centroid_4'].mean(axis=1).to_numpy() + 1)
        for row in meta:
            row['shared_genes'] = len(genes)
            row['matched_expression_fraction'] = float(x.loc[genes, row['gsm']].sum() / x[row['gsm']].sum())
            metadata.append(row)
            for name, s in refs.items():
                coef, _ = nnls(s.to_numpy() * scale[:, None], y[row['gsm']].to_numpy() * scale, maxiter=1000)
                p = pd.Series(coef / coef.sum(), index=s.columns)
                broad = (pd.Series({'Somatic': p['Somatic'], 'Spermatogonia': p['GC1'],
                    'Meiotic': p[[f'GC{i}' for i in range(2, 9)]].sum(),
                    'Postmeiotic': p[[f'GC{i}' for i in range(9, 13)]].sum()}) if name == 'centroid_13' else p)
                label = expected.get(row['stage'], 'transition_not_scored')
                weights.append({'gse': gse, 'gsm': row['gsm'], 'stage': row['stage'], 'model': name,
                    'expected_broad_class': label, 'largest_fitted_class': broad.idxmax(),
                    'shared_genes': len(genes), 'expected_class_weight': broad.get(label, np.nan),
                    **broad.to_dict()})
    frame = pd.DataFrame(weights)
    frame.to_csv(out / 'fitted_weights.csv', index=False)
    pd.DataFrame(metadata).to_csv(out / 'source_sample_evidence.csv', index=False)
    scored = frame[frame.expected_broad_class != 'transition_not_scored'].copy()
    scored['matches'] = scored.expected_broad_class == scored.largest_fitted_class
    summary = scored.groupby(['gse', 'model']).agg(profiles=('gsm', 'size'),
        largest_class_matches=('matches', 'sum'), median_expected_weight=('expected_class_weight', 'median'))
    summary.to_csv(out / 'summary.csv')
    notes = {'scope': 'Compatibility of 20 additional sorted WT expression profiles from three studies; 18 have a prespecified broad class, two preL transitions are unscored.',
        'fixed_before_fit': ['Only author-labeled WT populations; no KO profiles for validation.',
            'Previously generated pooled references and their existing gene universe; same weighting method and no retuning.',
            'FPKM and counts renormalized within the same shared genes; platform/genome build differences remain.'],
        'limitations': ['Reference expression profiles were previously assessed on another external dataset; these new profiles were not used for refinement.',
            'The reference gene universe was originally selected by availability in an earlier external table; this can limit generalization.',
            'Replicate labels do not prove unique donor animals or independent preparation across stages.',
            'A largest-class match is a diagnostic, not calibrated proportion accuracy or causal validation.',
            'CNOT4 data use mm9 symbols; H3mm15 gene mapping uses current org.Mm.eg.db.'],
        'reference_sha256': {name: hashlib.sha256((ROOT / f'results/stage_reference/pooled_{name}_reference.csv.gz').read_bytes()).hexdigest()
                             for name in ['original_4', 'centroid_4', 'centroid_13']}}
    (out / 'design.json').write_text(json.dumps(notes, indent=2) + '\n')
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5), layout='constrained')
    order = [r['gsm'] for r in metadata if r['stage'] in expected]
    labels = [r['gse'].removeprefix('GSE') + ':' + r['stage'] + ':' + r['gsm'][-2:]
              for r in metadata if r['stage'] in expected]
    for name, color in [('centroid_4', '#e0a458'), ('centroid_13', '#246b87')]:
        part = scored[scored.model == name].set_index('gsm').loc[order]
        ax.plot(range(len(order)), part.expected_class_weight, 'o', label=name, color=color)
    ax.set(xticks=range(len(labels)), xticklabels=labels, ylim=(0, 1),
           ylabel='Fitted expected-class weight\n(not cell proportions)',
           title='New sorted WT profiles: fixed reference compatibility')
    ax.tick_params(axis='x', rotation=75)
    ax.legend(frameon=False)
    fig.savefig(out / 'new_sorted_compatibility.png', dpi=180)
    print(summary.to_string())


if __name__ == '__main__':
    main()
