#!/usr/bin/env python3
"""Offline source audit and exploratory categorical-age Setx analysis."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AGES = ['P10', 'P12', 'P14', 'P16', 'P18']
GROUPS = ['Control', 'KO']
CELLS = [f'{a}_{g}' for a in AGES for g in GROUPS]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def check_inputs(root):
    for asset in json.loads((root / 'data/setx_interaction_inputs.json').read_text())['assets']:
        p = Path(asset['path'])
        if p.is_absolute() or '..' in p.parts or not (root / p).resolve().is_relative_to(root.resolve()):
            raise ValueError('Unsafe frozen input path')
        path = root / p
        if path.stat().st_size != asset['bytes'] or digest(path) != asset['sha256']:
            raise ValueError('Frozen input changed: ' + str(p))


def parse_member(name):
    pattern = (r'(GSM\d+)_d(10|12|14|16|18)_(\d+)_(WT|Mut)_GES14_(\d+)_'
               r'([ACGT]{6})_GRCm38\.trimmed\.genes\.results\.gz')
    m = re.fullmatch(pattern, name)
    if not m:
        raise ValueError('Unexpected Setx filename: ' + name)
    gsm, age, code, genotype, serial, barcode = m.groups()
    return dict(gsm=gsm, age='P'+age, group='Control' if genotype == 'WT' else 'KO',
                source_code=code, filename_genotype=genotype, library_serial=serial,
                barcode_token=barcode, member=name)


def audit_design(meta):
    if meta.gsm.duplicated().any() or len(meta) != 40:
        raise ValueError('Unexpected sample identity/count')
    labels = meta.age + '_' + meta.group
    if labels.value_counts().to_dict() != {c: 4 for c in CELLS}:
        raise ValueError('Unexpected age/genotype design')
    cell_design = np.column_stack([(labels == c).astype(float) for c in CELLS])
    barcodes = sorted(meta.barcode_token.unique())
    barcode_design = np.column_stack([(meta.barcode_token == b).astype(float) for b in barcodes])
    nesting = meta.assign(cell=labels).groupby('barcode_token').cell.nunique()
    repeated = meta[meta.source_code.duplicated(False)][['gsm', 'age', 'group', 'source_code']]
    return dict(samples=len(meta), cells=len(CELLS), unique_source_codes=int(meta.source_code.nunique()),
                repeated_source_code_rows=repeated.to_dict('records'),
                unique_barcode_tokens=len(barcodes), barcodes_nested_in_cells=bool((nesting == 1).all()),
                cell_design_rank=int(np.linalg.matrix_rank(cell_design)), cell_design_columns=10,
                cell_plus_barcode_rank=int(np.linalg.matrix_rank(np.column_stack([cell_design, barcode_design]))),
                cell_plus_barcode_columns=10 + len(barcodes),
                interpretation='Filename source-code meaning unconfirmed; barcode tokens are not proven batches. '
                'A model allowing arbitrary barcode effects cannot separately identify age/genotype cell effects. '
                'No animal pairing, relabeling, or technical-artifact diagnosis inferred.')


def prepare(root, out):
    check_inputs(root)
    meta = pd.read_csv(root / 'new_models_snapshot/sample_evidence.csv')
    meta = meta[meta.gse == 'GSE224875'].sort_values('gsm').reset_index(drop=True)
    archive_audit = pd.read_csv(root / 'new_models_snapshot/Setx_archive_audit.csv').set_index('gsm')
    source = json.loads((root / 'results/replicability_discovery/GSE224875_metadata.json').read_text())
    columns, design = {}, []
    genes = None
    with tarfile.open(root / 'data/raw/replicability_discovery/GSE224875_RAW.tar') as archive:
        for record in meta.to_dict('records'):
            gsm = record['gsm']
            entry = archive_audit.loc[gsm]
            parsed = parse_member(entry.member)
            if parsed['gsm'] != gsm or parsed['age'] != record['age'] or parsed['group'] != record['group']:
                raise ValueError('Filename versus frozen sample annotation disagreement')
            title = source['samples'][gsm]['title'][0]
            wanted = f"{'wt' if record['group']=='Control' else 'SetxKO'} testis, {record['age'][1:]}dpp, rep"
            if not title.startswith(wanted) or title != record['source_title']:
                raise ValueError('Source title versus frozen annotation disagreement')
            payload = archive.extractfile(entry.member).read()
            if hashlib.sha256(payload).hexdigest() != entry.member_sha256:
                raise ValueError('Frozen genome-member digest differs')
            frame = pd.read_csv(io.BytesIO(gzip.decompress(payload)), sep='\t')
            if frame.gene_id.isna().any() or frame.gene_id.duplicated().any():
                raise ValueError('Invalid gene IDs')
            if genes is None:
                genes = frame.gene_id.copy()
            if not genes.equals(frame.gene_id):
                raise ValueError('Gene ordering differs')
            values = frame.TPM.to_numpy(float)
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError('Invalid TPM')
            columns[gsm] = values
            design.append(parsed)
    matrix = pd.DataFrame(columns, index=pd.Index(genes, name='gene'))
    matrix.to_csv(out / 'TPM.tsv.gz', sep='\t', compression={'method':'gzip','mtime':0})
    design = pd.DataFrame(design)
    design.to_csv(out / 'sample_design.csv', index=False)
    audit = audit_design(design)
    audit.update(source_genes=len(genes), source_members_hash_matched=len(columns),
                 TPM_reconstructed='40 original GRCm38 gene RSEM members read; no cached matrix, piRNA or rounded counts used.')
    (out / 'design_audit.json').write_text(json.dumps(audit, indent=2) + '\n')


def summarize_candidates(effects, candidates, stable):
    if effects.duplicated(['variant','symbol']).any():
        raise ValueError('Ambiguous symbol reached candidate summary')
    required = ['main', 'pseudo_0p1', 'pseudo_1', 'exclude_ambiguous_pair']
    omissions = sorted(v for v in effects.variant.unique() if v.startswith('omit_'))
    if len(omissions) != 40:
        raise ValueError('Incomplete omission suite')
    required += omissions
    rows = []
    for c in candidates.to_dict('records'):
        symbol = c['gene']
        sub = effects[effects.symbol == symbol].set_index('variant')
        complete = set(required) <= set(sub.index)
        delta = sub.delta_P18
        main = sub.loc['main'] if 'main' in sub.index else None
        rows.append(dict(symbol=symbol, fixed_stable_priority=symbol in stable,
            earlier_direction=c['descriptive_direction'], variants_measured=len(sub), required_variants=len(required),
            all_variants_measured=complete,
            main_delta_P18=float(main.delta_P18) if main is not None else np.nan,
            main_conditional_BY=float(main.BY) if main is not None else np.nan,
            main_conditional_BH=float(main.BH) if main is not None else np.nan,
            minimum_delta_P18=float(delta.min()) if len(delta) else np.nan,
            maximum_delta_P18=float(delta.max()) if len(delta) else np.nan,
            all_variants_negative=bool(complete and (delta < 0).all()),
            all_variants_down_1p5fold=bool(complete and (delta <= -np.log2(1.5)).all())))
    return pd.DataFrame(rows)


def summarize(root, out):
    variants = pd.read_csv(out / 'variant_summary.csv')
    effects = pd.read_csv(out / 'candidate_effects.csv.gz')
    candidates = pd.read_csv(root / 'replicability_snapshot/candidate_stage_context.csv')
    stable = set(pd.read_csv(root / 'replicability_snapshot/stable_descriptive_priorities.csv').gene)
    summary = summarize_candidates(effects, candidates, stable)
    summary.to_csv(out / 'candidate_stability.csv', index=False)
    priority = summary[summary.fixed_stable_priority]
    overlap = pd.read_csv(out / 'all_gene_statistics.csv.gz')
    sets = {v:set(d.loc[d.BY < .05, 'gene']) for v,d in overlap.groupby('variant')}
    report = dict(scope='Exploratory post-hoc conditional library-level interaction; not biological confirmation.',
        fixed_candidates=len(candidates), fixed_stable_priorities=len(stable),
        main_tested_genes=int(variants.loc[variants.variant=='main','tested_genes'].iloc[0]),
        main_conditional_BY05=int(variants.loc[variants.variant=='main','conditional_BY05'].iloc[0]),
        all_four_variant_BY05_overlap=len(set.intersection(*sets.values())),
        fixed_priorities_complete=int(priority.all_variants_measured.sum()),
        fixed_priorities_negative_all_variants=int(priority.all_variants_negative.sum()),
        fixed_priorities_down_1p5fold_all_variants=int(priority.all_variants_down_1p5fold.sum()),
        fixed_priorities_main_conditional_BY05=int((priority.main_conditional_BY<.05).sum()),
        all_variant_count=44, omission_count=40,
        interpretation='The complete-gene family is adjusted before fixed candidates are annotated. '
        'P18-P10 deltas compare KO-Control differences; a negative delta is not itself a KO-Control decrease at either age. '
        'No subgroup FDR guarantee, animal-independence verification, barcode-effect separation or shared-mechanism claim.')
    (out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    order = sorted(stable)
    main = effects[effects.variant=='main'].set_index('symbol').reindex(order)
    fig, axes = plt.subplots(1,2,figsize=(10,11),gridspec_kw={'width_ratios':[2,1]})
    image = axes[0].imshow(main[[f'effect_{a}' for a in AGES]],aspect='auto',cmap='coolwarm',vmin=-3,vmax=3)
    axes[0].set_yticks(range(len(order)),order,fontsize=7);axes[0].set_xticks(range(5),AGES)
    axes[0].set_title('KO minus WT log2(TPM + 0.5) means',fontsize=10)
    p = priority.set_index('symbol').reindex(order)
    for i,row in enumerate(p.itertuples()):
        if np.isfinite(row.main_delta_P18):
            axes[1].plot([row.minimum_delta_P18,row.maximum_delta_P18],[i,i],color='0.45',lw=1)
            axes[1].plot(row.main_delta_P18,i,'o',color='tab:blue',ms=3)
    axes[1].set_ylim(len(order)-.5,-.5);axes[1].set_yticks([]);axes[1].axvline(0,color='0.4',lw=.8)
    axes[1].axvline(-np.log2(1.5),color='tab:red',ls=':',lw=.8)
    axes[1].set_title('P18 minus P10 interaction effect',fontsize=10)
    axes[1].set_xlabel('Full fit and sensitivity range\n44 variants; range is not a confidence interval',fontsize=8)
    fig.colorbar(image,ax=axes[0],shrink=.35,label='Log2 difference')
    fig.suptitle('Setx: fixed 36 descriptive priorities\nExploratory source-library fits; biological independence unconfirmed',fontsize=11)
    fig.subplots_adjust(left=.14,right=.98,top=.92,bottom=.08,wspace=.2)
    fig.savefig(out/'fixed_priorities_interaction.png',dpi=160);plt.close(fig)
    print(json.dumps(report,indent=2),flush=True)


def main():
    out = ROOT/'results/setx_interaction';out.mkdir(parents=True,exist_ok=True)
    prepare(ROOT,out)
    subprocess.run(['bash','envs/rscript.sh','code/analyze_setx_interaction.R'],cwd=ROOT,check=True)
    summarize(ROOT,out)


if __name__ == '__main__':
    main()
