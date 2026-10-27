#!/usr/bin/env python3
"""Symmetric post-hoc coverage of every fixed candidate in six earlier contrasts."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]


def unique_effects(frame,symbol_column):
    if not np.isfinite(frame.logFC).all():raise ValueError('Nonfinite external effect')
    symbols=frame[symbol_column]
    keep=symbols.notna() & symbols.astype(str).str.len().gt(0)
    mapped=frame.loc[keep].copy()
    mapped=mapped[~mapped[symbol_column].duplicated(keep=False)]
    return mapped.set_index(symbol_column).logFC


def cached_loo(root,spec):
    if spec['cohort'].startswith('Lrh1'):
        table=pd.read_csv(root/'results/followup/leave_one_out_effects.csv.gz')
        mapping=pd.read_csv(root/'results/followup/Lrh1_mapping.csv')
        keys=mapping.input_id.str.replace(r'\.[0-9]+$','',regex=True)
        if keys.duplicated().any():raise ValueError('Ambiguous Lrh1 identifier map')
        table=table[table.cohort==spec['cohort']].copy()
        table['symbol']=table.gene.map(dict(zip(keys,mapping.symbol)))
        expected=8;symbol='symbol'
    else:
        table=pd.read_csv(root/'results/data_expansion/leave_one_out_logFC.csv.gz')
        table=table[table.cohort==spec['cohort']].copy();expected=6;symbol='gene'
    if table.deleted.nunique()!=expected or table.duplicated(['deleted','gene']).any():raise ValueError('External omission inventory differs')
    return {str(v):unique_effects(part,symbol) for v,part in table.groupby('deleted')},expected


def transfer_summary(effect,omissions,genes,directions):
    """Missing effects/variants never become zero or pass stability."""
    full=effect.reindex(genes);loo=pd.DataFrame({v:e.reindex(genes) for v,e in omissions.items()})
    signed=full.to_numpy()*directions
    oriented=loo.to_numpy()*np.asarray(directions)[:,None]
    complete=np.isfinite(oriented).all(axis=1);available=full.notna().to_numpy()
    return pd.DataFrame({'gene':genes,'tested_exact_unique_symbol':available,'external_logFC':full.to_numpy(),
        'candidate_directional_external_log2FC':signed,'sign_matches':available & (signed>0),
        'directional_1p5fold':available & (signed>=np.log2(1.5)),
        'available_LOO_variants':np.isfinite(oriented).sum(axis=1),'expected_LOO_variants':len(omissions),
        'all_LOO_available':complete,'baseline_and_all_LOO_sign_match':available & (signed>0) & complete & (oriented>0).all(axis=1),
        'baseline_and_all_LOO_1p5fold':available & (signed>=np.log2(1.5)) & complete & (oriented>=np.log2(1.5)).all(axis=1)})


def main():
    contract=json.loads((ROOT/'data/replicability_context_contract.json').read_text())
    out=ROOT/'results/replicability';candidate=pd.read_csv(out/'candidate_stage_context.csv')
    all_genes=pd.read_csv(out/'primary_gene_summary.csv.gz').gene
    rows=[];summary=[];coverage=[];fingerprints={}
    for spec in contract['contrasts']:
        effect=unique_effects(pd.read_csv(ROOT/spec['source']),spec['symbol_column'])
        omissions,expected=cached_loo(ROOT,spec)
        part=transfer_summary(effect,omissions,candidate.gene,np.where(candidate.descriptive_direction=='up',1,-1))
        part['cohort']=spec['cohort'];part['model']=spec['model'];part['age']=spec['age'];part['role']=spec['role']
        part['candidate_direction']=candidate.descriptive_direction.to_numpy()
        part['primary_LOO_priority']=candidate.LOO_stable_at_least_3_units.to_numpy()
        rows.append(part)
        for subset,frame in [('all_candidates',part),('primary_LOO_priorities',part[part.primary_LOO_priority])]:
            for direction in ['all','up','down']:
                d=frame if direction=='all' else frame[frame.candidate_direction==direction]
                summary.append({'cohort':spec['cohort'],'model':spec['model'],'age':spec['age'],'role':spec['role'],
                    'subset':subset,'candidate_direction':direction,'fixed_genes':len(d),
                    'tested_exact_unique_symbol':int(d.tested_exact_unique_symbol.sum()),
                    'all_LOO_available':int(d.all_LOO_available.sum()),
                    **{k:int(d[k].sum()) for k in ['sign_matches','directional_1p5fold',
                        'baseline_and_all_LOO_sign_match','baseline_and_all_LOO_1p5fold']}})
        primary=effect.reindex(all_genes)
        coverage.append({'cohort':spec['cohort'],'source_tested_unique_symbols':len(effect),'primary_common_genes':len(all_genes),
            'primary_covered':int(primary.notna().sum()),'primary_not_covered':int(primary.isna().sum()),
            'primary_up_1p5fold':int((primary>=np.log2(1.5)).sum()),'primary_down_1p5fold':int((primary<=-np.log2(1.5)).sum()),
            'expected_LOO_variants':expected})
        fingerprints[spec['source']]=hashlib.sha256((ROOT/spec['source']).read_bytes()).hexdigest()
    for name in ['results/followup/leave_one_out_effects.csv.gz','results/followup/Lrh1_mapping.csv',
                 'results/data_expansion/leave_one_out_logFC.csv.gz']:
        fingerprints[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    pd.concat(rows,ignore_index=True).to_csv(out/'external_candidate_context.csv',index=False)
    pd.DataFrame(summary).to_csv(out/'external_context_summary.csv',index=False)
    pd.DataFrame(coverage).to_csv(out/'external_context_coverage.csv',index=False)
    (out/'external_context_design.json').write_text(json.dumps({'status':contract['analysis_status'],
        'contract_sha256':hashlib.sha256((ROOT/'data/replicability_context_contract.json').read_bytes()).hexdigest(),
        'source_sha256':fingerprints,'candidate_genes':len(candidate),'contrasts':6,'additional_independent_primary_votes':0,
        'scope':'Frozen previously fitted effects and LOO caches only; no external refit or new inferential test in this branch.',
        'limitations':contract['limitations']},indent=2)+'\n')
    print(pd.DataFrame(summary).query("subset=='primary_LOO_priorities' and candidate_direction=='all'").to_string(index=False),flush=True)


if __name__=='__main__':main()
