#!/usr/bin/env python3
"""Freeze prepared expression, original fits/cache and independent stage resources."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
COHORTS=['GSE236388','GSE283904','GSE68283','GSE68284','GSE101511','GSE294209','GSE38096','GSE42601']


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def input_paths():
    paths=[]
    for g in COHORTS:
        paths.append(f'results/de/{g}_de.csv')
        base='results/arrays' if g in COHORTS[-2:] else 'data/bulk'
        suffix='_log2_expression.tsv.gz' if base=='results/arrays' else '_matrix.tsv.gz'
        paths.extend([f'{base}/{g}{suffix}',f'{base}/{g}_samples.tsv'])
    paths+=['results/robustness/leave_one_out_logFC.csv.gz','results/all_cohort_comparison/common_logFC.csv.gz',
        'results/reference_limits/chen_stage_profiles.csv.gz','results/reference_limits/chen_profile_evidence.csv',
        'results/reference_limits/chen_preparation_audit.json','reference_limits_snapshot/transfer_design.json',
        'data/raw/ref/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx',
        'results_snapshot/stage_reference/provenance.json']
    context=json.loads((ROOT/'data/replicability_context_contract.json').read_text())
    paths+=[r['source'] for r in context['contrasts']]
    paths+=['results/data_expansion/leave_one_out_logFC.csv.gz','results/followup/leave_one_out_effects.csv.gz',
            'results/followup/Lrh1_mapping.csv']
    return sorted(set(paths))


def main():
    assets=[{'path':p,'bytes':(ROOT/p).stat().st_size,'sha256':digest(ROOT/p)} for p in input_paths()]
    spec={'scope':'Prepared expression re-fit and 40 primary count omissions; frozen historical results are comparators. Six external contexts use frozen prior fits/LOO without re-fitting. No raw preprocessing or independent biological validation.', 'assets':assets}
    (ROOT/'data/replicability_inputs.json').write_text(json.dumps(spec,indent=2)+'\n')
    print('Frozen',len(assets),'input files')


if __name__=='__main__':main()
