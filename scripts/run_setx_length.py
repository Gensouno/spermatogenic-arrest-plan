#!/usr/bin/env python3
"""Reconstruct RSEM estimates/lengths from frozen source TAR and fit conditional QL models."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import numpy as np
import pandas as pd
from check_project import verify_manifest
from run_setx_interaction import check_inputs

ROOT=Path(__file__).resolve().parents[1]


def main():
    check_inputs(ROOT);verify_manifest(ROOT,'setx_interaction_snapshot')
    for asset in json.loads((ROOT/'data/setx_length_inputs.json').read_text())['assets']:
        p=ROOT/asset['path']
        if p.stat().st_size!=asset['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=asset['sha256']:
            raise ValueError('Frozen length-sensitivity input differs: '+asset['path'])
    out=ROOT/'results/setx_length';out.mkdir(parents=True,exist_ok=True)
    inventory=pd.read_csv(ROOT/'new_models_snapshot/Setx_archive_audit.csv').sort_values('gsm')
    columns={k:{} for k in ['expected_count','effective_length']};genes=None
    with tarfile.open(ROOT/'data/raw/replicability_discovery/GSE224875_RAW.tar') as tar:
        for record in inventory.itertuples():
            payload=tar.extractfile(record.member).read()
            if hashlib.sha256(payload).hexdigest()!=record.member_sha256:raise ValueError('Raw genome member differs')
            d=pd.read_csv(io.BytesIO(gzip.decompress(payload)),sep='\t')
            if genes is None:genes=d.gene_id.copy()
            if not genes.equals(d.gene_id):raise ValueError('Original gene order differs')
            for k in columns:
                x=d[k].to_numpy(float)
                if not np.isfinite(x).all() or (x<0).any():raise ValueError('Invalid RSEM '+k)
                columns[k][record.gsm]=x
    for k,values in columns.items():
        pd.DataFrame(values,index=pd.Index(genes,name='gene')).to_csv(out/(k+'.tsv.gz'),sep='\t',compression={'method':'gzip','mtime':0})
    subprocess.run(['bash','envs/rscript.sh','code/analyze_setx_length.R'],cwd=ROOT,check=True)
    d=pd.read_csv(out/'all_gene_statistics.csv.gz');old=pd.read_csv(ROOT/'setx_interaction_snapshot/all_gene_statistics.csv.gz')
    original=old[old.variant=='main'][['gene','delta_P18','BY']].rename(columns={'delta_P18':'TPM_delta_P18','BY':'TPM_BY'})
    comparisons=[]
    for variant,part in d.groupby('variant',sort=True):
        matched=part.merge(original,on='gene',how='inner')
        comparisons.append(dict(variant=variant,matched_original_gene_IDs=len(matched),
            TPM_conditional_BY05_in_matched=int((matched.TPM_BY<.05).sum()),
            count_conditional_BY05_in_matched=int((matched.BY<.05).sum()),
            both_conditional_BY05_in_matched=int(((matched.BY<.05)&(matched.TPM_BY<.05)).sum())))
    pd.DataFrame(comparisons).to_csv(out/'matched_universe_comparison.csv',index=False)
    candidates=pd.read_csv(ROOT/'replicability_snapshot/candidate_stage_context.csv')
    stable=set(pd.read_csv(ROOT/'replicability_snapshot/stable_descriptive_priorities.csv').gene)
    context=candidates[['gene','descriptive_direction']].rename(columns={'gene':'symbol'})
    subsets=[];summaries=[]
    for variant,part in d.groupby('variant',sort=True):
        mapped=part[part.symbol.notna()]
        if mapped.symbol.duplicated().any():raise ValueError('Ambiguous symbol in length context')
        table=context.merge(mapped,on='symbol',how='left').merge(original,on='gene',how='left')
        table['variant']=variant;table['fixed_stable_priority']=table.symbol.isin(stable)
        subsets.append(table)
        p=table[table.fixed_stable_priority];comp=p[p.delta_P18.notna()&p.TPM_delta_P18.notna()]
        summary=dict(variant=variant,fixed_priorities=len(p),measured=int(p.gene.notna().sum()),matched_to_TPM=len(comp),
            negative_delta=int((p.delta_P18<0).sum()),negative_delta_1p5fold=int((p.delta_P18<=-np.log2(1.5)).sum()),
            same_delta_sign_as_TPM=int(((comp.delta_P18<0)==(comp.TPM_delta_P18<0)).sum()),
            conditional_BY05=int((p.BY<.05).sum()))
        summaries.append(summary)
    pd.concat(subsets,ignore_index=True).to_csv(out/'candidate_context.csv',index=False)
    pd.DataFrame(summaries).to_csv(out/'candidate_summary.csv',index=False)
    report=dict(scope='Same-source post-hoc length-aware estimated-count sensitivity, not independent biological replication.',
        original_gene_rows=len(genes),source_genome_members=40,variants=2,
        estimated_counts_rounded=False,quantification_inferential_replicates_available=False,
        candidate_summary=summaries,
        limitations=['No count-based single-sample omission suite in this branch.','Count and TPM filtering/effect estimands differ.','Animal independence, pooling and barcode/source limitations persist.','No new primary genetic-unit vote, subgroup FDR guarantee or causal mechanism inference.'])
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(pd.read_csv(out/'variant_summary.csv').to_string(index=False),flush=True)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
