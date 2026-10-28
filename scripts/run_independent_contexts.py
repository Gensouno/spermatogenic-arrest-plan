#!/usr/bin/env python3
"""Audit cached GEO assets, prepare gene counts and fit separate exploratory contexts."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COHORTS = ['Rif1_P14', 'Hormad_dsfb_P19', 'Mael_P16', 'Morc1_P14p5']


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def validate_counts(frame):
    if frame.index.has_duplicates or frame.index.isna().any():
        raise ValueError('Missing/duplicate gene identity')
    x = frame.to_numpy(float)
    if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
        raise ValueError('Counts must be finite nonnegative source integers; never round')
    if (x.sum(axis=0) <= 0).any(): raise ValueError('Empty library')


def gene_rows(frame):
    """Exclude repeat loci by deposited repeat annotations, not symbol heuristics."""
    return frame[['repName', 'repClass', 'repFamily']].isna().all(axis=1)


def prepare(root, out):
    spec = json.loads((root/'data/independent_contexts_inputs.json').read_text())
    for a in spec['assets']:
        p = root/a['path']
        if p.stat().st_size != a['bytes'] or digest(p) != a['sha256']:
            raise ValueError('Input identity differs: '+a['path'])
    samples=[]; qc=[]
    definitions=[('Rif1_P14','GSE294151_rif1_count.txt.gz','Symbol',
        ['con-1','con-2','con-3','scKO-1','scKO-2','scKO-3'],
        ['GSM8899502','GSM8899503','GSM8899504','GSM8899505','GSM8899506','GSM8899507']),
        ('Hormad_dsfb_P19','GSE274491_read_counts_all_control_dsfb.txt.gz','name',
        ['control_1','control_2','control_3','dsfb_1','dsfb_2','dsfb_3'],
        ['GSM8450234','GSM8450236','GSM8450238','GSM8450240','GSM8450242','GSM8450244']),
        ('Mael_P16','GSE222415_counts_rna_ribo_annotated.txt.gz','geneID',
        ['geneCount_mrna_het_'+str(i) for i in range(1,4)]+['geneCount_mrna_mael_'+str(i) for i in range(1,4)],
        ['GSM6923077','GSM6923078','GSM6923079','GSM6923080','GSM6923081','GSM6923082'])]
    for cohort, name, gene, columns, gsms in definitions:
        selected=[]; total=0
        for chunk in pd.read_csv(root/'data/raw/independent_discovery'/name,sep='\t',chunksize=200000):
            total+=len(chunk)
            if cohort=='Mael_P16': chunk=chunk.loc[gene_rows(chunk)]
            selected.append(chunk[[gene]+columns])
        frame=pd.concat(selected,ignore_index=True).set_index(gene)
        validate_counts(frame)
        frame.index.name='gene'
        frame.to_csv(out/(cohort+'.tsv.gz'),sep='\t',compression={'method':'gzip','mtime':0})
        qc.append(dict(cohort=cohort,source_rows=total,gene_rows=len(frame),excluded_repeat_rows=total-len(frame),
                       selected_RNA_libraries=6,noninteger_values=0,duplicate_gene_ids=0))
        for i,(column,gsm) in enumerate(zip(columns,gsms)):
            samples.append(dict(cohort=cohort,matrix_sample=column,gsm_alias=gsm,
                alias_scope='tentative group/suffix alias; not explicit source crosswalk' if cohort=='Rif1_P14' else 'exact GEO title/description',
                group='Control' if i<3 else 'Perturbation',library_count_sum=int(frame[column].sum())))
    morc=[];source_rows=None;gene_index=None
    for a in spec['assets']:
        if a.get('gsm') not in ['GSM'+str(i) for i in range(1547662,1547668)]:continue
        d=pd.read_csv(root/a['path'],sep='\t',header=None)
        source_rows=len(d);d=d[d[0].str.fullmatch(r'ENSMUSG\d+')].set_index(0)
        if gene_index is None:gene_index=d.index
        if not gene_index.equals(d.index):raise ValueError('Morc1 gene rows/order differ')
        morc.append(d[6].rename(a['gsm']))
    if len(morc)!=6:raise ValueError('Morc1 six file identities required')
    frame=pd.concat(morc,axis=1);validate_counts(frame);frame.index.name='gene'
    frame.to_csv(out/'Morc1_P14p5.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
    qc.append(dict(cohort='Morc1_P14p5',source_rows=source_rows,gene_rows=len(frame),excluded_repeat_rows=source_rows-len(frame),selected_RNA_libraries=6,noninteger_values=0,duplicate_gene_ids=0))
    for i,column in enumerate(frame.columns):
        samples.append(dict(cohort='Morc1_P14p5',matrix_sample=column,gsm_alias=column,alias_scope='exact GSM-specific source file',group='Control' if i<3 else 'Perturbation',library_count_sum=int(frame[column].sum())))
    pd.DataFrame(samples).to_csv(out/'samples.csv',index=False)
    pd.DataFrame(qc).to_csv(out/'input_qc.csv',index=False)


def summarize(root,out):
    effects=pd.read_csv(out/'priority_effects.csv.gz')
    priorities=pd.read_csv(root/'replicability_snapshot/stable_descriptive_priorities.csv').gene.tolist()
    rows=[]
    for cohort in COHORTS:
        e=effects[effects.cohort==cohort]
        required={'full'}|{'omit_'+s for s in pd.read_csv(out/'samples.csv').query('cohort == @cohort').matrix_sample}
        for symbol in priorities:
            d=e[e.symbol==symbol]
            if d.variant.duplicated().any(): raise ValueError('Duplicate symbol/fit')
            complete=set(d.variant)==required and np.isfinite(d.logFC).all()
            main=d[d.variant=='full']
            rows.append(dict(cohort=cohort,symbol=symbol,variants_measured=len(d),variants_required=7,
                complete=bool(complete),full_logFC=float(main.logFC.iloc[0]) if len(main) else np.nan,
                minimum_logFC=d.logFC.min(),maximum_logFC=d.logFC.max(),
                all_negative=bool(complete and (d.logFC<0).all()),
                all_down_1p5fold=bool(complete and (d.logFC<=-np.log2(1.5)).all())))
    stable=pd.DataFrame(rows);stable.to_csv(out/'priority_stability.csv',index=False)
    summary=[]
    for cohort,d in stable.groupby('cohort',sort=False):
        summary.append(dict(cohort=cohort,fixed_priorities=36,full_measured=int(d.full_logFC.notna().sum()),
            full_negative=int((d.full_logFC<0).sum()),full_down_1p5fold=int((d.full_logFC<=-np.log2(1.5)).sum()),
            complete_all_seven=int(d.complete.sum()),all_negative=int(d.all_negative.sum()),
            all_down_1p5fold=int(d.all_down_1p5fold.sum())))
    pd.DataFrame(summary).to_csv(out/'priority_summary.csv',index=False)
    common=stable.groupby('symbol').all_down_1p5fold.agg(['sum','count'])
    result=dict(contexts=4,libraries=24,fits=28,fixed_priorities=36,
        all_contexts_all_omissions_down_1p5fold=int(((common['sum']==4)&(common['count']==4)).sum()),
        primary_added_votes=0,scope='Exploratory separate-source contexts, no expanded primary test or subgroup FDR control.')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default='results/independent_contexts');args=p.parse_args()
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    prepare(ROOT,out)
    subprocess.run(['bash','envs/rscript.sh','code/analyze_independent_contexts.R',str(out)],cwd=ROOT,check=True)
    summarize(ROOT,out)


if __name__=='__main__': main()
