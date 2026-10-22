#!/usr/bin/env python3
"""Cross-fit composition-only predictions of matched log mean CPM contrasts.

Descriptive model check only: held-out genes remain biologically correlated;
weights do not identify a causal contribution or cell-intrinsic regulation.
"""
import itertools
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/composition_prediction'
    out.mkdir(parents=True,exist_ok=True)
    cohorts=pd.read_csv(ROOT/'results/import_qc/count_cohorts.csv').gse.tolist()+['GSE101511']
    results, metrics, coefficients=[],[],[]
    for refname in ['pooled','ST_only']:
        reference=pd.read_csv(ROOT/f'results/reference_calibration/{refname}_cpm.csv.gz',index_col=0).T
        for gse in cohorts:
            x=pd.read_csv(ROOT/f'data/bulk/{gse}_matrix.tsv.gz',sep='\t',index_col=0)
            meta=pd.read_csv(ROOT/f'data/bulk/{gse}_samples.tsv',sep='\t').set_index('matrix_sample').loc[x.columns]
            cpm=x.div(x.sum(axis=0),axis=1)*1e6
            tested=pd.read_csv(ROOT/f'results/de/{gse}_de.csv').gene
            genes=reference.index.intersection(tested).intersection(x.index)
            s=reference.loc[genes].to_numpy()
            y=cpm.loc[genes].to_numpy()
            weights=1/np.sqrt(s.mean(axis=1)+1)
            assignment=np.empty(len(genes),dtype=int)
            rng=np.random.default_rng(42)
            assignment[rng.permutation(len(genes))]=np.arange(len(genes))%5
            prediction=np.full_like(y,np.nan)
            for fold in range(5):
                tr,te=assignment!=fold,assignment==fold
                for j,sample in enumerate(x.columns):
                    coef,_=nnls(s[tr]*weights[tr,None],y[tr,j]*weights[tr])
                    p=coef/coef.sum()
                    prediction[te,j]=s[te]@p
                    row={'gse':gse,'reference':refname,'fold':fold,'sample':sample,'coefficient_sum':coef.sum()}
                    row.update(dict(zip(reference.columns,p)))
                    coefficients.append(row)
            if not np.isfinite(prediction).all():
                raise ValueError('Incomplete cross-fit predictions')
            wt=(meta.group=='WT').to_numpy()
            ko=(meta.group=='KO').to_numpy()
            observed=np.log2(y[:,ko].mean(axis=1)+.5)-np.log2(y[:,wt].mean(axis=1)+.5)
            predicted=np.log2(prediction[:,ko].mean(axis=1)+.5)-np.log2(prediction[:,wt].mean(axis=1)+.5)
            residual=observed-predicted
            frame=pd.DataFrame({'gse':gse,'reference':refname,'gene':genes,'fold':assignment,
                    'observed_log_mean_CPM_contrast':observed,'predicted_composition_contrast':predicted,
                    'residual_contrast':residual})
            results.append(frame)
            metrics.append({'gse':gse,'reference':refname,'genes':len(genes),
                    'spearman':spearmanr(observed,predicted).statistic,
                    'rmse':np.sqrt(np.mean(residual**2)),
                    'zero_prediction_rmse':np.sqrt(np.mean(observed**2)),
                    'mae':np.abs(residual).mean()})
    all_results=pd.concat(results,ignore_index=True)
    all_results.to_csv(out/'crossfit_gene_contrasts.csv.gz',index=False)
    pd.DataFrame(coefficients).to_csv(out/'crossfit_RNA_weights.csv',index=False)
    pd.DataFrame(metrics).to_csv(out/'prediction_metrics.csv',index=False)
    # Exact covariance accounting includes both cross terms. No ratio of these
    # terms is labelled a percentage of biological causation.
    decompositions=[]
    for refname in ['pooled','ST_only']:
        sub=all_results[all_results.reference==refname]
        o=sub.pivot(index='gene',columns='gse',values='observed_log_mean_CPM_contrast').dropna()
        p=sub.pivot(index='gene',columns='gse',values='predicted_composition_contrast').loc[o.index,o.columns]
        r=o-p
        for a,b in itertools.combinations(o.columns,2):
            cov=lambda x,y: float(np.cov(x,y,ddof=1)[0,1])
            total=cov(o[a],o[b]); cp=cov(p[a],p[b]); cr=cov(r[a],r[b])
            cross1=cov(p[a],r[b]); cross2=cov(r[a],p[b])
            assert np.isclose(total,cp+cr+cross1+cross2)
            decompositions.append({'reference':refname,'cohort_a':a,'cohort_b':b,'genes':len(o),
                 'observed_cov':total,'predicted_cov':cp,'residual_cov':cr,
                 'predictedA_residualB_cov':cross1,'residualA_predictedB_cov':cross2})
    pd.DataFrame(decompositions).to_csv(out/'covariance_terms.csv',index=False)
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for ax,gse in zip(axes.flat,cohorts):
        sub=all_results[(all_results.reference=='pooled')&(all_results.gse==gse)]
        xx=sub.predicted_composition_contrast; yy=sub.observed_log_mean_CPM_contrast
        ax.hexbin(xx,yy,gridsize=40,mincnt=1,bins='log',cmap='viridis')
        lim=max(np.quantile(np.abs(xx),.99),np.quantile(np.abs(yy),.99),.2)
        ax.plot([-lim,lim],[-lim,lim],color='red',lw=.8)
        ax.set(title=gse,xlabel='Predicted composition contrast',ylabel='Observed log mean CPM contrast')
    axes.flat[-1].axis('off')
    fig.savefig(out/'heldout_prediction.png',dpi=170)
    plt.close(fig)
    (out/'design.json').write_text(json.dumps({'folds':5,'seed':42,'pseudocount_CPM':.5,
          'references':['pooled','ST_only'],'cohorts':cohorts,
          'note':'Observed contrasts use log of group-mean CPM to match composition predictions, not limma logFC. Gene-level cross-fitting is not independent biological validation. No inferential p-values or causal shares.'},indent=2)+'\n')
    print(pd.DataFrame(metrics).to_string(index=False))


if __name__=='__main__':
    main()
