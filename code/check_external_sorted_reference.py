#!/usr/bin/env python3
"""Independent sorted-cell profile check; not a full-mixture calibration.

da Cruz et al. 2016 have a mixed 2C population, LZ and pachytene meiotic cells,
and round spermatids. 2C is not a pure spermatogonial or somatic reference.
"""
import hashlib
import json
import os
from pathlib import Path
import urllib.request

import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT=Path(__file__).resolve().parents[1]
URL='https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12864-016-2618-1/MediaObjects/12864_2016_2618_MOESM2_ESM.xlsx'


def main():
    path=ROOT/'data/raw/external_reference/12864_2016_2618_MOESM2_ESM.xlsx'
    path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        with urllib.request.urlopen(URL,timeout=25) as response:
            data=response.read()
        if not data.startswith(b'PK'):
            raise ValueError('Not an XLSX archive')
        path.write_bytes(data)
    source=pd.read_excel(path,sheet_name='Hoja1').set_index('Feature ID')
    source=source[['2C RPKM','LZ RPKM','PS RPKM','RS RPKM']]
    out=ROOT/'results/external_sorted_reference'
    out.mkdir(parents=True,exist_ok=True)
    source.to_csv(out/'source_RPKM_as_published.csv')
    numeric=source.apply(pd.to_numeric,errors='coerce')
    malformed=numeric.isna().any(axis=1)
    source.loc[malformed].to_csv(out/'excluded_nonnumeric_rows.csv')
    # Published strings such as '1,160,51' are ambiguous. Do not guess values
    # for influential genes (including Prm1/Prm2/Sycp3); record their exclusion.
    source=numeric.loc[~malformed].copy()
    if source.index.has_duplicates or source.isna().any().any() or (source<0).any().any():
        raise ValueError('Invalid sorted-cell abundance table')
    source.to_csv(out/'source_RPKM.csv')
    expectations={'LZ RPKM':'Meiotic','PS RPKM':'Meiotic','RS RPKM':'Postmeiotic'}
    weights, scores, contrasts=[],[],[]
    for refname in ['pooled','ST_only']:
        reference=pd.read_csv(ROOT/f'results/reference_calibration/{refname}_cpm.csv.gz',index_col=0).T
        common=reference.index.intersection(source.index)
        s=reference.loc[common].to_numpy()
        s=s/s.sum(axis=0)*1e6
        # RPKM -> relative abundance over the measured common-gene universe.
        y=source.loc[common].to_numpy()
        y=y/y.sum(axis=0)*1e6
        scale=1/np.sqrt(s.mean(axis=1)+1)
        for j,stage in enumerate(source.columns):
            coef,_=nnls(s*scale[:,None],y[:,j]*scale)
            p=coef/coef.sum()
            row={'reference':refname,'stage':stage,'shared_genes':len(common),
                 'expected_class':expectations.get(stage,'mixed_2C_not_scored'),
                 'largest_fitted_class':reference.columns[np.argmax(p)]}
            row.update(dict(zip(reference.columns,p)))
            weights.append(row)
            if stage in expectations:
                scores.append({'reference':refname,'stage':stage,
                    'expected_class':expectations[stage],
                    'expected_class_RNA_weight':p[list(reference.columns).index(expectations[stage])]})
        # Pure LZ and PS are BOTH meiotic. Their changes expose within-class
        # state differences that one average Meiotic profile cannot represent.
        effect=np.log2(source['PS RPKM']+.5)-np.log2(source['LZ RPKM']+.5)
        contrasts.append({'reference':refname,'shared_genes':len(common),
            'pure_meiotic_PS_vs_LZ_median_abs_log2_RPKM':float(effect.loc[common].abs().median()),
            'genes_abs_log2_RPKM_gt1':int((effect.loc[common].abs()>1).sum())})
    frame=pd.DataFrame(weights)
    frame.to_csv(out/'fitted_RNA_weights.csv',index=False)
    pd.DataFrame(scores).to_csv(out/'expected_class_weights.csv',index=False)
    pd.DataFrame(contrasts).to_csv(out/'within_meiotic_stage_contrast.csv',index=False)
    provenance={'url':URL,'doi':'10.1186/s12864-016-2618-1','sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'source_genes':len(source),'excluded_nonnumeric_genes':int(malformed.sum()),
      'normalization':'Renormalize each source and reference profile over shared genes to 1e6.',
      'scope':'External sorted-profile compatibility check, not known proportion mixture recovery or cell-number validation.',
      'limitations':['2C includes spermatogonia and somatic cells; never used as pure cell-class truth.',
        'Single expression profile per stage; no replicate-based uncertainty.',
        'Publication table contains a selected expressed-gene universe and rounded RPKM.',
        'Seven genes with ambiguous number strings excluded, including influential stage markers; this can affect compatibility.',
        'FACS-RNAseq and Drop-seq have assay and cell-state differences.',
        'No reannotation of legacy gene symbols; intersection is explicitly reported.',
        'Fit uses all common genes; expected-class weight is a compatibility diagnostic, not a certified proportion error.']}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    plot=frame[frame.reference=='pooled'].set_index('stage')[list(reference.columns)]
    plot.plot.bar(stacked=True,ax=ax,color=['#246b87','#e0a458','#8a5a83','#6b8e52'])
    ax.set(ylabel='Fitted RNA weight over common genes',xlabel=None,
           title='Independent sorted profiles: broad reference compatibility')
    ax.tick_params(axis='x',rotation=0)
    ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False)
    fig.savefig(out/'sorted_profile_compatibility.png',dpi=180)
    plt.close(fig)
    print(frame.to_string(index=False))
    print(pd.DataFrame(contrasts).to_string(index=False))


if __name__=='__main__':
    main()
