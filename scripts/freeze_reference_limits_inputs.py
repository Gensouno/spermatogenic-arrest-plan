#!/usr/bin/env python3
"""Freeze raw and pre-existing prepared inputs; never silently refresh a mismatch."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()


def main():
    rows={}
    for name in ['data/chen_umi_assets.json','data/leptotene_assets.json']:
        for r in json.loads((ROOT/name).read_text())['assets']:rows[r['path']]=r
    extra=['data/raw/ref/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.xlsx','results/reference_calibration/pooled_cpm.csv.gz',
        'data/ref/batches/counts.npz','data/ref/batches/batch_classes.csv','results/stage_reference/provenance.json',
        'data/ref/annotation/Mus_musculus.GRCm38.79.gtf.gz','data/ref/annotation/provenance.json','data/raw/labeled_reference/genes.tsv',
        'results/emptydrops_p15/type_profiles.csv.gz','emptydrops_p15_snapshot/profile_evidence.csv','leptotene_snapshot/sample_evidence.csv',
        'labeled_reference_snapshot/design.json']
    extra += [f'results/stage_reference/pooled_{m}_reference.csv.gz' for m in ['original_4','centroid_4','centroid_13']]
    for pmc in ['PMC6123400','PMC10541148']:
        path=f'data/raw/reference_limits/papers/{pmc}.xml';rows[path]={'path':path,'url':f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'}
    for path in extra:rows.setdefault(path,{'path':path,'provenance':'Previously cached project raw/prepared input; source reconstruction and prior manifest checks documented separately.'})
    for path,r in rows.items():
        p=ROOT/path
        if p.is_symlink() or not p.resolve().is_relative_to(ROOT):raise ValueError('Unsafe input path')
        sha=digest(p)
        if 'sha256' in r and (r['sha256']!=sha or r['bytes']!=p.stat().st_size):raise ValueError('Previously frozen source differs')
        r.update(bytes=p.stat().st_size,sha256=sha)
    value={'scope':'Raw source assets plus unchanged pre-existing prepared inputs. Newly generated Chen/technical/transfer results are deliberately excluded.',
        'assets':[rows[k] for k in sorted(rows)]}
    spec=ROOT/'data/reference_limits_inputs.json'
    if spec.exists() and json.loads(spec.read_text())!=value:raise ValueError('Frozen reference-limits inputs differ; not refreshed')
    if not spec.exists():spec.write_text(json.dumps(value,indent=2)+'\n')
    print('Frozen/verified',len(rows),'raw or pre-existing prepared inputs',flush=True)


if __name__=='__main__':main()
