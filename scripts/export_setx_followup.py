#!/usr/bin/env python3
"""Freeze compact source audit and count-model sensitivity, preserving all prior branches."""
import json
from pathlib import Path
import shutil
from check_setx_followup import SOURCE_TABLES,SOURCE_JSON,LENGTH_TABLES,check
from export_reference_limits import local_source_closure
from audit_setx_sources import digest

ROOT=Path(__file__).resolve().parents[1]


def main():
    check(ROOT)
    log=json.loads((ROOT/'data/raw/setx_source_audit/retrieval.json').read_text())
    for item in log:
        if item['name']=='preprint_metadata.json':item['response_validation']='Rejected: HTTP success but API returned non-numeric value error; not used as scientific evidence.'
    audit=dict(scope='Targeted lookup; URLs contain no credentials. Failed downloads and invalid responses are retained explicitly.',
        requests=len(log),http_failures=sum('error_type' in r for r in log),logical_response_rejections=1,records=log)
    (ROOT/'results/setx_source_audit/retrieval_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    sources=local_source_closure(ROOT,['scripts/audit_setx_sources.py','scripts/run_setx_length.py',
        'scripts/check_setx_followup.py','scripts/export_setx_followup.py','tests/test_setx_followup.py'])+['code/analyze_setx_length.R','code/analyze_setx_interaction.R']
    specs=['data/setx_source_audit_contract.json','data/setx_source_audit_inputs.json','data/setx_length_contract.json',
        'data/setx_length_inputs.json','envs/python.sh','envs/rscript.sh','envs/network.sh','docs/SETX_SOURCE_AND_LENGTH.md']
    for branch,files in [('setx_source_audit',SOURCE_TABLES+SOURCE_JSON+['retrieval_audit.json','reproduction.json']),
                         ('setx_length',LENGTH_TABLES+['summary.json','reproduction.json'])]:
        snap=ROOT/(branch+'_snapshot')
        if (snap/'manifest.json').exists():raise ValueError('Already frozen; create a new named snapshot for changes')
        snap.mkdir(exist_ok=True)
        for name in files:shutil.copy2(ROOT/'results'/branch/name,snap/name)
        manifest=dict(scope='Source reconciliation and conditional model sensitivity; not biological validation.',
            source_sha256={p:digest(ROOT/p) for p in sources},input_specification_sha256={p:digest(ROOT/p) for p in specs},
            files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(snap.iterdir()) if p.is_file()])
        (snap/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(branch,len(files),'files;',len(sources)+len(specs),'source/spec checks')


if __name__=='__main__':main()
