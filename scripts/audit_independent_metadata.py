#!/usr/bin/env python3
"""Retain relevant source assertions while excluding submitter contact fields."""
import argparse
import gzip
import json
from pathlib import Path
from run_independent_contexts import ROOT,digest

SERIES={'Series_geo_accession','Series_title','Series_summary','Series_overall_design','Series_pubmed_id','Series_supplementary_file'}
SAMPLE={'Sample_geo_accession','Sample_title','Sample_source_name_ch1','Sample_organism_ch1','Sample_characteristics_ch1','Sample_extract_protocol_ch1','Sample_description','Sample_data_processing','Sample_library_strategy','Sample_supplementary_file_1'}


def parse(text):
    series={};samples={};current=None
    for line in text.splitlines():
        if line.startswith('^SERIES = '): current=series
        elif line.startswith('^SAMPLE = '): current={};samples[line.split(' = ',1)[1]]=current
        elif line.startswith('^'):current=None
        elif current is not None and line.startswith('!') and ' = ' in line:
            key,value=line[1:].split(' = ',1)
            if key in SERIES|SAMPLE:current.setdefault(key,[]).append(value)
    return series,samples


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default=str(ROOT/'results/independent_contexts'));args=p.parse_args();out=Path(args.output)
    records=[]
    for a in json.loads((ROOT/'data/independent_contexts_inputs.json').read_text())['assets']:
        if not a['path'].endswith('_family.soft.gz'):continue
        path=ROOT/a['path']
        if digest(path)!=a['sha256']:raise ValueError('SOFT identity differs')
        series,samples=parse(gzip.decompress(path.read_bytes()).decode())
        records.append(dict(source_path=a['path'],sha256=a['sha256'],series=series,samples=samples))
    if len(records)!=9:raise ValueError('Expected nine selected source metadata records')
    (out/'source_metadata.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':main()
