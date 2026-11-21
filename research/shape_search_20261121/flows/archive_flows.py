"""Compare all saved ODE geometries with the campaign archive and its ODE subset."""
import argparse
import csv
import json
import os
import sys
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from compare import compare,write


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',required=True,type=Path)
    ap.add_argument('--cache',required=True,type=Path,help='Campaign cache root, containing the flows directory')
    ap.add_argument('--reference',required=True,type=Path,help='Full archived reference fingerprints')
    ap.add_argument('--old',type=Path,default=Path('results/attractor-novelty/reference'))
    a=ap.parse_args();m=json.loads((a.packet/'manifest.json').read_text())
    index={}
    for coll,tag in [('cases','flows'),('linked_cycles','links'),('knotted_periodic_skeletons','skeletons')]:
        for case in m[coll]:index[f"{tag}/{case['case']}"]=case['cloud_relative']
    with (a.packet/'periodic-catalog.csv').open() as stream:
        for case in csv.DictReader(stream):index['cycles/'+case['case']]=case['cloud_relative']
    index_file=a.packet/'archive-inputs.json';write(index_file,index)
    compare(index_file,a.reference,a.cache,a.packet/'archive-comparison.json.gz',False)
    names=set()
    for file in ('runs.json','atlas-runs.json','nearby-runs.json'):
        for row in json.loads((a.old/file).read_text()):
            if row.get('status','ok')=='ok' and row.get('success',True) and 'map' not in row:
                names.add('reference/'+row['name'])
    rows=[row for row in json.loads(a.reference.read_text())['rows'] if row['name'] in names]
    ode_input=a.cache/'flows'/'ode-reference-fingerprints.json'
    write(ode_input,dict(reference_count=len(rows),rows=rows))
    compare(index_file,ode_input,a.cache,a.packet/'ode-reference-comparison.json.gz',False)


if __name__=='__main__':
    main()
