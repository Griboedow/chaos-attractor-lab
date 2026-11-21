"""Freeze the period-doubling and dense-polygon audits as compact portable records."""
import argparse
import gzip
import json
import shutil
from pathlib import Path


def read(path):
    return json.loads(path.read_text())


def write(path,obj):
    data=json.dumps(obj,indent=None if path.suffix=='.gz' else 2,allow_nan=False).encode() + b'\n'
    path.write_bytes(gzip.compress(data,mtime=0) if path.suffix=='.gz' else data)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--packet',type=Path,required=True)
    ap.add_argument('--isotopies',type=Path)
    ap.add_argument('--braids',type=Path)
    ap.add_argument('--branch',type=Path,required=True)
    ap.add_argument('--cable',type=Path,required=True)
    ap.add_argument('--doublings',type=Path,nargs='+')
    ap.add_argument('--fine-links',type=Path,nargs='+')
    ap.add_argument('--long-return',type=Path)
    ap.add_argument('--jones',type=Path)
    a=ap.parse_args();r=a.packet
    values=read(a.branch/'report.json')
    assert values['trials'][-1]['fraction']==1 and len(values['trials'])==36
    write(r/'period-doubling.json',values)
    intermediate=r/'doubling-witnesses';intermediate.mkdir(exist_ok=True)
    link_evidence=r/'doubling-link-witnesses';link_evidence.mkdir(exist_ok=True)
    first=next(t for t in values['trials'] if t.get('witness_file'))
    shutil.copyfile(a.branch/first['witness_file'],intermediate/'case-first-daughter.json')
    shutil.copyfile(a.cable/'link-witness.json',link_evidence/'case-first-daughter.json')
    write(r/'doubling-cable.json',read(a.cable/'report.json'))
    if a.doublings:
        pairs=[]
        for folder in a.doublings:
            for p in read(folder/'report.json'):
                pairs.append(p);order=p['child_period_order']
                shutil.copyfile(folder/f'link-{order}.json',link_evidence/f'case-order-{order}.json')
        if {p['child_period_order'] for p in pairs}!={2,4,8,16}:raise ValueError('Incomplete doubled orbit collection')
        write(r/'doubling-pairs.json',pairs)
    if a.fine_links:
        fine=[p for folder in a.fine_links for p in read(folder/'report.json')]
        assert {x['order'] for x in fine}=={2,4,8,16}
        write(r/'doubling-link-refinement.json',fine)
    if a.long_return:
        long=read(a.long_return/'report.json')
        assert len(long['initial_condition_runs'])==3
        write(r/'long-doubling.json',long)
    if a.isotopies:
        dense={f.stem[5:]:read(f) for f in sorted(a.isotopies.glob('case-*.json'))}
        for name,rec in dense.items():
            if len(rec['views'])<2 or rec['passes'][-1][-1]!=len(rec['final_vertices']):raise ValueError(name)
        failures={f.stem[7:]:read(f)['reason'] for f in sorted(a.isotopies.glob('failed-*.json'))}
        write(r/'dense-polygon-isotopies.json.gz',dict(success=dense,unresolved=failures))
    if a.braids:
        recs={f.stem[5:]:read(f) for f in sorted(a.braids.glob('case-*.json'))}
        write(r/'two-strand-witnesses.json',recs)
    if a.jones:
        with gzip.open(a.jones,'rt') as stream:polys=json.load(stream)
        if len(polys['resolved'])+len(polys['unresolved'])!=174:raise ValueError('Incomplete Jones audit')
        shutil.copyfile(a.jones,r/'jones-polynomials.json.gz')


if __name__=='__main__':
    main()
