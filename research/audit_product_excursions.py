"""Rational certificate audit and deterministic replay of short excursions."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from audit_product_orientation import interval, point, add, sub, mul, isum, contains


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8'))
def unpack(p):
    with gzip.open(p,'rt',encoding='utf-8') as f:return json.load(f)


def singular(folder):
    r=read(folder/'report.json');assert r['status']=='passed'
    assert sha(folder/'input.json')==r['input_sha256']
    boxes=list(map(interval,r['box']));f=list(map(interval,r['field_at_center']))
    J=[[interval(v) for v in row] for row in r['derivative']]
    R=[[interval(v) for v in row] for row in r['remainder']]
    C=list(map(point,r['inverse_diagonal']));x=list(map(point,r['center']));eps=Q(r['radius'])
    for i in range(3):
        assert contains(boxes[i],add(x[i],(-eps,eps)))
        for j in range(3):assert contains(R[i][j],sub(point(int(i==j)),mul(C[i],J[i][j])))
        exact=add(sub(x[i],mul(C[i],f[i])),isum(mul(R[i][j],(-eps,eps)) for j in range(3)))
        k=interval(r['krawczyk'][i]);assert contains(k,exact)
        assert boxes[i][0]<k[0]<=k[1]<boxes[i][1]
        assert J[i][i][0]-sum(max(abs(v) for v in J[i][j]) for j in range(3) if j!=i)>0
    norm=max(sum(max(abs(v) for v in cell) for cell in row) for row in R)
    assert interval(r['contraction'])[1]>=norm and norm<1
    return dict(run=folder.name,report_sha256=sha(folder/'report.json'),contraction_upper=float(norm))


def excursion(folder):
    r=read(folder/'report.json');assert r['status']=='passed'
    lo=list(map(Q,r['threshold_lower']));hi=list(map(Q,r['threshold_upper']));rows=[];payloads=[]
    for row in r['cases']:
        p=folder/row['certificate'];assert sha(p)==row['certificate_sha256'];v=unpack(p)
        eps=Q(v['radius']);ss=[[interval(x) for x in s] for s in v['states']]
        assert len(ss)==r['steps']+1 and len(v['radii'])==r['steps']
        centers=[Q(float.fromhex(x)) for x in v['input_center']]
        for j in range(3):
            assert contains(ss[0][j],(centers[j]-eps,centers[j]+eps))
            assert lo[j]<ss[0][j][0]<=ss[0][j][1]<hi[j]
        assert all(interval(x)[0]>0 for x in v['radii'])
        assert v['escapes']
        for e in v['escapes']:
            b=ss[e['step']][e['axis']]
            margin=sub(b,point(hi[e['axis']])) if e['side']=='upper' else sub(point(lo[e['axis']]),b)
            assert margin[0]>0 and contains(interval(e['margin']),margin)
        assert row['first_escape']==v['escapes'][0]
        rows.append(dict(block=row['block'],first_escape=row['first_escape']['step'],largest_margin=row['largest_escape_margin'],minimum_normalization_lower=min(float(interval(x)[0]) for x in v['radii'])))
        payloads.append(v)
    return dict(run=folder.name,report_sha256=sha(folder/'report.json'),cases=rows),payloads


def main(a):
    a.out.mkdir(parents=True,exist_ok=False)
    s=[singular(p) for p in a.singularities]
    checks=[excursion(p) for p in a.excursions]
    for u,v in zip(checks[0][1],checks[1][1]):
        assert u['input_center']==v['input_center'] and u['params']==v['params'] and u['points']==v['points']
        for x,y in zip(u['states'],v['states']):
            for xx,yy in zip(x,y):
                aa,bb=interval(xx),interval(yy);assert max(aa[0],bb[0])<=min(aa[1],bb[1])
    replay=[]
    for folder in a.excursions:
        r=read(folder/'report.json');dest=a.replay/folder.name
        cmd=[sys.executable,str(folder/'certify_product_excursions.py'),'--source',str(a.source),'--out',str(dest),'--mpmath-path',str(a.mpmath_path),'--engine',r['engine'],'--digits',str(r['digits']),'--radius',r['initial_box_radius'],'--steps',str(r['steps'])]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        for row in r['cases']:assert unpack(folder/row['certificate'])==unpack(dest/row['certificate'])
        replay.append(dict(run=folder.name,payloads=len(r['cases']),exact=True))
    result=dict(status='passed',singularities=s,excursions=[x[0] for x in checks],cross_engine_state_intersections=True,replays=replay,scope='Exact rational audit of supplied interval arithmetic and escape inequalities; six deterministic payload replays. Transcendental enclosures depend on mpmath interval arithmetic, checked through two formulas sharing that library. No long-orbit shadowing assertion.')
    (a.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'audit_product_orientation.py']:(a.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out','replay','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    for name in ['singularities','excursions']:p.add_argument('--'+name,type=Path,nargs=2,required=True)
    main(p.parse_args())
