"""Certify short escapes from the portrait envelope, not long-orbit shadowing."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time
import numpy as np
from product_binary_interval import parameters, model, mv


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    if args.mpmath_path: sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    import verify_product_coverings as full
    full.mp=mp;mp.mp.dps=mp.iv.dps=args.digits;iv=mp.iv
    exact=lambda v:iv.mpf(Fraction(v).numerator)/Fraction(v).denominator
    enc=lambda v:dict(decimal=str(v),binary=[list(t) for t in v._mpi_])
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    report=json.loads((args.source/'report.json').read_text(encoding='utf-8'))
    chosen=sorted(report['rows'],key=lambda r:r['minimum_radius'])[:args.examples]
    low=list(map(exact,report['threshold_lower']));high=list(map(exact,report['threshold_upper']))
    epsilon=exact(args.radius);rows=[];inputs=[]
    for row in chosen:
        filename=f'block-{row["block"]:03}-minimum.npz';source=args.source/filename
        with np.load(source) as f:
            states=f['states'];raw=dict(params=f['params'].tolist(),points=list(map(list,itertools.product(*f['levels'][:,:2]))));start=int(f['start']);minimum_at=int(f['minimum_at'])
        assert len(states)>args.steps and minimum_at==32
        par=parameters(raw['params'],np.array(raw['points']),iv);other=full.Model(raw)
        center=list(map(exact,states[0]));box=[c+iv.mpf([-epsilon.b,epsilon.b]) for c in center]
        assert all(x.a>lo.b and x.b<hi.a for x,lo,hi in zip(box,low,high))
        states_record=[list(map(enc,box))];radii=[];escapes=[];worst=0.
        for t in range(args.steps):
            if args.engine=='full':
                natural,J,r=other.value_jacobian(box)
                mid=[iv.mpf((mp.mpf(v.a)+mp.mpf(v.b))/2) for v in box]
                assert all(m.a>=v.a and m.b<=v.b for m,v in zip(mid,box))
                value,_=other.field(mid);change=mv(J,[v-m for v,m in zip(box,mid)])
                mean=[a+b for a,b in zip(value,change)]
                box=[iv.mpf([max(a.a,b.a),min(a.b,b.b)]) for a,b in zip(natural,mean)]
            else: box,J,r=model(box,par,iv)
            assert r.a>0;radii.append(enc(r));states_record.append(list(map(enc,box)))
            for j in range(3):
                if box[j].a>high[j].b: escapes.append(dict(step=t+1,axis=j,side='upper',margin=enc(box[j]-high[j])))
                if box[j].b<low[j].a: escapes.append(dict(step=t+1,axis=j,side='lower',margin=enc(low[j]-box[j])))
            worst=max(worst,max(abs(float((mp.mpf(v.a)+mp.mpf(v.b))/2)-states[t+1,j]) for j,v in enumerate(box)))
        assert escapes
        filename=f'block-{row["block"]:03}-certificate.json.gz'
        payload=dict(input_center=[float(v).hex() for v in states[0]],radius=args.radius,start=start,params=raw['params'],points=raw['points'],states=states_record,radii=radii,escapes=escapes)
        with gzip.open(args.out/filename,'wt',encoding='utf-8',newline='\n') as handle:json.dump(payload,handle,separators=(',',':'))
        result=dict(block=row['block'],start=start,minimum_offset=minimum_at,input_sha256=sha(source),certificate=filename,certificate_sha256=sha(args.out/filename),first_escape=escapes[0],largest_escape_margin=max(float(mp.mpf(e['margin']['decimal'].split(',')[0].lstrip('['))) for e in escapes),maximum_midpoint_difference_from_binary64=worst)
        rows.append(result);inputs.append(dict(file=source.name,sha256=sha(source)));(args.out/source.name).write_bytes(source.read_bytes())
        print(json.dumps({k:v for k,v in result.items() if k not in ['first_escape']}),flush=True)
    names=[Path(__file__).name,'product_binary_interval.py','verify_product_coverings.py']
    for name in names:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    result=dict(status='passed',engine=args.engine,digits=args.digits,mpmath=mp.__version__,steps=args.steps,initial_box_radius=args.radius,threshold_lower=report['threshold_lower'],threshold_upper=report['threshold_upper'],excursion_report_sha256=sha(args.source/'report.json'),selection='Three smallest normalization minima in the completed exploratory long-horizon audit.',inputs=inputs,cases=rows,seconds=time.time()-begin,scope='Every exact trajectory from each tiny initial box starts inside the declared portrait envelope and leaves it within forty updates. This disproves forward invariance of that envelope and certifies finite excursions. It does not shadow the full long orbit, identify an attractor, or prove a crisis mechanism.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out','mpmath-path']:p.add_argument('--'+name,type=Path,required=name!='mpmath-path')
    p.add_argument('--engine',choices=['factorized','full'],default='factorized');p.add_argument('--digits',type=int,default=140)
    p.add_argument('--radius',default='1e-80');p.add_argument('--steps',type=int,default=40);p.add_argument('--examples',type=int,default=3)
    main(p.parse_args())
