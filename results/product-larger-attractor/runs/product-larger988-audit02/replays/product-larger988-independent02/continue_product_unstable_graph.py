"""Validated finite continuation of the exact unstable point used in the escape.

Every state is enclosed in an exact rational moving chart. Distances to a saved
numerical state are finite geometric observations, not invariant membership.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.spatial import cKDTree
from product_binary_interval import parameters,model,inverse,mm,mv


def main(args):
    begin=time.time();args.out.mkdir(parents=True,exist_ok=False)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    save=lambda p,v:p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'product_binary_interval.py','verify_product_coverings.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    sys.path.append(str(args.mpmath_path));import mpmath as mp
    import verify_product_coverings as full
    mp.mp.dps=mp.iv.dps=args.digits;iv=mp.iv;full.mp=mp
    exact=lambda s,ctx=iv:ctx.mpf(Fraction(s).numerator)/Fraction(s).denominator
    enc=lambda v:dict(decimal=str(v),binary=[list(t) for t in v._mpi_])
    def decode(v):
        ends=[]
        for sign,m,e,_ in v['binary']:ends.append(exact(Fraction((-1 if sign else 1)*m)*Fraction(2)**e))
        return iv.mpf([ends[0].a,ends[1].b])
    cv=lambda a,ctx=iv:[[exact(v,ctx) for v in row] for row in a]
    sym=lambda v:iv.mpf([-abs(v).b,abs(v).b])
    raw=json.loads(args.input.read_text());prior=json.loads(args.certificate.read_text())
    assert prior['status']=='passed' and prior['windows_sha256']==sha(args.input)
    par=parameters(raw['params'],np.array(raw['points']),mp.mp);pari=parameters(raw['params'],np.array(raw['points']),iv)
    independent=full.Model(raw)
    field=independent.value_jacobian if args.engine=='full' else lambda x:model(x,pari,iv)
    E=[[Fraction(s) for s in row] for row in prior['graph']['E']];Ei=inverse(E)
    EI,II,EM=cv(E),cv(Ei),cv(E,mp.mp)
    p=[exact(s) for s in prior['root']['center_decimal']];pm=[exact(s,mp.mp) for s in prior['root']['center_decimal']]
    K=list(map(decode,prior['root']['krawczyk']))
    eps=exact(prior['epsilon']);L=exact(prior['graph']['L'])
    nominal=[pm[i]+EM[i][0]*exact(prior['epsilon'],mp.mp) for i in range(3)]
    # Preserve the correlation of the two unknown stable coordinates. Mapping
    # their box through E and E^-1 before starting would lose this information.
    correction=mv(II,[K[i]-p[i]+p[i]+EI[i][0]*eps-iv.mpf(nominal[i]) for i in range(3)])
    delta=[correction[0],correction[1]+sym(L*eps),correction[2]+sym(L*eps)]
    with np.load(args.source) as f:reference=f['states']
    tree=cKDTree(reference);trail=[];failure=None
    for n in range(0,args.steps+1,4):
        off=mv(EI,delta);box=[iv.mpf(nominal[i])+off[i] for i in range(3)]
        dist,idx=tree.query(np.array(list(map(float,nominal))))
        nearest=[exact(float(v)) for v in reference[idx]]
        distance=iv.sqrt(sum((box[i]-nearest[i])**2 for i in range(3)))
        row=dict(iteration=n,center=[mp.nstr(v,args.digits+10) for v in nominal],
                 center_binary=[list(v._mpf_) for v in nominal],E=[[str(v) for v in r] for r in E],
                 inverse=[[str(v) for v in r] for r in Ei],chart_error=list(map(enc,delta)),box=list(map(enc,box)),
                 nearest_index=int(idx),nearest_state=reference[idx].tolist(),distance=enc(distance))
        trail.append(row)
        if n==args.steps:break
        try:
            outer=box;J=[[iv.mpf(i==j) for j in range(3)] for i in range(3)]
            next_nominal=nominal;y0=[iv.mpf(v) for v in nominal]
            Jnum=[[mp.mpf(i==j) for j in range(3)] for i in range(3)];radii=[]
            for phase in range(4):
                natural,jac,rr=field(outer);assert rr.a>0;radii.append(rr)
                J=mm(jac,J)
                mid=[iv.mpf((mp.mpf(v.a)+mp.mpf(v.b))/2) for v in outer]
                base,_,_=field(mid);change=mv(jac,[outer[i]-mid[i] for i in range(3)])
                mean=[base[i]+change[i] for i in range(3)]
                outer=[iv.mpf([max(a.a,b.a),min(a.b,b.b)]) for a,b in zip(natural,mean)]
                y0,_,_=field(y0)
                next_nominal,jnum,_=model(next_nominal,par,mp.mp);Jnum=mm(jnum,Jnum)
            # Keep the return-map eigenchart during the local departure. A QR
            # chart is useful only after that departure; using it immediately
            # destroys the near-diagonal stable-to-unstable return coupling.
            if n<80:En=E;Ein=Ei
            else:
                next_frame,_=np.linalg.qr(np.array(Jnum,dtype=float)@np.array(E,dtype=float))
                En=[[Fraction(float(v)) for v in r] for r in next_frame];Ein=inverse(En)
            EN,IN=cv(En),cv(Ein)
            center_error=mv(IN,[y0[i]-iv.mpf(next_nominal[i]) for i in range(3)])
            transition=mm(mm(IN,J),EI);prop=mv(transition,delta)
            newdelta=[center_error[i]+prop[i] for i in range(3)]
            row['step']=dict(radii=list(map(enc,radii)),derivative=[[enc(v) for v in r] for r in J],
                center_error=list(map(enc,center_error)),transition=[[enc(v) for v in r] for r in transition])
            if max(float(abs(v).b) for v in newdelta)>.01:
                failure=dict(iteration=n,reason='Declared interval-width stop (0.01 in chart coordinates); no orbit escape inferred.');break
            nominal=next_nominal;E,Ei,EI,II=En,Ein,EN,IN;delta=newdelta
        except (AssertionError,ValueError,ZeroDivisionError) as error:
            failure=dict(iteration=n,reason=type(error).__name__+': '+str(error));break
        if n%20==16:save(args.out/'progress.json',dict(iteration=n+4,seconds=time.time()-begin))
    save(args.out/'trajectory.json',trail)
    report=dict(status='completed' if failure is None else 'finite_prefix_only',engine=args.engine,digits=args.digits,
        requested_steps=args.steps,completed_steps=trail[-1]['iteration'],stop=failure,
        sources={str(p):sha(p) for p in [args.input,args.certificate,args.source]},
        trajectory_sha256=sha(args.out/'trajectory.json'),seconds=time.time()-begin,
        checkpoints=[dict(iteration=r['iteration'],distance=r['distance'],box=r['box'],nearest_index=r['nearest_index'])
                     for r in trail if r['iteration']>=80 and r['iteration']%4==0],
        scope='Validated finite prefix of the exact unstable-graph point from the prior certificate, with its true unknown stable coordinates enclosed. Distance to a saved numerical state proves neither recurrence to an invariant set nor attraction.')
    save(args.out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='checkpoints'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','certificate','source','mpmath-path','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--engine',choices=['factorized','full'],default='factorized');p.add_argument('--digits',type=int,default=120)
    p.add_argument('--steps',type=int,default=200);main(p.parse_args())
