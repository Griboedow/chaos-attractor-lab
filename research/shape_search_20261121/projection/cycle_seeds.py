"""Numerical near-returns for interval validation of least-period cycles."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.spatial import cKDTree
from model import evaluate
from reproduce import raw_parameters


def orbit_eval(x,raw,p):
    beta,inv,tube,count=raw[0];n=int(count)
    b=raw[2:2+n,:3].copy()
    orbit=[];mat=np.eye(3)
    out=np.zeros(3);j=np.zeros((3,3));w=np.zeros(n)
    cov=np.zeros((3,3));dp=np.zeros((3,3))
    for _ in range(p):
        orbit.append(x.copy())
        evaluate(x,b,beta,inv,tube,out,j,w,cov,dp,int(raw[1,3]),raw[2,3])
        mat=j@mat
        x=out.copy()
    return x,mat,np.array(orbit)


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    for case in args.cases:
        dest=args.out/f'seeds-{case}.json'
        if dest.exists(): continue
        rec=json.loads((args.source/f'case-{case}.json').read_text())
        with np.load(rec['cloud']) as f:z=f['points'].astype(float)
        raw=raw_parameters(rec)
        tree=cKDTree(z);rad=np.sqrt(np.mean(np.sum((z-z.mean(0))**2,axis=1)))
        found=[];orbits=[]
        for p in range(1,args.max_period+1):
            def f(x):
                y,j,_=orbit_eval(x,raw,p)
                return y-x,j-np.eye(3)
            errs=np.linalg.norm(z[p:]-z[:-p],axis=1)
            ids=np.argpartition(errs,min(1000,len(errs)-1))[:1000]
            ids=ids[np.argsort(errs[ids])]
            picked=[]
            for j in ids:
                if all(np.linalg.norm(z[j]-z[k])>.075*rad for k in picked):
                    picked.append(j)
                    if len(picked)>=args.seeds: break
            count=0
            for j in picked:
                try:
                    sol=root(f,z[j],jac=True,tol=1e-11)
                    x=sol.x; y,mat,orb=orbit_eval(x,raw,p)
                    err=float(np.max(abs(y-x)))
                    if not np.isfinite(err) or err>1e-9: continue
                    if any(np.linalg.norm(orb[d]-x)<1e-6 for d in range(1,p) if p%d==0): continue
                    if any(len(old)==p and np.min(np.linalg.norm(old-x,axis=1))<1e-5 for old in orbits): continue
                    near=float(tree.query(orb)[0].max()/rad)
                    if near>args.distance: continue
                    mult=np.linalg.eigvals(mat)
                except (ValueError,FloatingPointError,np.linalg.LinAlgError):
                    continue
                found.append(dict(period=p,seed=x.tolist(),orbit=orb.tolist(),residual=err,max_orbit_distance_to_sample_rms=near,multipliers=[[float(a.real),float(a.imag)] for a in mult],numerical_unstable=int(np.sum(abs(mult)>1))))
                orbits.append(orb)
                count+=1
                if count>=args.per_period: break
        tmp=dest.with_suffix('.tmp.json')
        tmp.write_text(json.dumps(dict(case=case,cycles=found),indent=2)+'\n',encoding='utf8')
        os.replace(tmp,dest)
        print(case,len(found),sorted(set(r['numerical_unstable'] for r in found)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--cases',type=int,nargs='+',required=True)
    p.add_argument('--max-period',type=int,default=12)
    p.add_argument('--seeds',type=int,default=10)
    p.add_argument('--per-period',type=int,default=3)
    p.add_argument('--distance',type=float,default=.07)
    main(p.parse_args())
