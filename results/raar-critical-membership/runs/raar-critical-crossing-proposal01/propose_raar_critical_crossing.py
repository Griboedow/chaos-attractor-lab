"""Numerical determinant-sign scout on images of certified covering exit faces.

Sampling guides a later interval calculation; it proves no critical point or
membership in an invariant set. Every sampled time and face is retained.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from numba import njit
from search_raar_circles import step


@njit(cache=True)
def follow(states, params, steps):
    x=states.copy()
    determinants=np.empty((steps+1,len(x)))
    axes=np.empty_like(determinants)
    saved=np.empty((steps+1,len(x),3))
    for k in range(steps+1):
        saved[k]=x
        for i in range(len(x)):
            y,J,axis,flag=step(x[i],params)
            assert flag==0
            determinants[k,i]=np.linalg.det(J)
            axes[k,i]=axis
            x[i]=y
    return saved,determinants,axes


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    raw=json.loads(args.input.read_text())
    params=np.array(raw['params'])
    records=[]
    for w,(c,B,r) in enumerate(zip(raw['centers'],raw['bases'],raw['radii'])):
        grid=np.linspace(-1,1,args.grid)
        q=np.array([[u,v,z] for u in [-1,1] for v in grid for z in grid])
        states=np.array(c)+(q*np.array(r))@np.array(B).T
        saved,ds,axes=follow(states,params,args.steps)
        ranges=np.stack([ds[:,:args.grid**2].min(1),ds[:,:args.grid**2].max(1),ds[:,args.grid**2:].min(1),ds[:,args.grid**2:].max(1)],axis=1)
        candidates=np.flatnonzero(((ranges[:,1]<0)&(ranges[:,2]>0))|((ranges[:,3]<0)&(ranges[:,0]>0)))
        np.savez_compressed(args.out/f'window-{w}.npz',states=saved,determinants=ds,axes=axes,normalized_starts=q,ranges=ranges)
        records.append({'window':w,'candidate_steps':candidates.tolist(),'candidate_ranges':ranges[candidates].tolist(),'minimum_sampled_axis':float(axes.min())})
    for source in [Path(__file__),Path(__file__).with_name('search_raar_circles.py')]:
        (args.out/source.name).write_bytes(source.read_bytes())
    (args.out/'input.json').write_bytes(args.input.read_bytes())
    report={'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'grid':args.grid,'steps':args.steps,'windows':records,'scope':'Numerical finite face sampling only; no interval sign or invariant-set membership certificate.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',type=Path,default=Path('results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'))
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--grid',type=int,default=9)
    p.add_argument('--steps',type=int,default=70)
    main(p.parse_args())
