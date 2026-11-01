"""Finite branch and reflection diagnostics of the retained RK45 clouds."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from numba import njit
from search_rk45_oscillator import oscillator_step


@njit(cache=True)
def branches(states, params):
    n=len(states); rejects=np.empty(n,dtype=np.int64); determinants=np.empty(n)
    images=np.empty_like(states); residual=0.
    for i in range(n):
        s=states[i].copy();s[2]=np.log(s[2])
        new,J,r,code=oscillator_step(s,params[1],params[0]);assert code==0
        reflected=s.copy();reflected[:2]*=-1
        other,_,_,code=oscillator_step(reflected,params[1],params[0]);assert code==0
        other[:2]*=-1
        residual=max(residual,np.max(np.abs(other-new)))
        rejects[i]=r;determinants[i]=np.linalg.det(J)
        new[2]=np.exp(new[2]);images[i]=new
    return rejects,determinants,images,residual


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[]
    for ident in [335,448,1890,508,2006,726]:
        source=args.gallery/f'case-{ident}.npz'
        with np.load(source) as f:states=f['states'];params=f['params']
        rejects,determinants,images,residual=branches(states,params)
        transitions=np.zeros((2,2),dtype=np.int64)
        sign=states[:,1]>0
        np.add.at(transitions,(sign[:-1].astype(int),sign[1:].astype(int)),1)
        uq,nq=np.unique(rejects,return_counts=True)
        row=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),points=len(states),
            reflection_max_residual=float(residual),velocity_sign_transition_counts=transitions.tolist(),
            minimum_absolute_velocity=float(np.min(np.abs(states[:,1]))),
            sign_labels=['negative','positive'],rejection_histogram={str(k):int(v) for k,v in zip(uq,nq)},
            negative_jacobian_fraction=float(np.mean(determinants<0)),
            determinant_quantiles=np.quantile(determinants,[0,.01,.1,.5,.9,.99,1]).tolist(),
            consecutive_state_max_residual=float(np.max(np.abs(images[:-1]-states[1:]))))
        np.savez_compressed(args.out/f'case-{ident}.npz',rejects=rejects,determinants=determinants)
        rows.append(row)
    # Also check every pair in the dense main portrait, without subsampling.
    dense=np.load(args.render/'case-335-full.npy',mmap_mode='r')
    alternating=bool(np.all(dense[1:,1]*dense[:-1,1]<0))
    report=dict(runs=rows,dense_rosettes=dict(points=len(dense),velocity_sign_alternates_every_step=alternating,
        minimum_absolute_velocity=float(np.min(np.abs(dense[:,1]))),source_sha256=hashlib.sha256((args.render/'case-335-full.npy').read_bytes()).hexdigest()),
        scope='Finite sign/branch diagnostics, not a certified band decomposition. Reflection is an exact symmetry of the mathematical map; both visible rosettes belong to one sampled orbit. No new symmetry or attractor-class claim.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for n in ['analyze_rk45_geometry.py','search_rk45_oscillator.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--gallery',type=Path,required=True);p.add_argument('--render',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
