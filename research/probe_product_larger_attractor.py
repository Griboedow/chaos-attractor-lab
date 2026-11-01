"""Trace both numerical unstable branches and test nearby finite-time attraction.

The exact local graph is known from a separate certificate. Here its unknown
stable coordinates are set to zero for proposals, and explicitly not certified.
High precision is used before the exponentially small launch becomes visible.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from numba import njit, prange, set_num_threads
from scipy.spatial import cKDTree
from search_product_alphabet import step


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(p, value):
    p.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


@njit(cache=True, parallel=True)
def follow(starts, beta, tau, levels, warm, count, keep):
    tails=np.empty((len(starts),keep,3)); rates=np.zeros((len(starts),10))
    minima=np.full(len(starts),np.inf); flags=np.zeros(len(starts),np.int64)
    for k in prange(len(starts)):
        x=starts[k].copy(); w=np.array([1.,.371,-.217]); w/=np.linalg.norm(w)
        for it in range(warm+count):
            x,J,r,_,_,flag=step(x,beta,tau,levels,2); minima[k]=min(minima[k],r)
            if flag or not np.isfinite(x).all(): flags[k]=flag or 3; break
            w=J@w; length=np.linalg.norm(w)
            if length<1e-280 or not np.isfinite(length): flags[k]=4; break
            w/=length
            if it>=warm: rates[k,(it-warm)//(count//10)]+=np.log(length)
            if it>=warm+count-keep: tails[k,it-(warm+count-keep)]=x
    return tails,rates/(count//10),minima,flags


@njit(cache=True, parallel=True)
def fan(starts,beta,tau,levels,steps):
    states=np.empty((len(starts),steps+1,3)); minima=np.full(len(starts),np.inf)
    for k in prange(len(starts)):
        x=starts[k].copy(); states[k,0]=x
        for it in range(steps):
            x,_,r,_,_,flag=step(x,beta,tau,levels,2)
            assert not flag and np.isfinite(x).all()
            states[k,it+1]=x;minima[k]=min(minima[k],r)
    return states,minima


def stats(values):
    return dict(zip(['median','q90','q99','maximum'],map(float,np.quantile(values,[.5,.9,.99,1]))))


def main(args):
    begin=time.time();args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    for name in [Path(__file__).name,'search_product_alphabet.py','product_binary_interval.py',
                 'follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py',
                 'search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    raw=json.loads(args.input.read_text());cert=json.loads(args.certificate.read_text())
    assert cert['status']=='passed' and cert['windows_sha256']==sha(args.input)
    with np.load(args.source) as f: reference=f['states'];beta,tau=f['params'];levels=f['levels']
    sys.path.append(str(args.mpmath_path));import mpmath as mp
    from product_binary_interval import parameters,model
    mp.mp.dps=100;par=parameters(raw['params'],np.array(raw['points']),mp.mp)
    exact=lambda s:mp.mpf(Fraction(s).numerator)/Fraction(s).denominator
    p=list(map(exact,cert['root']['center_decimal']));E=[[exact(v) for v in row] for row in cert['graph']['E']]
    # Positive magnitudes span one approximate fundamental interval for G^2:
    # lambda_u^2 ~= 131.36, preserving each branch's sign. This is a numerical
    # sampling convention; coverage of a true fundamental domain is not proved.
    magnitudes=np.geomspace(1e-27,1.4e-25,args.samples_per_branch)
    signed=np.concatenate([magnitudes,-magnitudes]);starts=[];prefix=[]
    for k,u in enumerate(signed):
        u=exact(float(u));x=[p[i]+E[i][0]*u for i in range(3)];trace=[]
        for it in range(89):
            trace.append([float(v) for v in x])
            if it<88:x,_,_=model(x,par,mp.mp)
        starts.append(trace[-1]);prefix.append(trace)
        if k%32==31:print(json.dumps(dict(stage='high_precision_launches',done=k+1,total=len(signed))),flush=True)
    starts=np.array(starts);prefix=np.array(prefix)
    states,minima=fan(starts,beta,tau,levels,args.fan_steps)
    np.savez_compressed(args.out/'unstable-fan.npz',signed_u=signed,starts=starts,prefix=prefix,states=states,minimum_radii=minima)
    print(json.dumps(dict(stage='fan_complete',seconds=time.time()-begin)),flush=True)
    # A declared local perturbation cohort, with deterministic spherical directions.
    n=args.directions;i=np.arange(n);z=1-2*(i+.5)/n;phi=i*np.pi*(3-np.sqrt(5))
    directions=np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
    radii=np.array([1e-6,1e-4,1e-2,.1]);root=np.array(list(map(float,p)))
    nearby=np.concatenate([root+r*directions for r in radii])
    tails,rates,mins,flags=follow(nearby,beta,tau,levels,args.warm,args.count,4096)
    assert not flags.any()
    np.savez_compressed(args.out/'nearby-cohorts.npz',starts=nearby,radii=radii,directions=directions,tails=tails,rates=rates,minimum_radii=mins,flags=flags)
    # Two long continuations, one per outgoing branch. They are proposals, not
    # enclosures of the exact unknown graph points after the finite certificate.
    chosen=np.array([args.samples_per_branch-1,2*args.samples_per_branch-1])
    long,rates_long,mins_long,flags_long=follow(starts[chosen],beta,tau,levels,20000,1000000,1000000)
    assert not flags_long.any()
    np.savez_compressed(args.out/'branch-long.npz',starts=starts[chosen],states=long,rates=rates_long,minimum_radii=mins_long,flags=flags_long)
    tree=cKDTree(reference);allfan=states.reshape(-1,3);fd=tree.query(allfan,workers=4)[0].reshape(states.shape[:2])
    nd=tree.query(tails.reshape(-1,3),workers=4)[0].reshape(tails.shape[:2])
    ld=tree.query(long.reshape(-1,3),workers=4)[0].reshape(long.shape[:2])
    # Matched sample counts for support comparison: reverse nearest-neighbor
    # distances from the first 200k reference points to each 1m branch tail,
    # and to a disjoint 1m reference tail. No alignment or rescaling is used.
    query=reference[:200000];control=reference[-1000000:]
    reverse=[stats(cKDTree(v).query(query,workers=4)[0]) for v in long]
    reverse_control=stats(cKDTree(control).query(query,workers=4)[0])
    # Coarse occupancy and frequency comparisons on a fixed physical grid.
    # Exact equality of invariant measures is not tested or asserted.
    origin=np.array([-.5,-.5,-.5]);width=.02
    def frequencies(x):
        v=np.floor((x-origin)/width).astype(np.int64)
        keys,counts=np.unique(v,axis=0,return_counts=True)
        return {tuple(k):float(c/len(x)) for k,c in zip(keys,counts)}
    def compare(a,b):
        keys=set(a)|set(b)
        return dict(total_variation=.5*sum(abs(a.get(k,0)-b.get(k,0)) for k in keys),
                    reference_mass_in_shared_cells=sum(v for k,v in a.items() if k in b),
                    branch_mass_in_shared_cells=sum(v for k,v in b.items() if k in a),
                    reference_cells=len(a),branch_cells=len(b))
    base=frequencies(reference[:1000000]);freq=[compare(base,frequencies(v)) for v in long]
    freq_control=compare(base,frequencies(control))
    np.savez_compressed(args.out/'distances.npz',fan=fd,nearby=nd,long=ld)
    cohorts=[]
    for k,r in enumerate(radii):
        sl=slice(k*n,(k+1)*n)
        cohorts.append(dict(radius=float(r),count=n,positive_every_block=int((rates[sl].min(1)>0).sum()),
                            leading_mean_range=[float(rates[sl].mean(1).min()),float(rates[sl].mean(1).max())],
                            tail_distance=stats(nd[sl].ravel()),minimum_normalization=float(mins[sl].min())))
    report=dict(status='completed',sources={str(p):sha(p) for p in [args.input,args.certificate,args.source]},
        parameters=raw['params'],samples_per_branch=args.samples_per_branch,launch_magnitudes=[float(magnitudes[0]),float(magnitudes[-1])],
        launch_digits=100,high_precision_steps=88,fan_steps=args.fan_steps,fan_distance=stats(fd.ravel()),
        fan_checkpoints={str(t+88):stats(fd[:,t]) for t in [0,4,16,64,256,1024,args.fan_steps]},
        perturbation_protocol=dict(directions=n,radii=radii.tolist(),warm=args.warm,count=args.count,tail=4096,blocks=10),
        nearby=cohorts,long_branch=dict(indices=chosen.tolist(),warm=20000,count=1000000,
          mean_leading=rates_long.mean(1).tolist(),minimum_leading_block=rates_long.min(1).tolist(),
          forward_distance=[stats(d) for d in ld],reverse_distance=reverse,reverse_control=reverse_control,
          grid_origin=origin.tolist(),grid_width=width,grid_comparison=freq,grid_control=freq_control),
        seconds=time.time()-begin,scope='Finite numerical evidence for an unstable-core-in-larger-regime hypothesis. Eigenline launches approximate the certified local graph but do not certify its global manifold, attraction, recurrence, transitivity, or historical novelty.')
    save(args.out/'report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','certificate','source','mpmath-path','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--samples-per-branch',type=int,default=128);p.add_argument('--fan-steps',type=int,default=2048)
    p.add_argument('--directions',type=int,default=64);p.add_argument('--warm',type=int,default=20000);p.add_argument('--count',type=int,default=40000)
    main(p.parse_args())
