"""Algebraic posterior identities and a declared simple-parameter sphere survey.

Rotation equivalence is exact for ideal regular polyhedra. Rounded stored sites
are separately compared numerically; they are not silently replaced in old runs.
"""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from numba import set_num_threads
from search_soft_sphere import screen, step, DEPENDENCIES
from search_soft_constellations import posterior
from search_raar_constellations import vertices


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p, x): p.write_text(json.dumps(x, indent=2)+'\n', encoding='utf-8', newline='\n')


def spin_mean(h, kind):
    t = np.tanh(h)
    if kind == 2:
        return t
    return (t + np.array([t[1]*t[2], t[0]*t[2], t[0]*t[1]]))/(1+np.prod(t))


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    set_num_threads(4)
    begin = time.time()
    for name in [Path(__file__).name]+DEPENDENCIES:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rng = np.random.default_rng(2026110130)
    identities = []
    for kind in [0, 2]:
        signs = vertices(kind)*np.sqrt(3)
        errors = []
        for _ in range(128):
            # Moderate fields avoid cancellation in the tetrahedral tanh quotient.
            h = rng.uniform(-2, 2, 3)
            logits = signs@h
            weights = np.exp(logits-logits.max()); weights /= weights.sum()
            errors.append(float(np.max(abs(weights@signs-spin_mean(h, kind)))))
        assert max(errors) < 1e-11
        identities.append(dict(kind=kind, samples=128, field_range=[-2,2], max_error=max(errors)))
    rows, canon, params, points, counts, kinds = [], [], [], [], [], []

    def add(kind, beta, tau, a, c, provenance):
        signs = np.rint(vertices(kind)*np.sqrt(3))
        sites = np.asarray(c)+a*signs
        index = len(params)
        params.append([beta,tau]);points.append(np.pad(sites,((0,12-len(sites)),(0,0))))
        counts.append(len(sites));kinds.append(kind)
        rows.append(dict(id=index,kind=kind,beta=float(beta),tau=float(tau),a=float(a),center=np.asarray(c).tolist(),selection=provenance))

    with np.load(args.hard_source) as hard, np.load(args.soft_source) as soft:
        for base, old in [(1086,155),(1504,251)]:
            p=hard['params'][base];kind=int(hard['kinds'][base]);s,t,u=p[5:]
            quat=[np.sqrt(1-s)*np.sin(2*np.pi*t),np.sqrt(1-s)*np.cos(2*np.pi*t),np.sqrt(s)*np.sin(2*np.pi*u),np.sqrt(s)*np.cos(2*np.pi*u)]
            Q=Rotation.from_quat(quat).as_matrix();c=Q.T@p[2:5];a=p[1]/np.sqrt(3)
            beta,tau=soft['params'][old];original=soft['points'][old,:soft['counts'][old]]
            ideal=p[2:5]+a*np.rint(vertices(kind)*np.sqrt(3))@Q.T
            canonical=c+a*np.rint(vertices(kind)*np.sqrt(3))
            mean_errors=[];step_errors=[]
            for _ in range(128):
                z=rng.uniform(-2,2,3);x=Q@z
                mean_errors.append(float(np.max(abs(posterior(z,canonical,tau)[0]-(c+a*spin_mean(2*a*(z-c)/tau,kind))))) if kind==2 else None)
                step_errors.append(float(np.max(abs(Q.T@step(x,beta,original,tau)[0]-step(z,beta,canonical,tau)[0]))))
            canon.append(dict(old_case=old,base=base,kind=kind,a=float(a),center=c.tolist(),rotation=Q.tolist(),site_reconstruction_error=float(np.max(abs(original-ideal))),orthogonality_error=float(np.max(abs(Q.T@Q-np.eye(3)))),sampled_rotation_step_error=max(step_errors),cube_mean_error=max(mean_errors) if kind==2 else None,scope='Floating-point residuals, not exact identities for the original rounded site arrays.'))
            add(kind,beta,tau,a,c,'Ideal canonical reconstruction of old Case '+str(old)+'; a changed binary64 input.')
            for decimals in [3,2,1]:
                add(kind,round(float(beta),max(2,decimals)),round(float(tau),max(3,decimals)),round(float(a),decimals),np.round(c,decimals),'Joint decimal simplification of old Case '+str(old)+f'; {decimals} geometry decimal places.')
    for beta,tau,a,c in itertools.product([.9,.93,.96,.99],[.01,.02,.04,.06,.1,.2],[.7,.8,.9],[[-.1,-.2,.2],[0,0,0],[0,0,.2],[.1,.1,.1],[-.1,-.1,.2]]):
        add(0,beta,tau,a,c,'Complete declared tetrahedron grid.')
    for beta,tau,a,c in itertools.product([.9,.93,.96,.99],[.002,.005,.01,.02,.04,.08],[.25,.3,.35],[[-.6,-.5,-.5],[-.5,-.5,-.5],[-.6,-.4,-.2],[0,0,0]]):
        add(2,beta,tau,a,c,'Complete declared cube grid.')
    params=np.array(params);points=np.array(points);counts=np.array(counts);kinds=np.array(kinds)
    protocol=dict(hard_source=str(args.hard_source),hard_source_sha256=sha(args.hard_source),soft_source=str(args.soft_source),soft_source_sha256=sha(args.soft_source),total=len(params),identities=identities,canonical_reconstructions=canon,cases=rows,warm=5000,observed=15000,windows=10,initial=[.317,-.211,.419],flag='Complete, every leading window greater than 0.005.',scope='Known posterior and reflection ingredients. Simple-parameter exploration, not a priority or asymptotic certificate.',source_sha256=sha(Path(__file__)))
    save(args.out/'protocol.json',protocol)
    parts=[[] for _ in range(6)]
    for offset in range(0,len(params),64):
        s=slice(offset,offset+64);values=screen(params[s],points[s],counts[s],5000,15000)
        for dest,value in zip(parts,values):dest.append(value)
        print(json.dumps(dict(completed=min(offset+64,len(params)),total=len(params),seconds=time.time()-begin)),flush=True)
    rates,tails,status,minima,stats,ratios=[np.concatenate(v) for v in parts]
    flags=(status[:,0]==0)&(rates.min(1)>.005)
    np.savez_compressed(args.out/'screen.npz',params=params,points=points,counts=counts,kinds=kinds,base_ids=np.full(len(params),-1),relative_tau=np.full(len(params),np.nan),rates=rates,tails=tails,status=status,minima=minima,soft_stats=stats,boundratios=ratios,flags=flags)
    save(args.out/'report.json',dict(protocol,complete=int((status[:,0]==0).sum()),flags=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),seconds=time.time()-begin))
    print(json.dumps(dict(complete=int((status[:,0]==0).sum()),flags=int(flags.sum()),seconds=time.time()-begin)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--hard-source',type=Path,required=True);p.add_argument('--soft-source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
