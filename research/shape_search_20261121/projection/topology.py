"""Finite-cloud alpha-complex persistence in equal spatial units."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
import gudhi


def barcode(z):
    tree=gudhi.AlphaComplex(points=z,precision='safe').create_simplex_tree()
    tree.compute_persistence(homology_coeff_field=2)
    groups={}
    for k in (0,1,2):
        intervals=tree.persistence_intervals_in_dimension(k)
        intervals=intervals[np.isfinite(intervals).all(axis=1)]
        # Default filtration values are squared ball radii in the data units.
        r=np.sqrt(np.maximum(intervals,0))
        lengths=r[:,1]-r[:,0]
        idx=np.argsort(lengths)[::-1][:6]
        groups[str(k)]=dict(bars=[[float(r[i,0]),float(r[i,1])] for i in idx],largest=float(lengths.max(initial=0)),counts=[int(np.sum(lengths>t)) for t in (.10,.20,.35)])
    return groups


def analyze(z,case,seed,sizes):
    z=np.asarray(z,dtype=float)
    centered=z-z.mean(axis=0)
    scale=np.sqrt(np.sum(centered**2,axis=1).mean())
    centered/=scale
    rng=np.random.default_rng(seed)
    half=len(z)//2
    reps=[]
    for side,data in enumerate((centered[:half],centered[half:])):
        for size in sizes:
            points=data[rng.choice(len(data),size,replace=False)]
            reps.append(dict(half=side+1,count=size,groups=barcode(points)))
    return dict(case=case,seed=seed,source_count=len(z),rms_radius=scale,coefficient_field=2,gudhi=gudhi.__version__,sizes=sizes,replicates=reps)


def save(path,doc):
    tmp=path.with_suffix('.tmp.json')
    tmp.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
    os.replace(tmp,path)


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    controls=args.out/'controls.json'
    if not controls.exists():
        rng=np.random.default_rng(7833)
        n=20000; theta=rng.uniform(-np.pi,np.pi,n);phi=rng.uniform(-np.pi,np.pi,n)
        sphere=rng.normal(size=(n,3));sphere/=np.linalg.norm(sphere,axis=1,keepdims=True)
        models={'circle':np.column_stack([np.cos(theta),np.sin(theta),np.zeros(n)]),
                'sphere':sphere,
                'ring_torus':np.column_stack([(1+.38*np.cos(phi))*np.cos(theta),(1+.38*np.cos(phi))*np.sin(theta),.38*np.sin(phi)]),
                'gaussian':rng.normal(size=(n,3))}
        # Analytic circle and sphere clouds are exactly co-spherical. Add tiny
        # deterministic jitter to resolve triangulation degeneracies in controls.
        control_rows=[]
        for j,(name,cloud) in enumerate(models.items()):
            cloud=cloud+1e-9*rng.normal(size=cloud.shape)
            control_rows.append(analyze(cloud,name,7834+j,args.sizes))
        save(controls,dict(generator_seed=7833,jitter_sigma=1e-9,models=control_rows))
    for path in sorted(args.source.glob('case-*.json')):
        r=json.loads(path.read_text())
        if not r.get('three_positive_passes',r['passes']): continue
        dest=args.out/f'persistence-{r["case"]}.json'
        if dest.exists(): continue
        with np.load(r['cloud']) as f: cloud=f['points']
        row=analyze(cloud,r['case'],28132+r['case'],args.sizes)
        save(dest,row)
        print(r['case'],[round(min(e['groups'][str(k)]['largest'] for e in row['replicates']),3) for k in (1,2)],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--sizes',nargs='+',type=int,default=[1536,3072])
    main(p.parse_args())
