"""Refine the single-torus three-rate search with independent deformations and extra sites."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pyopencl as cl
from third_growth import kernel_for,run_gpu


def mutate(data,parent_index,b,batch_size):
    rng=np.random.default_rng(np.random.SeedSequence([658334,b,batch_size]))
    which=rng.choice(parent_index,size=batch_size)
    mutant=data[:,which,:].copy()
    n=batch_size
    size=np.exp(rng.uniform(np.log(.025),np.log(.65),n))
    mutant[0,:,0]=1-np.clip((1-mutant[0,:,0])*np.exp(rng.normal(0,1.2,n)),.0003,.13)
    mutant[0,:,1]*=np.exp(rng.normal(0,.29,n))
    mutant[0,:,2]=np.clip(mutant[0,:,2]*np.exp(rng.normal(0,.36,n)),.025,.97)
    added=np.zeros(n,dtype='i1')
    for j in range(n):
        count=int(mutant[0,j,3])
        sites=mutant[2:2+count,j,:3]
        spread=np.sqrt(np.var(sites,axis=0).sum())
        sites+=rng.normal(0,size[j]*spread,(count,3))
        if count<8 and rng.random()<.50:
            extras=int(rng.integers(1,min(4,8-count)+1))
            origin=sites[rng.integers(count,size=extras)]
            mutant[2+count:2+count+extras,j,:3]=origin+rng.normal(0,.4*spread,(extras,3))
            mutant[0,j,3]+=extras
            added[j]=extras
    return mutant,which,added


def search(args):
    with np.load(args.deep) as f:
        ids=f['ids']; data=f['data']; kinds=f['kind']
    rows=json.loads(args.records.read_text())['cases']
    rows=[r for r in rows if r['min_cpu_third']>.001]
    parent_index=np.array([np.flatnonzero(ids==r['index'])[0] for r in rows])
    args.out.mkdir(parents=True,exist_ok=True)
    device=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    ctx=cl.Context([device]); queue=cl.CommandQueue(ctx);kernel=kernel_for(ctx,'single')
    for b in range(args.batches):
        path=args.out/f'evolved-{b:04d}.npz'
        if path.exists():continue
        mutant,which,added=mutate(data,parent_index,b,args.size)
        n=args.size
        rates=run_gpu(queue,kernel,mutant,'single',args.burn,args.window)
        score=np.nan_to_num(rates.min(axis=(1,2)),nan=-999.)
        # Stratify by number of inserted posterior sites to retain different geometries.
        keep=[]
        for nadd in range(5):
            matches=np.flatnonzero(added==nadd)
            keep.extend(matches[np.argsort(score[matches])[-8:]])
        keep=np.unique(keep)
        temp=path.with_suffix('.tmp.npz')
        np.savez_compressed(temp,index=keep,data=mutant[:,keep,:],parent=3000000+ids[which[keep]],kind=kinds[which[keep]],added=added[keep],rates=rates[keep],settings=[n,args.burn,args.window],all_min=score)
        os.replace(temp,path)
        if b%16==0: print(b,float(score.max()),int((score>.025).sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--deep',type=Path,required=True)
    p.add_argument('--records',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--size',type=int,default=4096)
    p.add_argument('--batches',type=int,default=256)
    p.add_argument('--burn',type=int,default=1024)
    p.add_argument('--window',type=int,default=384)
    search(p.parse_args())
