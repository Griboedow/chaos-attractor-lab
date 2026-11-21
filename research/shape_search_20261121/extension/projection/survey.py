"""GPU search of posterior reflections onto ellipsoidal images of a torus."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pyopencl as cl

BASE=Path(__file__).resolve().parents[2]/'projection'
sys.path.insert(0,str(BASE))
from reproduce import raw_parameters
from third_growth import run_gpu


def transformed_kernel(context,save=False):
    src=(BASE/'screen.cl').read_text()
    before='#ifdef SMOOTH_LINK\n        torus_frame a='
    after='        x=beta*(x+p)+(1-2*beta)*m;'
    assert src.count(before)==src.count(after)==1
    src=src.replace(before,'''        float3 scale=(float3)(data[2*stride+i].w,data[3*stride+i].w,data[4*stride+i].w);
        v/=scale; dv/=scale; dv2/=scale; dv3/=scale;
        float warp=data[5*stride+i].w;
        if (warp!=0) {
            dv.z-=2*warp*(v.x*dv.x-v.y*dv.y); dv2.z-=2*warp*(v.x*dv2.x-v.y*dv2.y); dv3.z-=2*warp*(v.x*dv3.x-v.y*dv3.y);
            v.z-=warp*(v.x*v.x-v.y*v.y);
        }
'''+before)
    src=src.replace(after,'''        if (warp!=0) {
            dp.z+=2*warp*(p.x*dp.x-p.y*dp.y); dp2.z+=2*warp*(p.x*dp2.x-p.y*dp2.y); dp3.z+=2*warp*(p.x*dp3.x-p.y*dp3.y);
            p.z+=warp*(p.x*p.x-p.y*p.y);
        }
        p*=scale; dp*=scale; dp2*=scale; dp3*=scale;
'''+after)
    return cl.Program(context,src).build(options=['-DSECOND_TANGENT','-DTHIRD_TANGENT']+(['-DSAVE_CLOUD'] if save else [])).screen


def get_parents(study):
    parents=[]
    for row in json.loads((study/'three'/'manifest.json').read_text()):
        rec=json.loads((study/'three'/row['record']).read_text())
        parents.append((row['case'],raw_parameters(rec).astype('f4')))
    return parents


def draw(batch,n,parents,warp=False):
    rng=np.random.default_rng(np.random.SeedSequence([489824,20261121,batch,n]))
    selected=rng.integers(0,len(parents),n)
    data=np.asarray([p[1] for p in parents])[selected].transpose(1,0,2).copy()
    amp=np.exp(rng.uniform(np.log(.02),np.log(.83),size=n))
    scale=np.exp(rng.normal(0,1,(n,3))*rng.uniform(.09,.65,(n,1)))
    scale=np.clip(scale,.24,3.5)
    if warp:
        scale=np.clip(scale*np.sqrt(data[2:5,:,3].T),.24,3.5)
        data[5,:,3]=rng.choice([-1.,1.],n)*np.exp(rng.uniform(np.log(.035),np.log(2.3),n))
    data[2:5,:,3]=scale.T
    data[0,:,0]=1-np.clip((1-data[0,:,0])*np.exp(rng.normal(0,1,n)),.0003,.25)
    data[0,:,1]*=np.exp(rng.normal(0,.34,n))
    data[0,:,2]=np.clip(data[0,:,2]*np.exp(rng.normal(0,.43,n)),.024,.95)
    mask=np.arange(8)[:,None]<data[0,:,3][None,:]
    coords=data[2:10,:,:3]
    # Some families co-deform their sites with the surface; all keep Euclidean
    # likelihoods in physical coordinates, breaking affine conjugacy in general.
    covary=rng.random(n)<.55
    coords*=np.where(covary[:,None],scale,1)[None,:,:]
    coords+=((rng.normal(size=(8,n,3))*.65+rng.normal(size=(1,n,3))*.65)*amp[None,:,None]*mask[:,:,None]).astype('f4')
    return data,np.array([parents[i][0] for i in selected])


def main(args):
    args.cache.mkdir(parents=True,exist_ok=True)
    parents=get_parents(args.study)
    if args.parent_shortlist:
        with np.load(args.parent_shortlist) as f:parents=[(int(i),f['raw'][:,k,:]) for k,i in enumerate(f['ids'])]
    dev=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    ctx=cl.Context([dev]);q=cl.CommandQueue(ctx);kernel=transformed_kernel(ctx)
    for b in range(args.first,args.first+args.batches):
        dest=args.cache/f'screen-{b:04d}.npz'
        if dest.exists():continue
        data,parents_of=draw(b,args.size,parents,args.warp)
        rates=run_gpu(q,kernel,data,'single',1024,384)
        score=np.nan_to_num(rates.min(axis=(1,2)),nan=-999.)
        axis=data[2:5,:,3];anis=np.max(axis,0)/np.min(axis,0)
        keep=[]
        for lo,hi in [(1,1.3),(1.3,1.7),(1.7,2.4),(2.4,4),(4,100)]:
            ix=np.flatnonzero((anis>=lo)&(anis<hi))
            keep.extend(ix[np.argsort(score[ix])[-9:]])
        keep=np.unique(keep)
        raw=np.ascontiguousarray(data.transpose(1,0,2))
        digests=np.frombuffer(b''.join(hashlib.blake2b(row,digest_size=16).digest() for row in raw),dtype='V16')
        np.savez_compressed(dest,ids=(2000000 if args.warp else 0)+b*args.size+keep,data=data[:,keep,:],parents=parents_of[keep],rates=rates[keep],scores=score,hashes=digests)
        if b%16==0:print(b,float(score.max()),int((score>.01).sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--first',type=int,default=0);p.add_argument('--batches',type=int,default=256);p.add_argument('--size',type=int,default=4096)
    p.add_argument('--warp',action='store_true');p.add_argument('--parent-shortlist',type=Path)
    main(p.parse_args())
