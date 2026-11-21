"""Rescreen fast flags and save consecutive, unjoined native iterates in the cache."""
import argparse
import os
from pathlib import Path
import numpy as np
import pyopencl as cl
from screen import parameters


def eigenvalues(metrics):
    c = np.empty((metrics.shape[1],3,3),float)
    c[:,0,0],c[:,1,1],c[:,2,2] = metrics[4:7]
    c[:,0,1]=c[:,1,0]=metrics[7]
    c[:,0,2]=c[:,2,0]=metrics[8]
    c[:,1,2]=c[:,2,1]=metrics[9]
    return np.linalg.eigvalsh(np.nan_to_num(c))


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    dev = next(d for p in cl.get_platforms() for d in p.get_devices() if d.type & cl.device_type.GPU)
    ctx=cl.Context([dev])
    q=cl.CommandQueue(ctx)
    source=Path(__file__).with_name('screen.cl').read_text()
    kernels={}
    totals=np.zeros(3,dtype=int)
    for path in sorted(args.survey.glob('batch-*.npz')):
        batch=int(path.stem.split('-')[1])
        target=args.out/f'deep-{batch:05d}.npz'
        if target.exists():
            continue
        with np.load(path) as f:
            metrics=f['metrics']
            sizes=f['settings']
            geometry=str(f['geometry']) if 'geometry' in f else 'single'
            eigen=eigenvalues(metrics)
            mask=(metrics[:4].min(axis=0)>.015)&(eigen[:,0]>.003*eigen.sum(axis=1))&(eigen[:,1]>.04*eigen.sum(axis=1))&(eigen.sum(axis=1)>.03)
            if geometry!='single':
                mask &= (metrics[12]>.025)&(metrics[12]<.975)
        indices=np.flatnonzero(mask)
        linked=geometry!='single'
        if geometry not in kernels:
            flags=['-DSAVE_CLOUD']+(['-DLINKED'] if linked else [])+(['-DSMOOTH_LINK'] if geometry=='smooth' else [])
            kernels[geometry]=cl.Program(ctx,source).build(options=flags).screen
        kernel=kernels[geometry]
        data,kind=parameters(batch,int(sizes[0]),geometry)
        count=len(indices)
        if count==0:
            raise RuntimeError(f'No spatial flags in batch {batch}')
        data=np.ascontiguousarray(data[:,indices,:])
        metrics=np.empty((14 if linked else 12,count),np.float32)
        points=np.empty((args.tail,count,4),np.float32)
        ib=cl.Buffer(ctx,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=data)
        ob=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,metrics.nbytes)
        pb=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,points.nbytes)
        kernel(q,(count,),None,ib,ob,np.int32(count),np.int32(args.burn),np.int32(args.window),pb,np.int32(args.tail))
        cl.enqueue_copy(q,metrics,ob).wait()
        cl.enqueue_copy(q,points,pb).wait()
        keep=(metrics[:4].min(axis=0)>.02)&np.isfinite(points).all(axis=(0,2))
        # Require both rings to contribute in the deeper linked-surface run.
        if linked:
            keep &= (metrics[12]>.025)&(metrics[12]<.975)
        offset={'single':0,'linked':1000000,'smooth':2000000}[geometry]
        tmp=target.with_suffix('.tmp.npz')
        np.savez(tmp,ids=offset+batch*int(sizes[0])+indices[keep],data=data[:,keep,:],kinds=kind[indices][keep],points=points[:,keep,:3].transpose(1,0,2),metrics=metrics[:,keep],all_ids=offset+batch*int(sizes[0])+indices,all_metrics=metrics,settings=[args.burn,args.window,args.tail])
        os.replace(tmp,target)
        totals+=[int(sizes[0]),count,keep.sum()]
        if batch%8==0: print(batch,totals.tolist(),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--survey',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--burn',type=int,default=8192)
    p.add_argument('--window',type=int,default=1024)
    p.add_argument('--tail',type=int,default=768)
    main(p.parse_args())
