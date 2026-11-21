"""Test both leading tangent directions on all surviving spatial GPU flags."""
import argparse
import os
from pathlib import Path
import numpy as np
import pyopencl as cl


def main(args):
    dev=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type & cl.device_type.GPU)
    ctx=cl.Context([dev]);q=cl.CommandQueue(ctx)
    args.out.mkdir(parents=True,exist_ok=True)
    linked=args.geometry!='single'
    flags=['-DSECOND_TANGENT','-DSAVE_CLOUD']+(['-DLINKED'] if linked else [])+(['-DSMOOTH_LINK'] if args.geometry=='smooth' else [])
    src=Path(__file__).with_name('screen.cl').read_text()
    kernel=cl.Program(ctx,src).build(options=flags).screen
    tally=np.zeros(2,dtype=int)
    for path in sorted(args.deep.glob('deep-*.npz')):
        dest=args.out/path.name
        if dest.exists(): continue
        with np.load(path) as f:
            ids=f['ids']; kinds=f['kinds'];raw=np.ascontiguousarray(f['data'])
        n=len(ids)
        if n==0: continue
        m=np.empty((18 if linked else 16,n),np.float32)
        z=np.empty((args.tail,n,4),np.float32)
        a=cl.Buffer(ctx,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=raw)
        b=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,m.nbytes)
        c=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,z.nbytes)
        kernel(q,(n,),None,a,b,np.int32(n),np.int32(args.burn),np.int32(args.window),c,np.int32(args.tail))
        cl.enqueue_copy(q,m,b).wait();cl.enqueue_copy(q,z,c).wait()
        keep=(np.min(m[-4:],axis=0)>.025)&(np.min(m[:4],axis=0)>.025)&np.isfinite(z).all(axis=(0,2))
        if linked: keep&=(m[12]>.025)&(m[12]<.975)
        tmp=dest.with_suffix('.tmp.npz')
        np.savez(tmp,ids=ids[keep],kinds=kinds[keep],metrics=m[:,keep],all_metrics=m,all_ids=ids,points=z[:,keep,:3].transpose(1,0,2),data=raw[:,keep,:],settings=[args.burn,args.window,args.tail])
        os.replace(tmp,dest)
        tally+=[n,keep.sum()]
    print(args.geometry,tally.tolist(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--deep',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--geometry',choices=('single','linked','smooth'),required=True)
    p.add_argument('--burn',type=int,default=8192)
    p.add_argument('--window',type=int,default=2048)
    p.add_argument('--tail',type=int,default=1536)
    main(p.parse_args())
