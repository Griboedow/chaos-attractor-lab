"""Targeted three-vector OpenCL survey near earlier high second-growth configurations."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pyopencl as cl
from screen import parameters


def kernel_for(context, geometry):
    options=['-DSECOND_TANGENT','-DTHIRD_TANGENT']
    if geometry!='single':
        options.append('-DLINKED')
    if geometry=='smooth':
        options.append('-DSMOOTH_LINK')
    source=Path(__file__).with_name('screen.cl').read_text()
    return cl.Program(context,source).build(options=options).screen


def run_gpu(queue,kernel,data,geometry,burn,window):
    size=data.shape[1]
    m=np.empty((20 if geometry=='single' else 22,size),np.float32)
    a=cl.Buffer(queue.context,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=np.ascontiguousarray(data))
    b=cl.Buffer(queue.context,cl.mem_flags.WRITE_ONLY,m.nbytes)
    kernel(queue,(size,),None,a,b,np.int32(size),np.int32(burn),np.int32(window))
    cl.enqueue_copy(queue,m,b).wait()
    offset=12 if geometry=='single' else 14
    return np.stack([m[:4].T,m[offset:offset+4].T,m[offset+4:offset+8].T],axis=2)


def anchors(root,geometry):
    prefix={'single':'','linked':'linked-','smooth':'smooth-'}[geometry]
    offset={'single':0,'linked':1000000,'smooth':2000000}[geometry]
    audit=json.loads((root/(prefix+'two-rate-audit.json')).read_text())
    rows=sorted(audit['cases'],key=lambda v:np.mean(v['cpu_qr'],axis=0)[2],reverse=True)[:14]
    entries=[]
    for row in rows:
        i=row['case']-offset
        raw,kind=parameters(i//4096,4096,geometry)
        entries.append((row['case'],raw[:,i%4096,:].copy(),int(kind[i%4096])))
    return entries


def draw(entries,geometry,batch,n):
    rng=np.random.default_rng(np.random.SeedSequence([810271,{'single':1,'linked':2,'smooth':3}[geometry],batch,n]))
    ids=rng.integers(len(entries),size=n)
    data=np.stack([entries[i][1] for i in ids],axis=1)
    parent=np.array([entries[i][0] for i in ids])
    kind=np.array([entries[i][2] for i in ids])
    # Half of the proposals follow near-neutral parents closely. The remainder
    # deform all sites and vary their common center over a wider range.
    amplitude=np.exp(rng.uniform(np.log(.009),np.log(1.0),size=n))
    data[0,:,0]=1-np.exp(rng.uniform(np.log(.0003),np.log(.13),size=n))
    data[0,:,1]*=np.exp(rng.normal(0,.65,size=n))
    cap=.98 if geometry=='single' else .49
    data[0,:,2]=np.clip(data[0,:,2]*np.exp(rng.normal(0,.42,size=n)),.025,cap)
    if geometry=='smooth':
        data[2,:,3]*=np.exp(rng.normal(0,.75,size=n))
    shift=rng.normal(size=(n,3))*amplitude[:,None]
    shear=np.eye(3)[None]+rng.normal(size=(n,3,3))*amplitude[:,None,None]*.48
    for j in range(n):
        count=int(data[0,j,3])
        sites=data[2:2+count,j,:3]
        center=sites.mean(0)
        sites[:]=(sites-center)@shear[j]+center+shift[j]+rng.normal(size=(count,3))*amplitude[j]*.40
    return data,parent,kind


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    device=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    context=cl.Context([device]);queue=cl.CommandQueue(context)
    for geometry in args.geometries:
        kernel=kernel_for(context,geometry)
        entries=anchors(args.source,geometry)
        for batch in range(args.batches):
            path=args.out/f'{geometry}-{batch:04d}.npz'
            if path.exists():continue
            data,parent,kind=draw(entries,geometry,batch,args.size)
            rates=run_gpu(queue,kernel,data,geometry,args.burn,args.window)
            score=np.nan_to_num(rates.min(axis=(1,2)),nan=-999.)
            keep=np.unique(np.r_[np.flatnonzero(score>.006),np.argsort(score)[-12:]])
            temp=path.with_suffix('.tmp.npz')
            np.savez_compressed(temp,index=keep,data=data[:,keep,:],parent=parent[keep],kind=kind[keep],rates=rates[keep],settings=[args.size,args.burn,args.window],threshold=.006,all_third_min=rates[:,:,2].min(1),all_min=score)
            os.replace(temp,path)
            if batch%8==0:print(geometry,batch,float(score.max()),int((score>.006).sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--geometries',nargs='+',choices=('single','linked','smooth'),default=['single','smooth','linked'])
    p.add_argument('--size',type=int,default=4096)
    p.add_argument('--batches',type=int,default=64)
    p.add_argument('--burn',type=int,default=1024)
    p.add_argument('--window',type=int,default=384)
    main(p.parse_args())
