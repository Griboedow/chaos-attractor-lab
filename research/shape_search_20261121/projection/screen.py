"""Reproducible OpenCL survey of a posterior mean and a ring torus."""
import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pyopencl as cl
from scipy.spatial.transform import Rotation

SEED = 202611210
BASES = (
    [[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],
    [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],
    [[x,y,z] for x in (-1,1) for y in (-1,1) for z in (-1,1)],
)
KINDS = ('tetrahedron','octahedron','cuboid','random six','perturbed eight')


def parameters(batch, size, geometry='single'):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED,batch,size])))
    data = np.zeros((10,size,4),dtype=np.float32)
    kind = rng.integers(0,len(KINDS),size=size)
    data[0,:,0] = 1 - np.exp(rng.uniform(np.log(.002),np.log(.35),size))
    data[0,:,1] = 1 / np.exp(rng.uniform(np.log(.018),np.log(1.6),size))
    data[0,:,2] = rng.uniform(.07,.92,size)
    data[1,:,:3] = [.317,-.211,.419]
    axes = np.exp(rng.uniform(np.log(.30),np.log(2.10),(size,1,3)))
    center = rng.uniform(-.95,.95,(size,1,3))
    if geometry in ('linked','smooth'):
        data[0,:,2] = .08 + (.46-.08)*(data[0,:,2]-.07)/(.92-.07)
        data[1,:,3] = 1 if geometry=='linked' else 2
        center[:,:,0] += .5
    rotations = Rotation.random(size,random_state=rng).as_matrix()
    for k in range(len(KINDS)):
        ids = np.flatnonzero(kind==k)
        if k < 3:
            points = np.asarray(BASES[k],float)
            points /= np.sqrt(np.mean(np.sum(points**2,axis=1)))
            points = np.broadcast_to(points,(len(ids),*points.shape))
        elif k == 3:
            points = rng.normal(size=(len(ids),6,3))
            points -= points.mean(axis=1,keepdims=True)
            points /= np.sqrt((points**2).sum(axis=2).mean(axis=1))[:,None,None]
        else:
            points = np.asarray(BASES[2],float)[None,:,:]/np.sqrt(3)
            points = points + rng.normal(0,.20,(len(ids),8,3))
        points = (points * axes[ids]) @ rotations[ids] + center[ids]
        n = points.shape[1]
        data[0,ids,3] = n
        data[2:2+n,ids,:3] = points.transpose(1,0,2)
    if geometry=='smooth':
        data[2,:,3]=np.exp(rng.uniform(np.log(.015),np.log(2.5),size))
    return data, kind


def screen(args):
    devices = [d for p in cl.get_platforms() for d in p.get_devices(device_type=cl.device_type.GPU)]
    ctx = cl.Context([devices[args.gpu]])
    q = cl.CommandQueue(ctx)
    linked=args.geometry!='single'
    options=(['-DLINKED'] if linked else [])+(['-DSMOOTH_LINK'] if args.geometry=='smooth' else [])
    kernel = cl.Program(ctx,Path(__file__).with_name('screen.cl').read_text(encoding='utf8')).build(options=options).screen
    args.cache.mkdir(parents=True,exist_ok=True)
    print(f'GPU {q.device.name}, {args.size} cases per batch, {args.burn}+4*{args.window}',flush=True)
    print('batch cases flags seconds',flush=True)
    times = []
    for batch in range(args.first,args.first+args.batches):
        target = args.cache / f'batch-{batch:05d}.npz'
        if target.exists():
            continue
        start = time.perf_counter()
        data,kind = parameters(batch,args.size,args.geometry)
        metrics = np.empty((14 if linked else 12,args.size),np.float32)
        input_b = cl.Buffer(ctx,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=data)
        output_b = cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,metrics.nbytes)
        kernel(q,(args.size,),None,input_b,output_b,np.int32(args.size),np.int32(args.burn),np.int32(args.window))
        cl.enqueue_copy(q,metrics,output_b).wait()
        flags = np.flatnonzero(np.min(metrics[:4],axis=0)>.015)
        tmp = target.with_suffix('.tmp.npz')
        np.savez(tmp,metrics=metrics,kind=kind.astype('u1'),flags=flags.astype('i4'),geometry=args.geometry,settings=np.array([args.size,args.burn,args.window,SEED]))
        os.replace(tmp,target)
        dt = time.perf_counter()-start
        times.append(dt)
        print(batch,args.size,len(flags),round(dt,3),flush=True)
    if times:
        print(json.dumps(dict(batches=len(times),median_seconds=float(np.median(times)))),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--size',type=int,default=4096)
    p.add_argument('--first',type=int,default=0)
    p.add_argument('--batches',type=int,default=128)
    p.add_argument('--burn',type=int,default=512)
    p.add_argument('--window',type=int,default=192)
    p.add_argument('--gpu',type=int,default=0)
    p.add_argument('--geometry',choices=('single','linked','smooth'),default='single')
    screen(p.parse_args())
