"""Rescreen each targeted batch's strongest third rates on GPU and independent CPU."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pyopencl as cl

from model import replay
from third_growth import kernel_for,run_gpu


def read_candidates(paths):
    arrays=[]
    for path in paths:
        batch=int(path.stem.rsplit('-',1)[1])
        with np.load(path) as f:
            arrays.append((f['index']+batch*int(f['settings'][0]),f['data'],f['parent'],f['kind'],f['all_min'],f['rates']))
    ids,raw,anchors,kinds,short,windows=zip(*arrays)
    return np.concatenate(ids),np.concatenate(raw,axis=1),np.concatenate(anchors),np.concatenate(kinds),np.concatenate(short),np.concatenate(windows)


def study(args):
    args.out.mkdir(parents=True,exist_ok=True)
    args.records.mkdir(parents=True,exist_ok=True)
    device=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    ctx=cl.Context([device]);queue=cl.CommandQueue(ctx)
    for geom in args.geometries:
        chunks=sorted(args.survey.glob(geom+'-*.npz'))
        if len(chunks)!=args.batches:raise ValueError(f'Expected {args.batches} {geom} batches, found {len(chunks)}')
        ids,raw,parents,kinds,short,short_windows=read_candidates(chunks)
        kernel=kernel_for(ctx,geom)
        rates=run_gpu(queue,kernel,raw,geom,args.burn,args.window)
        scores=np.nan_to_num(rates.min(axis=(1,2)),nan=-999.)
        cand=np.unique(np.r_[np.flatnonzero(scores>.001),np.argsort(scores)[-16:]])
        path=args.out/f'deep-{geom}.npz'
        np.savez_compressed(path,ids=ids,data=raw,parents=parents,kind=kinds,short_windows=short_windows,deep_windows=rates,chosen=cand,settings=[args.burn,args.window])
        ranked=[]
        for i in sorted(cand,key=lambda i:-scores[i]):
            row=raw[:,i,:].astype(float)
            points,long_windows,_,_,_=replay(row,row[1,:3].copy(),args.cpu_burn,args.cpu_window,4,4096,True)
            ranked.append(dict(index=int(ids[i]),parent=int(parents[i]),gpu_windows=rates[i].tolist(),cpu_windows=long_windows.tolist(),min_cpu_third=float(long_windows[:,2].min()),pca_fraction=(np.linalg.eigvalsh(np.cov(points.T))/np.var(points,axis=0,ddof=1).sum()).tolist()))
        doc=dict(geometry=geom,device=device.name,focus='Parameter perturbations near 14 parent spectra with highest third exponent.',short_burn=1024,short_window=384,short_tested=len(short),short_all_above_0006=int(np.sum(short>.006)),gpu_deep_tested=len(rates),gpu_deep_burn=args.burn,gpu_deep_window=args.window,gpu_deep_all_above_0001=int(np.sum(scores>.001)),cpu_tested=len(ranked),cpu_burn=args.cpu_burn,cpu_window=args.cpu_window,cpu_all_above_0001=int(sum(np.min(v['cpu_windows'])>.001 for v in ranked)),deep_cache_relative='projection/three-deep/'+path.name,deep_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),cases=ranked)
        (args.records/f'three-rate-{geom}.json').write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
        print(geom,len(short),int(np.sum(short>.006)),len(rates),float(scores.max()),[(v['index'],round(v['min_cpu_third'],5)) for v in ranked[:8]],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--survey',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--records',type=Path,required=True)
    p.add_argument('--batches',type=int,default=96)
    p.add_argument('--geometries',nargs='+',choices=('single','linked','smooth'),default=['single','smooth','linked'])
    p.add_argument('--burn',type=int,default=32768)
    p.add_argument('--window',type=int,default=8192)
    p.add_argument('--cpu-burn',type=int,default=100000)
    p.add_argument('--cpu-window',type=int,default=65536)
    study(p.parse_args())
