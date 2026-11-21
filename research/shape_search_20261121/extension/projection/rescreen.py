"""GPU deep rescreen, rotation/scale invariant ranking against 31 earlier cases."""
import argparse
import json
from pathlib import Path
import numpy as np
import pyopencl as cl
from survey import transformed_kernel
from select_shapes import descriptors


def fingerprint(z):
    # Equal-length consecutive samples; retain PCA variance, Radial and chord
    # quantiles, projection occupancy, local covariance and spatial extent.
    return descriptors(z)


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    dev=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    ctx=cl.Context([dev]);q=cl.CommandQueue(ctx);kernel=transformed_kernel(ctx,True)
    ids=[];pars=[];anchors=[];sets=[];all_deep=[];bases=[];short_flags=0;short_total=0
    for path in sorted(args.screen.glob('screen-*.npz')):
        with np.load(path) as f:
            if int(path.stem.split('-')[1])>=args.batches:continue
            ids.extend(f['ids']);pars.append(f['data']);anchors.extend(f['parents'])
            sets.append(f['hashes']);a=f['scores'];short_flags+=int(np.sum(a>.01));short_total+=len(a)
    data=np.concatenate(pars,axis=1);ids=np.array(ids);anchors=np.array(anchors)
    all_keys=np.concatenate(sets);unique=len(np.unique(all_keys))
    if unique!=short_total:raise ValueError('Duplicate parameter storage')
    candidate=[];features=[];spectra=[]
    for lo in range(0,len(ids),640):
        d=np.ascontiguousarray(data[:,lo:lo+640,:]);n=d.shape[1];tail_n=1536
        m=np.empty((20,n),'f4');z=np.empty((tail_n,n,4),'f4')
        ib=cl.Buffer(ctx,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=d)
        ob=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,m.nbytes);zb=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,z.nbytes)
        kernel(q,(n,),None,ib,ob,np.int32(n),np.int32(24576),np.int32(4096),zb,np.int32(tail_n))
        cl.enqueue_copy(q,m,ob).wait();cl.enqueue_copy(q,z,zb).wait()
        rates=np.stack((m[:4].T,m[12:16].T,m[16:20].T),axis=2)
        score=np.nan_to_num(rates.min(axis=(1,2)),nan=-999.)
        all_deep.append(score)
        for j in np.flatnonzero(score>.005):
            feat=fingerprint(z[:,j,:3])
            if np.isfinite(feat).all():
                candidate.append(lo+j);features.append(feat);spectra.append(rates[j])
    old=json.loads((args.study/'three'/'manifest.json').read_text())
    for row in old:
        with np.load(row['cloud']) as f:z=f['points']
        for shift in (0,32000,64000,128000):
            bases.append(fingerprint(z[shift:shift+1536]))
    f=np.array(features);ref=np.array(bases);pool=np.vstack((f,ref))
    spread=np.maximum(np.quantile(pool,.9,axis=0)-np.quantile(pool,.1,axis=0),.05)
    # Downweight sparse histogram bins and local random-neighbor estimates.
    spread[14:29]*=1.7;spread[29:]*=1.6
    normalized=f/spread;controls=ref/spread
    dist=np.sqrt(np.min(np.sum((normalized[:,None,:]-controls[None,:,:])**2,axis=2),axis=1)/f.shape[1])
    ix=np.array(candidate);qr=np.array(spectra);score=qr.min(axis=(1,2))
    ratio=data[2:5,ix,3].max(0)/data[2:5,ix,3].min(0)
    chosen=[int(np.argmax(score)),int(np.argmax(dist)),int(np.argmax(ratio))]
    novelty=np.minimum(dist,np.min([np.linalg.norm(normalized-normalized[j],axis=1)/np.sqrt(f.shape[1]) for j in chosen],axis=0))
    while len(chosen)<min(42,len(f)):
        novelty[chosen]=-1
        j=int(np.argmax(novelty*np.clip(score/.02,.5,1.6)));chosen.append(j)
        novelty=np.minimum(novelty,np.linalg.norm(normalized-normalized[j],axis=1)/np.sqrt(f.shape[1]))
    target=args.out/'shortlist.npz';k=ix[chosen]
    np.savez_compressed(target,ids=ids[k],raw=data[:,k,:],parents=anchors[k],rates=qr[chosen],distances=dist[chosen],ratio=ratio[chosen],spread=spread,refs=ref)
    stats=dict(unique=unique,short_above_001=short_flags,deep_tested=len(ids),deep_above_0005=int(sum(np.sum(x>.005) for x in all_deep)),ranked=len(f),shortlist=len(chosen))
    (args.out/'stage.json').write_text(json.dumps(stats,indent=2)+'\n',encoding='utf8')
    print(stats,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',type=Path,required=True);p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--batches',type=int,default=384)
    main(p.parse_args())
