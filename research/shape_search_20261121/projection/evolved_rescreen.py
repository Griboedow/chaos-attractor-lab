"""Confirm third-generation growth and keep varied folded spatial distributions."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pyopencl as cl
from third_growth import kernel_for,run_gpu
from model import replay
from select_shapes import descriptors


def study(args):
    paths=sorted(args.survey.glob('evolved-*.npz'))
    all_ids=[];raw=[];parents=[];kinds=[];added=[];short_counts=np.zeros(3,dtype=int)
    for path in paths:
        with np.load(path) as f:
            batch=int(path.stem.split('-')[1]);n=int(f['settings'][0])
            all_ids.append(f['index']+batch*n);raw.append(f['data'])
            parents.append(f['parent']);kinds.append(f['kind']);added.append(f['added'])
            short_counts += [n,np.sum(f['all_min']>.025),np.sum(f['all_min']>.05)]
    ids=np.concatenate(all_ids);raw=np.concatenate(raw,axis=1)
    parents=np.concatenate(parents);kinds=np.concatenate(kinds);added=np.concatenate(added)
    dev=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    ctx=cl.Context([dev]);queue=cl.CommandQueue(ctx);kernel=kernel_for(ctx,'single')
    # Bounded work chunks keep the Windows display watchdog responsive.
    deep=np.concatenate([run_gpu(queue,kernel,raw[:,i:i+1536,:],'single',32768,8192) for i in range(0,len(ids),1536)])
    score=np.nan_to_num(deep.min(axis=(1,2)),nan=-999.)
    best=[]
    for n in range(5):
        group=np.flatnonzero(added==n)
        best.extend(group[np.argsort(score[group])[-26:]])
    rng=np.random.default_rng(406119)
    more=np.setdiff1d(np.flatnonzero(score>.015),best)
    best=np.unique(np.r_[best,rng.choice(more,size=min(54,len(more)),replace=False)])
    cases=[];features=[];eligible=[]
    for i in sorted(best,key=lambda i:-score[i]):
        data=raw[:,i,:].astype(float)
        cloud,rates,_,_,_=replay(data,data[1,:3].copy(),100000,65536,4,8192,True)
        entry=dict(index=int(ids[i]),parent=int(parents[i]),added=int(added[i]),gpu_min_third=float(deep[i,:,2].min()),cpu_windows=rates.tolist(),min_cpu_third=float(rates[:,2].min()))
        cases.append(entry)
        if rates.min()>.005:
            eligible.append(len(cases)-1)
            features.append(descriptors(cloud))
    f=np.asarray(features)
    u=np.clip((f-np.median(f,0))/np.maximum(np.quantile(f,.9,axis=0)-np.quantile(f,.1,axis=0),.04),-3,3)
    chosen=[]
    for n in range(5):
        group=[i for i,j in enumerate(eligible) if cases[j]['added']==n]
        if group:chosen.append(max(group,key=lambda i:cases[eligible[i]]['min_cpu_third']))
    nearest=np.min([np.linalg.norm(u-u[k],axis=1) for k in chosen],axis=0)
    while len(chosen)<min(18,len(eligible)):
        nearest[chosen]=-1
        k=int(np.argmax(nearest));chosen.append(k)
        nearest=np.minimum(nearest,np.linalg.norm(u-u[k],axis=1))
    chosen_ids={eligible[i] for i in chosen}
    for i,row in enumerate(cases):row['chosen']=i in chosen_ids
    args.out.mkdir(parents=True,exist_ok=True)
    deep_path=args.out/'deep-evolved.npz'
    np.savez_compressed(deep_path,ids=ids,data=raw,parents=parents,kind=kinds,added=added,deep_windows=deep)
    # The full control table is kept as one small JSON per 48 cases.
    detail=args.measurements/'three'/'evolved-audit'
    detail.mkdir(parents=True,exist_ok=True)
    for j in range(0,len(cases),48):
        (detail/f'cpu-{j//48:02d}.json').write_text(json.dumps(cases[j:j+48],indent=2)+'\n',encoding='utf8')
    doc=dict(geometry='evolved',short_tested=int(short_counts[0]),short_all_above_0025=int(short_counts[1]),short_all_above_005=int(short_counts[2]),gpu_deep_tested=len(ids),gpu_deep_burn=32768,gpu_deep_window=8192,gpu_deep_all_above_0005=int((score>.005).sum()),cpu_tested=len(cases),cpu_burn=100000,cpu_window=65536,cpu_all_above_0005=len(eligible),deep_sha256=hashlib.sha256(deep_path.read_bytes()).hexdigest(),deep_cache_relative='projection/three-deep/deep-evolved.npz',cases=[row for row in cases if row['chosen']])
    (args.measurements/'three-rate-evolved.json').write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
    print(len(ids),len(eligible),[(r['index'],round(r['min_cpu_third'],3),r['added']) for r in doc['cases']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--survey',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--measurements',type=Path,required=True)
    study(p.parse_args())
