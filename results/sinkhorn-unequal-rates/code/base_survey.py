"""Declared fixed-relaxation matrix-scaling screen; all starts are recorded."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit,prange,set_num_threads
import model

@njit(cache=True,parallel=True)
def survey(matrices,starts,omegas,warm,observed,keep):
    count=len(starts);d=starts.shape[1]
    stats=np.zeros((count,8));ends=np.zeros_like(starts);tails=np.full((count,keep,d),np.nan)
    for lane in prange(count):
        L=matrices[lane];la,lb=model.targets(L);z=starts[lane].copy();w=omegas[lane]
        v=np.ones(d)/np.sqrt(d);total=0.;first=0.;last=0.;maximum=0.;captured=0
        for iteration in range(warm+observed):
            z,v=model.step(z,v,L,la,lb,w)
            norm=np.linalg.norm(v)
            if norm>1e-300:v/=norm
            else:v=np.ones(d)/np.sqrt(d)
            maximum=max(maximum,np.max(np.abs(z)))
            if iteration>=warm:
                growth=np.log(max(norm,1e-300));total+=growth
                if iteration<warm+observed//2:first+=growth
                else:last+=growth
            if iteration>=warm+observed-keep:tails[lane,iteration-(warm+observed-keep)]=z
            if np.max(np.abs(z))<1e-11:
                captured=iteration+1;break
        ends[lane]=z
        stats[lane]=np.array([total/observed,first/(observed//2),last/(observed-observed//2),maximum,np.max(np.abs(z)),captured,iteration+1,np.max(np.abs(z-model.step(z,v,L,la,lb,w)[0]))])
    return stats,ends,tails

def inputs(n,count):
    rng=np.random.default_rng(201212564+n)
    matrices=[];starts=[];omegas=[];records=[]
    for i in range(count):
        contrast=float(rng.choice([.5,2.,5.,10.,20.,40.]))
        kind=i%4
        if kind==0:L=rng.normal(size=(n,n))*contrast
        elif kind==1:
            L=np.eye(n)*contrast+rng.normal(scale=.1,size=(n,n))
        elif kind==2:
            points=rng.normal(size=(n,2));L=-contrast*np.sum((points[:,None]-points[None,:])**2,axis=2)
        else:L=rng.uniform(-contrast,contrast,size=(n,n))
        L-=np.max(L)
        omega=float(rng.choice([1.,1.5,1.8,1.9,1.95,1.98,1.99,1.995,1.999]))
        amplitude=float(rng.choice([.01,.3,3.,30.,300.]))
        z=rng.normal(size=2*n-2)*amplitude
        matrices.append(L);starts.append(z);omegas.append(omega)
        la,lb=model.targets(L)
        # log cross ratio, evaluated without products or overflowing exponentials.
        cross=max(float(np.max(L[i1]-L[i2])-np.min(L[i1]-L[i2])) for i1 in range(n) for i2 in range(n))
        contraction=float(np.tanh(cross/4))
        records.append(dict(case=i,shape=[n,n],family=kind,contrast=contrast,start_scale=amplitude,
            omega=omega,log_kernel=L.tolist(),initial_state=z.tolist(),
            sufficient_global_omega_upper=2/(1+contraction),
            coordinate_box=model.bounds(L,la,lb,omega).tolist()))
    return np.array(matrices),np.array(starts),np.array(omegas),records

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--count',type=int,default=192)
    ap.add_argument('--warm',type=int,default=10000);ap.add_argument('--observed',type=int,default=20000)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    for n in [2,3,4]:
        matrices,starts,omegas,records=inputs(n,args.count)
        start=time.perf_counter();stats,ends,tails=survey(matrices,starts,omegas,args.warm,args.observed,512)
        names=['finite_rate','first_half_rate','last_half_rate','max_abs_state','final_abs_state','capture_step','completed_steps','final_step_residual']
        for record,values in zip(records,stats):
            record.update({key:float(value) for key,value in zip(names,values)})
            record['growth_status']='not measured after early fixed-point capture' if record['capture_step'] else 'finite-window estimate'
            record['tail_samples_retained']=max(0,int(record['completed_steps'])-(args.warm+args.observed-512))
        result=dict(size=n,warm=args.warm,observed=args.observed,keep=512,runs=records,
            elapsed_seconds=time.perf_counter()-start,
            source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(model.__file__)]})
        (args.output/f'n{n}.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        np.savez_compressed(args.output/f'n{n}-endpoints-and-last512.npz',endpoints=ends,last512=tails)
        remaining=[r for r in records if not r['capture_step']]
        print(json.dumps(dict(size=n,cases=len(records),captured=len(records)-len(remaining),positive=[dict(case=r['case'],rate=r['finite_rate'],last=r['last_half_rate'],norm=r['final_abs_state'],omega=r['omega']) for r in remaining if r['finite_rate']>1e-5],elapsed=result['elapsed_seconds'])),flush=True)
