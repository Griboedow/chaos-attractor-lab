"""Declared rectangular Sinkhorn screen in three and four state dimensions."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import set_num_threads
import model
from base_survey import survey

def inputs(m,n,count):
    seed=201212600+10*m+n;rng=np.random.default_rng(seed);matrices=[];starts=[];omegas=[];rows=[]
    for case in range(count):
        contrast=float(rng.choice([1.,2.,5.,10.,20.,40.,80.]));family=case%4
        if family==0:L=rng.normal(size=(m,n))*contrast
        elif family==1:L=np.eye(m,n)*contrast+rng.normal(scale=.1,size=(m,n))
        elif family==2:
            a=rng.normal(size=(m,2));b=rng.normal(size=(n,2));L=-contrast*np.sum((a[:,None]-b[None,:])**2,axis=2)
        else:L=rng.uniform(-contrast,contrast,size=(m,n))
        L-=L.max();w=float(rng.choice([1.,1.5,1.8,1.9,1.95,1.98,1.99,1.995,1.999]))
        amplitude=float(rng.choice([.01,.3,3.,30.,300.]));z=rng.normal(size=m+n-2)*amplitude
        la,lb=model.targets(L);matrices.append(L);starts.append(z);omegas.append(w)
        rows.append(dict(case=case,shape=[m,n],family=family,contrast=contrast,omega=w,start_scale=amplitude,
            log_kernel=L.tolist(),initial_state=z.tolist(),coordinate_box=model.bounds(L,la,lb,w).tolist()))
    return np.array(matrices),np.array(starts),np.array(omegas),rows,seed

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--count',type=int,default=384);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    set_num_threads(4);warm=20000;observed=40000;keep=2048
    for m,n in [(2,3),(3,2),(2,4),(4,2)]:
        L,z,w,rows,seed=inputs(m,n,args.count);t=time.perf_counter();stats,ends,tails=survey(L,z,w,warm,observed,keep)
        names=['finite_rate','first_half_rate','last_half_rate','max_abs_state','final_abs_state','capture_step','completed_steps','final_step_residual']
        for r,values in zip(rows,stats):
            r.update(dict(zip(names,[float(x) for x in values])))
            r['growth_status']='Not assessed after early numerical capture' if r['capture_step'] else 'Finite-window estimate'
            r['tail_samples_retained']=max(0,int(r['completed_steps'])-(warm+observed-keep))
        result=dict(shape=[m,n],seed=seed,warm=warm,observed=observed,keep=keep,elapsed_seconds=time.perf_counter()-t,runs=rows,
            sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('base_survey.py'),Path('model.py'))})
        (args.output/f'{m}x{n}.json').write_text(json.dumps(result,indent=2)+'\n')
        np.savez_compressed(args.output/f'{m}x{n}-endpoints-and-last2048.npz',endpoints=ends,last2048=tails)
        flags=[r for r in rows if not r['capture_step'] and r['finite_rate']>1e-5]
        print(json.dumps(dict(shape=[m,n],cases=len(rows),captures=sum(bool(r['capture_step']) for r in rows),flags=[dict(case=r['case'],rate=r['finite_rate'],last=r['last_half_rate'],omega=r['omega']) for r in flags],elapsed=result['elapsed_seconds'])),flush=True)
