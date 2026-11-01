"""Finite screen of a declared unequal-rate extension, including escaping runs."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import njit,prange,set_num_threads
import model,model_dual

@njit(cache=True,parallel=True)
def survey(Ls,zs,ws,warm,observed,keep):
    count,d=zs.shape;stats=np.zeros((count,7));ends=np.zeros_like(zs);tails=np.full((count,keep,d),np.nan)
    for lane in prange(count):
        L=Ls[lane];la,lb=model.targets(L);z=zs[lane].copy();wx,wy=ws[lane];v=np.ones(d)/np.sqrt(d)
        total=0.;first=0.;last=0.;status=0;maximum=0.
        for k in range(warm+observed):
            z,v=model_dual.step(z,v,L,la,lb,wx,wy);norm=np.linalg.norm(v)
            if norm>1e-250:v/=norm
            else:v=np.ones(d)/np.sqrt(d)
            maximum=max(maximum,np.max(np.abs(z)))
            if not np.isfinite(z).all() or maximum>1e5:status=2;break
            if np.max(np.abs(z))<1e-11:status=1;break
            if k>=warm:
                val=np.log(max(norm,1e-300));total+=val
                if k<warm+observed//2:first+=val
                else:last+=val
            if k>=warm+observed-keep:tails[lane,k-warm-observed+keep]=z
        ends[lane]=z;stats[lane]=np.array([total/observed,first/(observed//2),last/(observed//2),maximum,status,k+1,np.max(np.abs(z))])
    return stats,ends,tails

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--count',type=int,default=768)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    for m,n in [(2,2),(2,3),(3,2)]:
        seed=201212700+10*m+n;rng=np.random.default_rng(seed);Ls=[];zs=[];ws=[];rows=[]
        for case in range(args.count):
            c=float(rng.choice([1.,2.,5.,10.,20.,40.]));L=rng.normal(size=(m,n))*c;L-=L.max()
            while True:
                wx=float(rng.choice([1.1,1.5,1.8,1.95,2.02,2.1,2.3,2.6,3.,4.]))
                wy=float(rng.choice([.1,.3,.5,.7,.9,1.1,1.5,1.8,1.95,1.99]))
                if abs((1-wx)**(m-1)*(1-wy)**(n-1))<.98:break
            if case%2:wx,wy=wy,wx
            amplitude=float(rng.choice([.001,.1,1.,5.,20.]));z=rng.normal(size=m+n-2)*amplitude
            Ls.append(L);zs.append(z);ws.append([wx,wy]);rows.append(dict(case=case,shape=[m,n],contrast=c,wx=wx,wy=wy,log_kernel=L.tolist(),initial_state=z.tolist(),start_scale=amplitude,determinant=(1-wx)**(m-1)*(1-wy)**(n-1)))
        t=time.perf_counter();stats,ends,tails=survey(np.array(Ls),np.array(zs),np.array(ws),10000,20000,2048)
        for r,s in zip(rows,stats):
            status=int(s[4]);r.update(status=['completed','numerical capture','escape threshold'][status],completed_steps=int(s[5]),max_abs_state=float(s[3]),final_abs_state=float(s[6]),
                rates=None if status else dict(mean=float(s[0]),first=float(s[1]),last=float(s[2])))
        result=dict(seed=seed,shape=[m,n],warm=10000,observed=20000,escape_threshold=1e5,capture_threshold=1e-11,runs=rows,elapsed_seconds=time.perf_counter()-t,
            sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'))})
        (args.output/f'{m}x{n}.json').write_text(json.dumps(result,indent=2)+'\n')
        np.savez_compressed(args.output/f'{m}x{n}-endpoints-and-last2048.npz',endpoints=ends,last2048=tails)
        print(json.dumps(dict(shape=[m,n],counts={s:sum(r['status']==s for r in rows) for s in ('completed','numerical capture','escape threshold')},flags=[dict(case=r['case'],wx=r['wx'],wy=r['wy'],rates=r['rates'],radius=r['final_abs_state']) for r in rows if r['rates'] and r['rates']['mean']>.005],elapsed=result['elapsed_seconds'])),flush=True)
