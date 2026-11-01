"""Three-row kernels with comparable row masses, avoiding nearly absent rows."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import set_num_threads
from survey_dual import survey

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--count',type=int,default=1024)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4);rng=np.random.default_rng(201212901)
    Ls=[];zs=[];ws=[];rows=[]
    for case in range(args.count):
        contrast=float(rng.choice([1.,2.,4.,8.,12.]));odds=np.array([rng.uniform(-2,2),-contrast*rng.uniform(.2,1.2),contrast*rng.uniform(.2,1.2)])
        masses=rng.normal(scale=.5,size=3);masses-=np.log(np.sum(np.exp(masses)))
        L=np.c_[-np.logaddexp(0,odds),-np.logaddexp(0,-odds)]+masses[:,None]
        while True:
            wx=float(rng.choice([2.002,2.005,2.01,2.02,2.04,2.06,2.1,2.15,2.2]))
            wy=float(rng.choice([1.05,1.1,1.2,1.35,1.5,1.65,1.8,1.9]))
            if abs((1-wx)**2*(1-wy))<.98:break
        amplitude=float(rng.choice([.001,.1,1.,5.,20.]));z=rng.normal(size=3)*amplitude
        Ls.append(L);zs.append(z);ws.append([wx,wy]);rows.append(dict(case=case,shape=[3,2],contrast=contrast,wx=wx,wy=wy,log_kernel=L.tolist(),initial_state=z.tolist(),log_odds=odds.tolist(),log_row_masses=masses.tolist(),determinant=(1-wx)**2*(1-wy)))
    t=time.perf_counter();stats,ends,tails=survey(np.array(Ls),np.array(zs),np.array(ws),10000,30000,2048)
    for r,s in zip(rows,stats):
        status=int(s[4]);r.update(status=['completed','numerical capture','escape threshold'][status],completed_steps=int(s[5]),max_abs_state=float(s[3]),final_abs_state=float(s[6]),rates=None if status else dict(mean=float(s[0]),first=float(s[1]),last=float(s[2])))
    result=dict(seed=201212901,shape=[3,2],warm=10000,observed=30000,escape_threshold=1e5,capture_threshold=1e-11,runs=rows,elapsed_seconds=time.perf_counter()-t,
        sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'),Path('survey_dual.py'))})
    (args.output/'3x2.json').write_text(json.dumps(result,indent=2)+'\n');np.savez_compressed(args.output/'3x2-endpoints-and-last2048.npz',endpoints=ends,last2048=tails)
    flags=[r for r in rows if r['rates'] and r['rates']['mean']>.005];print(json.dumps(dict(cases=len(rows),counts={s:sum(r['status']==s for r in rows) for s in ('completed','numerical capture','escape threshold')},flags=[dict(case=r['case'],wx=r['wx'],wy=r['wy'],rate=r['rates']['mean']) for r in flags],elapsed=result['elapsed_seconds'])),flush=True)
