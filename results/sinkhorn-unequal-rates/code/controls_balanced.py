"""Finite basin and parameter controls for the two selected filament cases."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import set_num_threads
from survey_dual import survey

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--pilot',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4);data=json.loads(args.input.read_text());pilot=json.loads(args.pilot.read_text());rng=np.random.default_rng(201212902)
    rows=[];Ls=[];zs=[];ws=[]
    for case in (150,987):
        r=next(r for r in data['runs'] if r['case']==case);p=next(r for r in pilot['runs'] if r['case']==case);end=np.array(r['blocks'][-1]['checkpoint']);original=np.array(p['initial_state']);L=np.array(r['log_kernel']);base=np.array([r['wx'],r['wy']])
        jobs=[('original',original,base),('ordinary-original',original,np.ones(2))]
        for scale in (1e-6,1e-3,.03,.3):
            for i in range(8):jobs.append((f'endpoint-offset-{scale}-{i}',end+scale*rng.normal(size=3),base))
        for scale in (1e-6,.1,1.,5.):
            for i in range(4):jobs.append((f'fresh-scale-{scale}-{i}',scale*rng.normal(size=3),base))
        for dx in (-.001,0.,.001):
            for dy in (-.001,0.,.001):
                if dx or dy:jobs.append((f'rate-offset-{dx}-{dy}',end,base+np.array([dx,dy])))
        for label,z,w in jobs:
            rows.append(dict(case=case,label=label,log_kernel=L.tolist(),initial_state=z.tolist(),wx=float(w[0]),wy=float(w[1])));Ls.append(L);zs.append(z);ws.append(w)
    t=time.perf_counter();stats,ends,tails=survey(np.array(Ls),np.array(zs),np.array(ws),20000,200000,2048)
    for r,s,end in zip(rows,stats,ends):
        status=int(s[4]);r.update(status=['completed','numerical capture','escape threshold'][status],completed_steps=int(s[5]),max_abs_state=float(s[3]),final_abs_state=float(s[6]),endpoint=end.tolist(),rates=None if status else dict(mean=float(s[0]),first=float(s[1]),last=float(s[2])))
    result=dict(runs=rows,warm=20000,observed=200000,seed=201212902,elapsed_seconds=time.perf_counter()-t,parent_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),pilot_sha256=hashlib.sha256(args.pilot.read_bytes()).hexdigest(),sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'),Path('survey_dual.py'))})
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');np.savez_compressed(args.output/'last2048-states.npz',states=tails)
    print(json.dumps([dict(case=c,counts={s:sum(r['status']==s for r in rows if r['case']==c) for s in ('completed','numerical capture','escape threshold')},negative_completed=[r['label'] for r in rows if r['case']==c and r['rates'] and r['rates']['mean']<=0]) for c in (150,987)]),flush=True)
