"""Recover attracting periodic companions signaled by a multistart orbit."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from model import replay
from reproduce import raw_parameters
from cycle_seeds import orbit_eval


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    for file in sorted(args.source.glob('case-*.json')):
        rec=json.loads(file.read_text())
        detected=[(i,e) for i,e in enumerate(rec['cohort']) if e['period'] is not None]
        if not detected: continue
        dest=args.out/f'seeds-{rec["case"]}.json'
        if dest.exists(): continue
        raw=raw_parameters(rec)
        with np.load(rec['cloud']) as f: z=f['points'].astype(float)
        rad=np.sqrt(np.sum((z-z.mean(0))**2,axis=1).mean());tree=cKDTree(z)
        seeds=[]
        for i,e in detected:
            p=e['period']
            tail,*_=replay(raw,np.asarray(e['initial']),rec['cohort_burn'],rec['cohort_window_length'],4,max(2*p,64),False)
            state=tail[-1];y,j,orbit=orbit_eval(state,raw,p)
            existing=next((s for s in seeds if s['period']==p and np.min(np.linalg.norm(np.asarray(s['orbit'])-state,axis=1))<1e-7),None)
            if existing:
                existing['cohort_indices'].append(i)
                continue
            vals=np.linalg.eigvals(j)
            seeds.append(dict(period=p,seed=state.tolist(),orbit=orbit.tolist(),numerical_unstable=int(np.sum(abs(vals)>1)),max_orbit_distance_to_sample_rms=float(tree.query(orbit)[0].max()/rad),cohort_indices=[i]))
        dest.write_text(json.dumps(dict(case=rec['case'],cycles=seeds),indent=2)+'\n',encoding='utf8')
        print(rec['case'],[(s['period'],s['cohort_indices']) for s in seeds],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
