"""Resolve long near-returns at the already sampled e immediately before chaotic windows."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from band_continuation import dynamics
from field import integrate,ring_tangent
from orbits import refine


def candidates(z,dt,base_t,rms):
    sols=[]
    for axis in np.argsort(np.ptp(z,axis=0))[::-1]:
        level=np.median(z[:,axis])
        ids=np.flatnonzero((z[:-1,axis]<level)&(z[1:,axis]>=level))
        a=(level-z[ids,axis])/(z[ids+1,axis]-z[ids,axis])
        cross=z[ids]+a[:,None]*(z[ids+1]-z[ids]);times=(ids+a)*dt
        for lag in range(1,min(6000,len(times)-50)):
            t=np.mean(times[-40:]-times[-lag-40:-lag])
            if t>base_t*66:break
            if min(abs(t-n*base_t)/(n*base_t) for n in (8,16,32,64))>.015:continue
            error=np.max(np.linalg.norm(cross[-40:]-cross[-lag-40:-lag],axis=1))/rms
            if error<.006:sols.append((t,error,int(axis),cross[-1-lag].copy()))
    return sorted(sols,key=lambda v:(round(v[0]/base_t),v[1]))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',required=True,type=Path)
    ap.add_argument('--branch',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    m=json.loads((a.packet/'manifest.json').read_text())
    p=np.asarray(next(r for r in m['stable_cycles'] if r['case']==28756)['parameters'])
    row=next(r for r in json.loads(a.branch.read_text())['trials'] if r['fraction']==.95)
    p[8]=row['e'];rms=.891
    results=[];all_start=[]
    for index,z0 in enumerate([row['initial'],[-2.61,-.72,1.515],[1.837,-2.08,-.72]]):
        step=.008 if index==0 else .00625
        rates,z=integrate(p,np.array(z0),step,round(15000/step),round(750/step),8,4,ring_tangent)
        guesses=candidates(z,step*4,row['primary_period'],rms)
        all_start.append(dict(initial=z0,step=step,rates=rates.tolist(),best_returns=[g[:3] for g in guesses[:6]]))
        if not results:
            for t,err,axis,loc in guesses[:6]:
                got=refine(axis,loc,t,rms,p,1)
                if got is None:continue
                seed,T=got;stats,sol=dynamics(seed,T,p,True)
                if stats['closure']/rms>1e-5 or stats['modulus']>=.995:continue
                phase=np.linspace(0,T,19,endpoint=False)
                gaps=[np.max(np.linalg.norm(sol.sol(phase)[:3]-sol.sol((phase+T/k)%T)[:3],axis=0))/rms for k in range(2,65)]
                if min(gaps)<.001:continue
                points=sol.sol(np.linspace(0,T,65536,endpoint=False))[:3].T
                np.savez_compressed(a.out/'child.npz',points=points)
                results.append(dict(**stats,period_ratio=T/row['primary_period'],refined_start=seed.tolist(),
                               min_phase_gap=float(min(gaps)),axis=axis,ranked_seed_error=err))
                break
        summary=dict(fraction=.95,e=float(p[8]),primary_period=row['primary_period'],parameters=p.tolist(),
                initial_condition_runs=all_start,accepted_cycles=results)
        (a.out/'report.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(index,len(guesses),[(v['period_ratio'],v['dominant']) for v in results],flush=True)


if __name__=='__main__':
    main()
