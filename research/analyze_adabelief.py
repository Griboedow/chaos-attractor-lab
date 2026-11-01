"""Saved-evidence analysis and targeted follow-ups; never overwrites input data."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from adabelief_map import step,jacobian
from check_adabelief import run,brief


def followups(source,output):
    saved=np.load(source/'neighborhood/run-022.npz')
    reports=[]
    for i,(z,b,c) in enumerate([(saved['end'],.522,.15),([.73,0,0],.53,.21)]):
        data=run(z,b,c,steps=2000000,keep=200000)
        np.savez_compressed(output/f'run-{i:03d}.npz',**data)
        reports.append(brief(data));print(json.dumps(reports[-1]),flush=True)
    b,c=.53,.21;k=(1-b)/(1+b);s=(k/2-1e-14)**2
    a=np.sqrt(s-1e-12/(1-c))/(2*b/(1+b));z=np.array([a,-k*a,s])
    multipliers=np.linalg.eigvals(jacobian(step(z,b,c),b,c)@jacobian(z,b,c))
    return dict(runs=reports,cycle_state=z.tolist(),cycle_multiplier_moduli=abs(multipliers).tolist(),
                cycle_residual=float(abs(step(step(z,b,c),b,c)-z).max()),
                continuation_source='neighborhood/run-022.npz; new 20000-step tangent warmup')


def analysis(source,output):
    report={}
    for name in ['ensemble','neighborhood','epsilon','long']:
        data=json.loads((source/name/'summary.json').read_text())['runs']
        rows=[]
        for r in data:
            rows.append(dict(id=r['id'],parameters=r['parameters'],leading=r['spectrum'][0],
                             minimum_window=r['window_min'][0]))
        report[name]=rows
    orbit=np.load(source/'long/run-000.npz')['tail']
    scale=np.ptp(orbit,axis=0)
    lag_rms=[]
    # Compare the last 8192 states to shifted states, in fixed full-state units.
    for lag in range(1,4097):
        difference=(orbit[-8192:]-orbit[-8192-lag:-lag])/scale
        lag_rms.append(float(np.sqrt(np.mean(difference**2))))
    report['period_screen']=dict(maximum_lag=4096,comparison_states=8192,
       best_lag=int(np.argmin(lag_rms)+1),minimum_rms=min(lag_rms),scale=scale.tolist())
    np.savez_compressed(output/'period-screen.npz',lags=np.arange(1,4097),rms=lag_rms)
    tree=cKDTree(orbit[::2]/scale)
    control=tree.query(orbit[1::2]/scale)[0]
    support=[]
    for i in range(32):
        data=np.load(source/f'ensemble/run-{i:03d}.npz')
        distance=tree.query(data['tail']/scale)[0]
        support.append(dict(id=i,median=float(np.median(distance)),p95=float(np.quantile(distance,.95)),
                            maximum=float(distance.max())))
    report['support_comparison']=dict(baseline_median=float(np.median(control)),
        baseline_p95=float(np.quantile(control,.95)),scale=scale.tolist(),runs=support,
        note='Distances to sampled reference support, not distances to an exact invariant set or basin probabilities.')
    return report


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--followups',action='store_true')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    report=followups(a.source,a.output) if a.followups else analysis(a.source,a.output)
    (a.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['ensemble','neighborhood','epsilon','long','support_comparison']},indent=2))


if __name__=='__main__': main()
