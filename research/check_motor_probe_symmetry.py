"""Numerical weighted permutation closures of the distinct five-cluster cycle."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
from math import lcm
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment
import motor_rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('Choose an unused output')
    root=Path(__file__).resolve().parents[1]
    source=root/'results/motor-activity-continuation/probe-destinations/case-16/periodic/cycle.npz'
    analysis=source.with_name('analysis.json')
    data=np.load(source);meta=json.loads(analysis.read_text());groups=meta['groups']
    sizes=np.array(list(map(len,groups)));period=float(data['time'][-1]);p=data['parameters']
    initial=data['full_states'][0]
    def rows(y): return np.stack((y[:8],y[8:16],np.r_[y[16:],-sum(y[16:])]),axis=1)
    a=np.array([rows(initial)[g].mean(axis=0) for g in groups])
    flow=solve_ivp(lambda t,y:motor_rows.rhs(t,y,p,8),(0,period),initial,
                   method='DOP853',rtol=5e-13,atol=5e-15,max_step=.01,dense_output=True)
    if not flow.success: raise RuntimeError(flow.message)
    records=[]
    for divisor in (2,3,4,5,6,8,10,12):
        target=np.array([rows(flow.sol(period/divisor))[g].mean(axis=0) for g in groups])
        cost=np.sum((target[:,None,:]-a[None,:,:])**2,axis=2)
        cost[sizes[:,None]!=sizes[None,:]]=1e10
        r,c=np.array(linear_sum_assignment(cost));seen=set();order=1
        for j in range(len(groups)):
            if j in seen: continue
            count=0;k=j
            while k not in seen: seen.add(k);count+=1;k=int(c[k])
            order=lcm(order,count)
        records.append(dict(time_divisor=divisor,time=period/divisor,permutation=c.tolist(),
            preserves_weights=bool(np.all(sizes==sizes[c])),permutation_order=order,
            maximum_closure_error=float(np.max(abs(target-a[c])))))
    sources=[Path(__file__),Path(motor_rows.__file__)]
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_path=source.relative_to(root).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_analysis_sha256=hashlib.sha256(analysis.read_bytes()).hexdigest(),
        source_code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        groups=groups,multiplicities=sizes.tolist(),period=period,tests=records,
        qualification='Numerical symmetry clues for future validated shooting; no exact symmetry or least-period theorem.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for p in sources: args.output.with_name(args.output.stem+'.'+p.name).write_bytes(p.read_bytes())
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print([(x['time_divisor'],x['permutation'],x['maximum_closure_error']) for x in records])


if __name__=='__main__':main()
