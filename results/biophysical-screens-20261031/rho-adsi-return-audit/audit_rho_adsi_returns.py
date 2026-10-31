"""Audit numerical periodic returns in completed homogeneous ADSI histories."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from rho_adsi import independent_field


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    records=[]
    for path in sorted(args.source.glob('case-*.npz')):
        data=np.load(path);p=data['parameters'];events=data['section_states'];times=data['section_times']
        assert len(events)>101
        periods=np.diff(times[-101:]);period=float(np.median(periods));initial=events[-1]
        endpoints={}
        for method in ['Radau','DOP853']:
            sol=solve_ivp(lambda t,x:independent_field(x,p),(0,period),initial,method=method,
                          rtol=2e-12,atol=2e-14,max_step=.005)
            assert sol.success
            endpoints[method]=sol.y[:,-1]
        records.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            section_events=len(events),last_return_count=100,median_period_tau=period,
            period_min_tau=float(periods.min()),period_max_tau=float(periods.max()),
            maximum_consecutive_section_distance=float(np.max(np.linalg.norm(np.diff(events[-101:],axis=0),axis=1))),
            return_residuals={m:float(np.linalg.norm(x-initial)) for m,x in endpoints.items()},
            independent_solver_discrepancy=float(np.linalg.norm(endpoints['Radau']-endpoints['DOP853'])),
            interpretation='Numerical periodic-return evidence; not a rigorous periodic-orbit or least-period certificate.'))
    report=dict(source=str(args.source.resolve()),time='tau=.025*t, physical t in seconds',
        method='Last100 upward section returns; independent one-period Radau and DOP853 at rtol2e-12,atol2e-14,maxstep.005tau.',
        runs=records)
    (out/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf8')
    for name in ['audit_rho_adsi_returns.py','rho_adsi.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
