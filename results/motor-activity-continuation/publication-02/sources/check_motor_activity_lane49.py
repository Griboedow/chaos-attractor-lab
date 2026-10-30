"""Classify the already subthreshold activity-cohort lane with full Floquet blocks."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import motor_rows
import motor_postcritical_cycles as cycles


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output')
    args.output.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]
    source=root/'results/motor-rows/activity-path/N8-lifetime-history.npz'
    data=np.load(source);index=int(np.flatnonzero(data['source_lane_indices']==49)[0])
    state=data['final_states'][index].copy();parameters=data['parameters'][index].copy()
    groups=[[0,3],[1,6],[2,5],[4,7]]
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_path=source.relative_to(root).as_posix(),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_lane=49,initial_state=state.tolist(),
        parameters=parameters.tolist(),groups=groups,original_final_windows=data['lle_windows'][-1,index].tolist(),
        projection_error=float(np.max(abs(cycles.lift(cycles.restrict(state,groups),groups)-state))),
        source_sha256_by_name={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(motor_rows.__file__),Path(cycles.__file__)]})
    for p in [Path(__file__),Path(motor_rows.__file__),Path(cycles.__file__)]:
        (args.output/p.name).write_bytes(p.read_bytes())
    def section(t,y):return y[16]
    section.direction=1
    sol=solve_ivp(lambda t,y:motor_rows.rhs(t,y,parameters,8),(0,200),state,method='DOP853',
        rtol=2e-12,atol=2e-14,max_step=.025,events=section)
    if not sol.success:raise RuntimeError(sol.message)
    times=sol.t_events[0];points=sol.y_events[0]
    lags=[]
    for lag in range(1,min(12,len(points)-1)+1):
        lags.append(dict(lag=lag,period=float(times[-1]-times[-1-lag]),
            return_error=float(np.max(abs(points[-1]-points[-1-lag])))))
    record['lags']=lags
    candidates=[x for x in lags if x['return_error']<1e-7]
    record['status']='periodic_return_found' if candidates else 'no_close_return'
    (args.output/'sample.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(args.output/'sample.npz',parameters=parameters,time=sol.t,states=sol.y.T,
        section_time=times,section_states=points,initial_state=state,final_state=sol.y[:,-1])
    print(record['status'],lags,flush=True)
    if candidates:
        cycles.run(args.output/'sample.npz',groups,candidates[0]['lag'],args.output/'periodic')


if __name__=='__main__':main()
