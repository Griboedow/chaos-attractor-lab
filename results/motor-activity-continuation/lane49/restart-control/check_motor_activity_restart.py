"""Separate tangent loss from physical departures in the original lane49."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
import motor_rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output')
    args.output.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]
    source=root/'results/motor-rows/activity-path/N8-lifetime-history.npz'
    old=np.load(source);i=int(np.flatnonzero(old['source_lane_indices']==49)[0])
    x0=old['final_states'][i].copy();p=old['parameters'][i].copy()
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_path=source.relative_to(root).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        parameters=p.tolist(),source_lane=49,initial_state=x0.tolist(),
        tangent_cases=[],physical_departures=[],status='running',
        settings=dict(method='DOP853',rtol=2e-11,atol=2e-13,max_step=.03),
        limitations=['Floating-point diagnostic, not an interval certificate.',
            'A fresh tangent changes only the diagnostic. Physical-departure cases deliberately change the stated initial state.',
            'Departures from a saddle do not prove chaos or a final attracting regime.'])
    sources=[Path(__file__),Path(motor_rows.__file__)]
    record['source_snapshot_sha256']={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in sources}
    for q in sources:(args.output/q.name).write_bytes(q.read_bytes())
    start=time.perf_counter()
    pairs=[(0,3),(1,6),(2,5),(4,7)]

    def contrasts(z):
        xx=np.r_[z[16:],-sum(z[16:])]
        return [float(np.linalg.norm([z[a]-z[b],z[8+a]-z[8+b],xx[a]-xx[b]])) for a,b in pairs]

    def save():
        record['wall_seconds']=time.perf_counter()-start
        (args.output/'analysis.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

    for label,tangent in [('carried',old['final_tangents'][i].copy()),('fresh_sine',np.sin(np.arange(23)+1.))]:
        tangent=tangent/np.linalg.norm(tangent);state=x0.copy();logs=0.;windows=[]
        initial=tangent.copy();times=[];states=[]
        def rhs(t,z):
            f,j=motor_rows.field_jac(z[:23],p,8)
            return np.r_[f,j@z[23:]]
        for block in range(5):
            block_log=0.
            for chunk in range(50):
                sol=solve_ivp(rhs,(0,2.),np.r_[state,tangent],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.03)
                if not sol.success:raise RuntimeError(sol.message)
                state=sol.y[:23,-1];v=sol.y[23:,-1];norm=np.linalg.norm(v)
                block_log+=np.log(norm);tangent=v/norm
                times.append(block*100.+(chunk+1)*2.);states.append(state.copy())
            logs+=block_log;windows.append(float(block_log/100.))
        case=dict(label=label,observed_tau=500.,initial_tangent=initial.tolist(),
            initial_tangent_pair_norms=contrasts(initial),initial_state_pair_norms=contrasts(x0),
            window_ftle=windows,mean_ftle=float(logs/500.),final_state_pair_norms=contrasts(state))
        record['tangent_cases'].append(case)
        np.savez_compressed(args.output/(label+'.npz'),time=np.array(times),states=np.array(states),
            initial_state=x0,initial_tangent=initial,final_state=state,final_tangent=tangent)
        save();print('tangent',label,case['mean_ftle'],windows,flush=True)
    for epsilon in (1e-7,1e-9):
        state=x0.copy();state[0]+=epsilon;state[3]-=epsilon
        def threshold(t,z):return np.linalg.norm([z[0]-z[3],z[8]-z[11],z[16]-z[19]])-1e-3
        threshold.direction=1
        times=np.arange(0.,1000.01,.1)
        sol=solve_ivp(lambda t,z:motor_rows.rhs(t,z,p,8),(0,1000.),state,method='DOP853',
            rtol=2e-11,atol=2e-13,max_step=.03,t_eval=times,events=threshold)
        if not sol.success:raise RuntimeError(sol.message)
        name='physical-1e-7' if epsilon==1e-7 else 'physical-1e-9'
        np.savez_compressed(args.output/(name+'.npz'),time=sol.t,states=sol.y.T,
            initial_state=state,threshold_times=sol.t_events[0])
        case=dict(label=name,half_c_contrast=epsilon,perturbed_rows=[0,3],
            initial_delta_from_archived_state=(state-x0).tolist(),initial_pair_norms=contrasts(state),
            departure_threshold_full_pair_norm=1e-3,threshold_times=sol.t_events[0].tolist(),
            final_pair_norms=contrasts(sol.y[:,-1]),observed_tau=1000.)
        record['physical_departures'].append(case);save()
        print('physical',epsilon,'firstdeparture',case['threshold_times'][:1],flush=True)
    record['status']='complete';save()


if __name__=='__main__':main()
