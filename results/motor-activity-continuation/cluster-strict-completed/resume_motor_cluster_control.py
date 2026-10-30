"""Resume all conditional QR accumulators from a verified complete block."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
import motor_rows
import motor_postcritical_cycles as cycles
import check_motor_activity_clusters as parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--observed-total',type=float,default=10000.)
    parser.add_argument('--wall',type=float,default=900.)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output')
    frozen=json.loads((args.input/'checkpoint.json').read_text())
    if frozen['status']!='passed' or not all(frozen['checks'].values()):
        raise ArithmeticError('Verified complete checkpoint required')
    for name,expected in frozen['hashes'].items():
        if sha(args.input/name)!=expected:
            raise ArithmeticError('Checkpoint hash mismatch: '+name)
    old=json.loads((args.input/'analysis.json').read_text());data=np.load(args.input/'history.npz')
    observed=float(old['observed_tau']);elapsed=float(old['elapsed_tau'])
    if elapsed!=200.+observed or args.observed_total<=observed or (args.observed_total-observed)%500:
        raise ValueError('A complete500-block start and later500-block target are required')
    p=np.asarray(old['parameters']);w=np.asarray(old['weights']);groups=old['groups'];d=old['quotient_dimension']
    paired=[j for j,g in enumerate(groups) if len(g)>1]
    state=data['final_state'].copy();tangent=data['final_quotient_tangent'].copy()
    bases=data['final_contrast_bases'].copy();logs=data['log_sums'].copy()
    times=list(data['time']);states=list(data['states'].copy());windows=list(old['windows'])
    settings=old['settings'];step=settings['max_step'];rtol=settings['rtol'];atol=settings['atol']
    args.output.mkdir(parents=True)
    for name in ['analysis.json','history.npz','checkpoint.json']:
        (args.output/('input-'+name)).write_bytes((args.input/name).read_bytes())
    sources=[Path(__file__),Path(parent.__file__),Path(cycles.__file__),Path(motor_rows.__file__)]
    for source in sources:
        (args.output/source.name).write_bytes(source.read_bytes())
    record=dict(status='running',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        input_folder=str(args.input),input_hashes={name:sha(args.input/name) for name in ['analysis.json','history.npz','checkpoint.json']},
        source_hashes={q.name:sha(q) for q in sources},source_lane=old['source_lane'],groups=groups,weights=w.tolist(),parameters=p.tolist(),
        settings=dict(max_step=step,rtol=rtol,atol=atol,qr_interval=1.,new_warmup=0.,observed_total_requested=args.observed_total),
        exact_state_tangent_basis_logs_restored=True,inherited_observed_tau=observed,windows=windows,
        caveat='Finite-time conditional QR continuation; a complete-block checkpoint was used rather than a partial wall-budget endpoint.')
    def rhs(t,z):
        x=z[:d];f,j=cycles.field_jac(x,p,w)
        parts=[f,j@z[d:2*d]]
        for k,a in enumerate(paired):
            parts.append((parent.contrast_matrix(x,p,w,a)@z[2*d+9*k:2*d+9*k+9].reshape(3,3)).ravel())
        return np.concatenate(parts)
    initial=np.r_[state,tangent,bases.ravel()]
    controls=[]
    for h in [step,step/2]:
        result=solve_ivp(rhs,(0.,1.),initial,method='DOP853',rtol=rtol,atol=atol,max_step=h)
        if not result.success:
            raise RuntimeError(result.message)
        controls.append(result.y[:,-1])
    error=float(np.max(abs(controls[0]-controls[1])))
    record['short_step_control_maximum_state_and_tangent_error']=error
    if error>1e-10:
        raise ArithmeticError('Resume step control failed')
    started=time.perf_counter();block_logs=np.zeros_like(logs);block_duration=0.
    def save():
        record.update(observed_tau=observed,elapsed_tau=elapsed,wall_seconds=time.perf_counter()-started,
            quotient_leading_ftle=float(logs[0]/observed),contrast_ftle=(logs[1:].reshape(-1,3)/observed).tolist(),
            partial_block_duration=block_duration)
        np.savez_compressed(args.output/'history.npz',time=np.asarray(times),states=np.asarray(states),parameters=p,
            final_state=state,final_quotient_tangent=tangent,final_contrast_bases=bases,log_sums=logs)
        record['history_sha256']=sha(args.output/'history.npz')
        (args.output/'analysis.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    save()
    while observed<args.observed_total:
        if time.perf_counter()-started>args.wall:
            record['status']='wall budget incomplete';save();raise SystemExit(2)
        z=np.r_[state,tangent,bases.ravel()]
        sol=solve_ivp(rhs,(0.,1.),z,method='DOP853',rtol=rtol,atol=atol,max_step=step)
        if not sol.success:
            raise RuntimeError(sol.message)
        end=sol.y[:,-1];state=end[:d];tangent=end[d:2*d]
        norm=np.linalg.norm(tangent);tangent/=norm;increment=np.zeros_like(logs);increment[0]=np.log(norm)
        for k,a in enumerate(paired):
            matrix=end[2*d+9*k:2*d+9*k+9].reshape(3,3);q,r=np.linalg.qr(matrix)
            bases[k]=q;increment[1+3*k:4+3*k]=np.log(abs(np.diag(r)))
        logs+=increment;block_logs+=increment;observed+=1.;elapsed+=1.;block_duration+=1.
        times.append(elapsed);states.append(state.copy())
        if block_duration==500.:
            windows.append(dict(duration=500.,end_tau=elapsed,quotient_leading_ftle=float(block_logs[0]/500.),
                contrast_ftle=(block_logs[1:].reshape(-1,3)/500.).tolist()))
            block_logs=np.zeros_like(logs);block_duration=0.;save()
            print('resumed',observed,record['quotient_leading_ftle'],[v[0] for v in record['contrast_ftle']],flush=True)
    record['status']='complete';save()


if __name__=='__main__':
    main()
