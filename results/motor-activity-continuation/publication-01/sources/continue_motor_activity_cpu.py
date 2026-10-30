"""Continue the frozen full-23D motor QR trajectory without a restart.

Immutable input and per-block records preserve all states, tangents and logs.
The added velocity-square quadrature checks the weighted-volume identity.
"""
from __future__ import annotations
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import time

for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_name] = '1'
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from motor_rows import rhs, field_jac, unpack, source_force_constraint_field

ROOT = Path(__file__).resolve().parents[1]
N,D = 8,23


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, record):
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def geometry(x,p):
    c,s,positions = unpack(x,N)
    velocity = -p[0]*positions+p[1]*(s-s.mean())
    phi = (1-1/N)*(p[1]*c.sum()+p[0]*np.dot(positions,positions)/2)
    return float(phi),float(np.dot(velocity,velocity)),velocity


def endpoint_control(x,p):
    f,j = field_jac(x,p,N)
    eye = np.eye(D)
    complex_j = np.column_stack([np.imag(rhs(0,x.astype(complex)+1e-28j*v,p,N))/1e-28 for v in eye])
    c,s,positions = unpack(x,N)
    phi,speed2,velocity = geometry(x,p)
    phidot = (1-1/N)*(p[1]*f[:N].sum()+p[0]*np.dot(positions,velocity))
    constant = -2*N+(N-1)*(p[1]-p[0])
    errors = dict(jacobian_complex_step=float(np.max(abs(j-complex_j))),
        independent_force_field=float(np.max(abs(f-source_force_constraint_field(x,p,N)))),
        divergence_formula=float(abs(np.trace(j)-(-2*N-p[0]*(N-1)+p[1]*(1-1/N)*c.sum()))),
        weighted_divergence_identity=float(abs(np.trace(j)+phidot-constant+(1-1/N)*speed2)))
    return dict(errors=errors,passed=bool(max(errors.values())<1e-9),phi=phi)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/motor-activity-continuation/cpu')
    parser.add_argument('--additional',type=float,default=16000.)
    parser.add_argument('--block',type=float,default=1000.)
    parser.add_argument('--wall',type=float,default=3000.)
    parser.add_argument('--cutoff-utc',default='2026-10-30T12:50:00+00:00')
    args=parser.parse_args()
    out=args.output.resolve()
    if out.exists(): raise FileExistsError('Preserve existing continuation evidence')
    if min(args.additional,args.block,args.wall)<=0: raise ValueError('Positive durations required')
    cutoff=datetime.fromisoformat(args.cutoff_utc)
    began=time.perf_counter()
    source=ROOT/'results/motor-rows/activity-path/cpu-N8-row69-extended'
    original_path=source.with_name('cpu-N8-row69.json')
    previous=json.loads(source.with_suffix('.json').read_text())
    original=json.loads(original_path.read_text())
    with np.load(source.with_suffix('.npz')) as loaded:
        old={name:loaded[name].copy() for name in loaded.files}
    assert previous['status']=='complete' and previous['observed_tau']==4000.
    original_bytes=original_path.read_bytes()
    historical_bytes=original_bytes
    expected_original_hash=previous['continuation']['input_metadata_exact_sha256']
    if hashlib.sha256(historical_bytes).hexdigest()!=expected_original_hash:
        historical_bytes=original_bytes.replace(b'\r\n',b'\n').replace(b'\n',b'\r\n')
    assert hashlib.sha256(historical_bytes).hexdigest()==expected_original_hash
    p=old['parameters'];x=old['final_state'].copy();basis=old['final_basis'].copy()
    logs=old['log_sums'].copy();trace=float(old['integrated_trace']);flowlog=float(old['normalized_flow_log'])
    initial_elapsed=float(old['elapsed_original']);initial_observed=float(previous['observed_tau'])
    next_dt=float(old['next_qr_interval']);count=previous['qr_chunks'];rejected=previous['rejected_qr_chunks']
    for key in ('final_state','final_basis','log_sums','parameters'):
        assert np.array_equal(old[key],np.asarray(previous[key]))
    warm=float(original['settings']['warm'])
    assert warm==200.  # The extended record's warm=0 is not the original observation origin.
    warm_index=int(np.argmin(abs(old['time']-warm)))
    assert abs(float(old['time'][warm_index])-warm)<1e-8
    flow_initial=float(np.linalg.norm(rhs(0,old['states'][warm_index],p,N)))
    flow_error=abs(np.log(np.linalg.norm(rhs(0,x,p,N))/flow_initial)/initial_observed-previous['telescopic_flow_exponent'])
    assert flow_error<1e-12
    control=endpoint_control(x,p)
    assert control['passed'] and previous['short_control']['passed']
    out.mkdir(parents=True)
    shutil.copyfile(source.with_suffix('.json'),out/'input.json')
    shutil.copyfile(source.with_suffix('.npz'),out/'input.npz')
    shutil.copyfile(original_path,out/'original-cpu.json')
    (out/'original-cpu-historical-bytes.json').write_bytes(historical_bytes)
    for file in (Path(__file__),ROOT/'research/motor_rows.py',ROOT/'research/motor_rows_activity_path.py'):
        shutil.copyfile(file,out/file.name)
    input_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in
        (source.with_suffix('.json'),source.with_suffix('.npz'),original_path)}
    source_hashes={name:sha(out/name) for name in
        ('continue_motor_activity_cpu.py','motor_rows.py','motor_rows_activity_path.py')}
    phi0=geometry(x,p)[0];constant=-2*N+(N-1)*(p[1]-p[0]);theta=1-1/N
    report=dict(status='running',started_utc=datetime.now(timezone.utc).isoformat(),N=N,dimension=D,
        parameters=p.tolist(),source_hashes=source_hashes,input_hashes=input_hashes,
        initial_endpoint_control=control,inherited_short_control=previous['short_control'],
        settings=dict(additional_requested_tau=args.additional,target_observed_tau=initial_observed+args.additional,
            block_tau=args.block,rtol=1e-11,atol=1e-13,max_step=.03,qr_max=.2,
            wall_budget_seconds=args.wall,cutoff_utc=cutoff.isoformat(),warmup=0,tangent_reset=False),
        restoration=dict(state_basis_logs_match=True,original_observed_tau=initial_observed,
            original_elapsed_tau=initial_elapsed,original_warmup_tau=warm,
            original_flow_sample_index=warm_index,original_flow_sample_tau=float(old['time'][warm_index]),
            original_flow_norm=flow_initial,recovered_flow_rate_error=float(flow_error),
            original_metadata_current_sha256=sha(original_path),
            original_metadata_historical_sha256=expected_original_hash,
            original_metadata_CRLF_reconstructed=historical_bytes!=original_bytes,
            pitfall_avoided='The extended settings warm=0 describes its added segment, not the original observation start.'),
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,blas_threads=1),
        windows=copy.deepcopy(previous['windows']),new_blocks=[],
        max_local_trace_error=previous['max_local_trace_error'],
        max_local_flow_identity_error=previous['max_local_flow_identity_error'],
        max_local_condition=previous['max_local_condition'],max_new_weighted_volume_error=0.,
        qualification='Continuation of one frozen numerical trajectory in original tau; no lane reselection, warmup, tangent reset, asymptotic-chaos theorem or priority claim.')
    times=old['time'].tolist();states=list(old['states']);section_times=old['section_time'].tolist()
    sections=list(old['section_states']);checkpoints=list(old['checkpoint'])
    added=0.;duration_correction=0.;new_trace=0.;speed_integral=0.;new_flow=0.
    new_logs=np.zeros(D);block_index=0

    def state_record():
        observed=initial_observed+added;phi=geometry(x,p)[0]
        report.update(wall_seconds=time.perf_counter()-began,added_tau=added,observed_tau=observed,
            elapsed_original=initial_elapsed+added,qr_chunks=count,rejected_qr_chunks=rejected,
            final_state=x.tolist(),final_basis=basis.tolist(),log_sums=logs.tolist(),integrated_trace=trace,
            normalized_flow_log=flowlog,next_qr_interval=next_dt,spectrum=(logs/observed).tolist(),
            mean_divergence=trace/observed,sum_trace_error=abs(logs.sum()-trace)/observed,
            flow_exponent=flowlog/observed,
            telescopic_flow_exponent=float(np.log(np.linalg.norm(rhs(0,x,p,N))/flow_initial)/observed),
            upward_X1_returns=len(sections),new_segment=dict(duration=added,integrated_trace=new_trace,
                integrated_speed_squared=speed_integral,initial_phi=phi0,final_phi=phi,
                volume_constant=constant,identity_residual=new_trace+phi-phi0-constant*added+theta*speed_integral,
                qr_trace_residual=float(new_logs.sum()-new_trace),normalized_flow_log=new_flow))
        return report

    def final_archive():
        arrays=dict(parameters=p,time=np.asarray(times),states=np.asarray(states),section_time=np.asarray(section_times),
            section_states=np.asarray(sections),checkpoint=np.asarray(checkpoints),final_state=x,final_basis=basis,
            log_sums=logs,integrated_trace=trace,normalized_flow_log=flowlog,elapsed_original=initial_elapsed+added,
            next_qr_interval=next_dt,added_tau=added,integrated_speed_squared_new=speed_integral,
            integrated_trace_new=new_trace,log_sums_new=new_logs)
        preservation={name:arrays[name][:len(old[name])].tobytes()==old[name].tobytes()
            for name in ('time','states','section_time','section_states','checkpoint')}
        report['inherited_history_prefix_bitwise_checks']=preservation
        assert all(preservation.values())
        np.savez_compressed(out/'continuation.npz',**arrays)
        report['continuation_npz_sha256']=sha(out/'continuation.npz')
        report['final_endpoint_control']=endpoint_control(x,p)
        report['source_inputs_unchanged']=all(sha(ROOT/name)==digest for name,digest in input_hashes.items())
        assert report['source_inputs_unchanged'] and report['final_endpoint_control']['passed']
        sample=np.asarray(states)[warm_index+1:];positions=np.c_[sample[:,16:],-sample[:,16:].sum(axis=1)]
        pairs=[]
        for i in range(N):
            for j in range(i+1,N):
                distances=np.sqrt((sample[:,i]-sample[:,j])**2+(sample[:,N+i]-sample[:,N+j])**2+(positions[:,i]-positions[:,j])**2)
                pairs.append(dict(rows=[i,j],minimum=float(distances.min()),maximum=float(distances.max()),
                                  rms=float(np.sqrt(np.mean(distances**2)))))
        report['sampled_pair_distances_zero_based']=pairs
        report['sampled_domain']=dict(max_disk_squared=float(np.max(sample[:,:N]**2+sample[:,N:2*N]**2)),
            max_abs_sum_X=float(np.max(abs(positions.sum(axis=1)))))
        write_json(out/'analysis.json',state_record())

    def augmented(t,z):
        ff,jj=field_jac(z[:D],p,N)
        return np.r_[ff,(jj@z[D:D+D*D].reshape(D,D)).ravel(),np.trace(jj),geometry(z[:D],p)[1]]

    def section(t,z): return z[2*N]
    section.direction=1
    write_json(out/'progress.json',state_record())
    print(json.dumps(dict(status='restored',observed=initial_observed,elapsed=initial_elapsed,
        flow_error=float(flow_error),requested_additional=args.additional)),flush=True)
    try:
        while added<args.additional-1e-8:
            block_index+=1;length=min(args.block,args.additional-added);spent=0.;spent_correction=0.
            wlogs=np.zeros(D);wtrace=wspeed=wflow=0.;wphi=geometry(x,p)[0]
            flow_start=float(np.linalg.norm(rhs(0,x,p,N)))
            first_state=len(states);first_section=len(sections)
            while spent<length-1e-11:
                if time.perf_counter()-began>=args.wall or datetime.now(timezone.utc)>=cutoff:
                    report['status']='wall_or_absolute_cutoff';break
                dt=min(next_dt,length-spent);f0=rhs(0,x,p,N);nf=float(np.linalg.norm(f0))
                if nf<1e-14: report['status']='near_equilibrium';break
                phi_start=geometry(x,p)[0]
                sol=solve_ivp(augmented,(0,dt),np.r_[x,basis.ravel(),0.,0.],method='DOP853',
                    rtol=1e-11,atol=1e-13,max_step=min(.03,dt),events=section)
                if not sol.success: raise ArithmeticError(sol.message)
                end=sol.y[:D,-1];matrix=sol.y[D:D+D*D,-1].reshape(D,D)
                new_basis,triangular=np.linalg.qr(matrix);diag=np.diag(triangular)
                growth=np.log(abs(diag));chunk_trace=float(sol.y[-2,-1]);chunk_speed=float(sol.y[-1,-1])
                error=float(abs(growth.sum()-chunk_trace));condition=float(np.linalg.cond(matrix))
                weighted_error=abs(chunk_trace+geometry(end,p)[0]-phi_start-constant*dt+theta*chunk_speed)
                if np.any(diag==0) or error>1e-8 or condition>1e10 or weighted_error>1e-8:
                    next_dt=dt/2;rejected+=1
                    if next_dt<1e-6: raise ArithmeticError('QR/volume control failed below minimum interval')
                    continue
                propagated=matrix@(basis.T@(f0/nf));f1=rhs(0,end,p,N)
                flow_error_local=float(np.linalg.norm(propagated*nf-f1)/(1+np.linalg.norm(f1)))
                fg=float(np.log(np.linalg.norm(propagated)))
                report['max_local_trace_error']=max(report['max_local_trace_error'],error)
                report['max_local_flow_identity_error']=max(report['max_local_flow_identity_error'],flow_error_local)
                report['max_local_condition']=max(report['max_local_condition'],condition)
                report['max_new_weighted_volume_error']=max(report['max_new_weighted_volume_error'],weighted_error)
                x=end;basis=new_basis;count+=1;next_dt=min(.2,next_dt*1.05)
                adjusted=dt-duration_correction;total=added+adjusted;duration_correction=(total-added)-adjusted;added=total
                adjusted=dt-spent_correction;total=spent+adjusted;spent_correction=(total-spent)-adjusted;spent=total
                logs+=growth;trace+=chunk_trace;flowlog+=fg;new_logs+=growth;new_trace+=chunk_trace
                speed_integral+=chunk_speed;new_flow+=fg;wlogs+=growth;wtrace+=chunk_trace;wspeed+=chunk_speed;wflow+=fg
                elapsed=initial_elapsed+added
                for hit,z in zip(sol.t_events[0],sol.y_events[0]):
                    if hit>1e-10: section_times.append(elapsed-dt+float(hit));sections.append(z[:D])
                times.append(elapsed);states.append(x.copy())
            if spent>0:
                block=dict(index=block_index,end_tau=initial_elapsed+added,end_observed_tau=initial_observed+added,
                    duration=spent,complete=abs(spent-length)<1e-8,spectrum=(wlogs/spent).tolist(),
                    cumulative_spectrum=(logs/(initial_observed+added)).tolist(),mean_divergence=wtrace/spent,
                    sum_trace_error=abs(wlogs.sum()-wtrace)/spent,flow_exponent=wflow/spent,
                    telescopic_flow_exponent=float(np.log(np.linalg.norm(rhs(0,x,p,N))/flow_start))/spent,
                    upward_X1_returns=len(sections)-first_section,integrated_speed_squared=wspeed,
                    weighted_volume_residual=wtrace+geometry(x,p)[0]-wphi-constant*spent+theta*wspeed,
                    wall_seconds=time.perf_counter()-began)
                filename=f'block-{block_index:02d}'
                np.savez_compressed(out/(filename+'.npz'),time=np.asarray(times[first_state:]),
                    states=np.asarray(states[first_state:]),section_time=np.asarray(section_times[first_section:]),
                    section_states=np.asarray(sections[first_section:]),final_state=x,final_basis=basis,
                    log_sums=logs,integrated_trace=trace,normalized_flow_log=flowlog,next_qr_interval=next_dt,
                    added_tau=added,integrated_speed_squared_new=speed_integral,integrated_trace_new=new_trace,
                    block_log_sums=wlogs,block_integrated_trace=wtrace,block_integrated_speed_squared=wspeed)
                block['npz_sha256']=sha(out/(filename+'.npz'));write_json(out/(filename+'.json'),block)
                report['windows'].append(block);report['new_blocks'].append(block)
                checkpoints.append(np.r_[initial_elapsed+added,x,basis.ravel(),logs,trace,flowlog,next_dt])
                write_json(out/'progress.json',state_record())
                print(json.dumps(dict(block=block_index,observed=initial_observed+added,
                    leading_spectrum=report['spectrum'][:8],weighted_residual=block['weighted_volume_residual'],
                    wall=round(time.perf_counter()-began,1))),flush=True)
            if report['status']!='running': break
        if report['status']=='running': report['status']='complete'
    except Exception as error:
        report.update(status='error',error=repr(error))
        raise
    finally:
        final_archive()
    if report['status']!='complete': raise SystemExit(2)


if __name__=='__main__': main()
