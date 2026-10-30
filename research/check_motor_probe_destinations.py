"""Full-space continuations and complete Floquet checks of four declared probes."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import datetime
import hashlib
import json
from pathlib import Path
import time
import traceback
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment, minimize_scalar
import motor_rows
import motor_postcritical_cycles as cycles

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/motor-activity-continuation/full-space-perturbations/outcomes.npz'
REFERENCE=ROOT/'results/motor-activity-continuation/lane49/periodic/cycle.npz'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(state):
    return np.stack((state[:8],state[8:16],np.r_[state[16:],-sum(state[16:])]),axis=1)


def groups_at(state):
    row=rows(state);todo=set(range(8));groups=[]
    while todo:
        first=min(todo)
        group=[j for j in sorted(todo) if np.linalg.norm(row[j]-row[first])<1e-10]
        groups.append(group);todo-=set(group)
    return groups


def reference_comparison(target_cycle):
    reference=np.load(REFERENCE)
    p=reference['parameters']
    target=rows(target_cycle['full_states'][0])
    period=float(reference['time'][-1])
    flow=solve_ivp(lambda t,y:motor_rows.rhs(t,y,p,8),(0,period),reference['full_states'][0],
                   method='DOP853',rtol=2e-12,atol=2e-14,max_step=.015,dense_output=True)
    if not flow.success: raise RuntimeError(flow.message)
    def distance(t,details=False):
        ref=rows(flow.sol(float(t)%period))
        cost=np.sum((target[:,None,:]-ref[None,:,:])**2,axis=2)
        a,b=linear_sum_assignment(cost)
        value=float(np.mean(cost[a,b]))
        return (value,b.tolist()) if details else value
    grid=np.linspace(0,period,513)[:-1]
    center=float(grid[np.argmin([distance(t) for t in grid])])
    refined=minimize_scalar(distance,bounds=(center-period/256,center+period/256),method='bounded',
                            options=dict(xatol=1e-13,maxiter=200))
    value,permutation=distance(refined.x,True)
    return dict(reference_path=REFERENCE.relative_to(ROOT).as_posix(),reference_sha256=sha(REFERENCE),
        reference_period=period,target_period=float(target_cycle['time'][-1]),
        phase=float(refined.x%period),target_row_to_reference_row=permutation,
        full_row_rms_distance=float(np.sqrt(value)),optimizer_success=bool(refined.success),
        comparison_scope='Floating-point phase optimization and optimal row assignment; not a validated orbit equivalence theorem.')


def run_case(case,output,duration):
    out=Path(output)/f'case-{case}';out.mkdir(parents=True)
    began=time.perf_counter()
    old=np.load(SOURCE);metadata=json.loads(SOURCE.with_suffix('.json').read_text())
    initial=old['final_states'][case].copy();p=old['parameters'][case].copy()
    record=dict(case=case,source_lane=int(old['source_lane_indices'][case]),parameters=p.tolist(),
        source_npz_sha256=sha(SOURCE),source_json_sha256=sha(SOURCE.with_suffix('.json')),
        source_case=metadata['protocol']['cases'][case],source_lle_windows=old['lle_windows'][case].tolist(),
        sample_scope='Unprojected full 23D state continuation from the archived probe endpoint.',
        settings=dict(duration=duration,method='DOP853',rtol=2e-12,atol=2e-14,max_step=.025,
                      section='X0=0, upward',group_tolerance=1e-10,close_return_lags=list(range(1,17)),
                      shooting_trigger='Smallest lag with final physical max error <1e-5 and recent median max error <1e-4.'),
        status='sampling',limitations=['Floating-point conditional destination diagnosis, not an interval proof or basin measure.',
            'Any cluster projection is measured and declared only after full-space sampling.',
            'An apparently periodic endpoint is not assumed stable; the complete 23D monodromy is requested.'])
    def save():
        record['wall_seconds']=time.perf_counter()-began
        (out/'analysis.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    save()
    try:
        def event(t,y): return y[16]
        event.direction=1.;event.terminal=False
        sol=solve_ivp(lambda t,y:motor_rows.rhs(t,y,p,8),(0,duration),initial,
            method='DOP853',rtol=2e-12,atol=2e-14,max_step=.025,
            t_eval=np.linspace(0,duration,int(duration*20)+1),events=event)
        if not sol.success: raise RuntimeError(sol.message)
        section_time=sol.t_events[0];section_states=sol.y_events[0]
        final=sol.y[:,-1];groups=groups_at(final)
        sample=out/'sample.npz'
        np.savez_compressed(sample,time=sol.t,full_states=sol.y.T,parameters=p,
            section_time=section_time,section_states=section_states,initial_state=initial,final_state=final)
        returns=[]
        for lag in range(1,min(16,len(section_states)-1)+1):
            errors=[float(np.max(abs(rows(section_states[-1-k])-rows(section_states[-1-k-lag]))))
                    for k in range(min(8,len(section_states)-lag))]
            returns.append(dict(lag=lag,period=float(section_time[-1]-section_time[-1-lag]),
                final_physical_max_error=errors[0],recent_median_max_error=float(np.median(errors)),
                recent_max_error=max(errors),recent_errors=errors,
                triggers_shooting=bool(errors[0]<1e-5 and np.median(errors)<1e-4)))
        record.update(status='sampled',section_count=len(section_time),groups=groups,returns=returns,
            initial_groups=groups_at(initial),sample_sha256=sha(sample),
            endpoint_projection_error=float(np.max(abs(cycles.lift(cycles.restrict(final,groups),groups)-final))))
        save()
        eligible=[r for r in returns if r['triggers_shooting']]
        if not eligible:
            record['status']='sampled_no_close_return';save();return record
        chosen=eligible[0];record.update(status='shooting',selected_lag=chosen['lag']);save()
        with (out/'shooting.log').open('w',encoding='utf-8') as log, contextlib.redirect_stdout(log):
            cycles.run(sample,groups,chosen['lag'],out/'periodic')
        periodic=json.loads((out/'periodic/analysis.json').read_text())
        cycle=np.load(out/'periodic/cycle.npz')
        record.update(status='periodic_destination_numerically_resolved',periodic_analysis_sha256=sha(out/'periodic/analysis.json'),
            period=periodic['period'],full_unstable_count=periodic['full_unstable_count_at_1e_minus_6'],
            full_neutral_count=periodic['full_neutral_count_at_1e_minus_6'],
            full_multipliers=periodic['full_multipliers'],closure_error=periodic['full_closure_error'],
            reference_comparison=reference_comparison(cycle))
        save()
    except Exception as error:
        record.update(status='failed',error=repr(error),traceback=traceback.format_exc());save()
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--duration',type=float,default=500.)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('Preserve evidence: choose a new output directory')
    if not 200<=args.duration<=500 or args.workers<1: raise ValueError('Require 200..500 tau and positive workers')
    args.output.mkdir(parents=True)
    sources=[Path(__file__),Path(motor_rows.__file__),Path(cycles.__file__)]
    protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases=[16,17,18,19],
        source=SOURCE.relative_to(ROOT).as_posix(),source_sha256=sha(SOURCE),duration=args.duration,
        source_snapshot_sha256={p.name:sha(p) for p in sources},
        selection='All four physical perturbations of original lane49, selected for follow-up after the near-zero final GPU window. No cases are removed by the outcome.',
        immutable_original_data=True,scope='New conditional numerical follow-up, separate from the frozen activity cohort and local interval theorems.')
    for source in sources:(args.output/source.name).write_bytes(source.read_bytes())
    (args.output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run_case,case,args.output.resolve(),args.duration):case for case in protocol['cases']}
        for future in as_completed(futures):
            case=futures[future]
            try: record=future.result()
            except Exception as error: record=dict(case=case,status='worker_failed',error=repr(error))
            results.append(record)
            (args.output/'summary.json').write_text(json.dumps(dict(protocol=protocol,cases=sorted(results,key=lambda r:r['case'])),indent=2)+'\n',encoding='utf-8')
            print('case',case,record['status'],'period',record.get('period'),'unstable',record.get('full_unstable_count'),flush=True)


if __name__=='__main__':main()
