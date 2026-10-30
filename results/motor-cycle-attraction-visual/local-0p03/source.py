"""Integrate declared full-space starts near the certified five-cluster cycle.

Numerical illustration of local attraction; not a validated basin enclosure.
Every declared start and both solver settings are retained.
"""
import argparse,datetime,hashlib,json,math
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree
import motor_rows

ROOT=Path(__file__).resolve().parents[1]
GROUPS=[[0,3],[1,2],[4],[5,6],[7]]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,data):Path(p).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
def physical(z):
    z=np.asarray(z)
    return np.concatenate((z,-np.sum(z[...,16:],axis=-1,keepdims=True)),axis=-1)
def midpoint(v):return math.ldexp(float(v['midpoint'][0]),v['midpoint'][1])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/motor-cycle-attraction-visual')
    parser.add_argument('--amplitude',type=float,default=.001)
    args=parser.parse_args();folder=args.output.resolve()
    folder.mkdir(parents=True,exist_ok=False)
    certificate_path=ROOT/'results/motor-five-cycle-proof/recheck-01.json'
    certificate=json.loads(certificate_path.read_text())
    assert certificate['status']=='attracting_five_cluster_cycle_certified'
    period=midpoint(certificate['positive_period'])
    q=np.array([midpoint(v) for v in certificate['initial_state_box']])
    weights=np.array([len(g)/8 for g in GROUPS])
    x=np.r_[q[10:],-np.dot(weights[:4],q[10:])/weights[4]]
    full=np.zeros((8,3))
    for i,g in enumerate(GROUPS):full[g]=[q[i],q[5+i],x[i]]
    initial=np.r_[full[:,0],full[:,1],full[:7,2]]
    with localcontext() as context:
        context.prec=80
        parameters=np.array([float(Decimal(129)/55),float(1+(Decimal(19)/20)/(Decimal(32001)/500000).sqrt()),1.])
    field=lambda t,z:motor_rows.rhs(t,z,parameters,8)
    reference=solve_ivp(field,(0.,period),initial,method='DOP853',rtol=3e-13,atol=3e-15,
                        max_step=.0125,dense_output=True)
    assert reference.success
    closure=float(np.max(abs(reference.y[:,-1]-initial)))
    assert closure<1e-9
    phases=np.arange(6)/6
    rng=np.random.default_rng(2026103031)
    directions=rng.normal(size=(6,24))
    directions[:,16:]-=directions[:,16:].mean(axis=1,keepdims=True)
    directions/=np.linalg.norm(directions,axis=1,keepdims=True)
    amplitude=args.amplitude
    if not 0<amplitude<10:raise ValueError('Amplitude must be positive and bounded')
    starts=np.stack([reference.sol(p*period)+amplitude*d[:-1] for p,d in zip(phases,directions)])
    times=np.linspace(0.,6*period,3001)
    protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scope='Six predeclared generic full23D starts; no cluster restriction or projection during evolution.',
        parameters=parameters.tolist(),certificate=certificate_path.relative_to(ROOT).as_posix(),certificate_sha256=sha(certificate_path),
        seed=2026103031,phases_as_period_fractions=phases.tolist(),full24D_perturbation_norm=amplitude,
        initial_states=starts.tolist(),full24D_zero_sum_X_directions=directions.tolist(),horizon_periods=6,
        period=period,reference_closure_error=closure,source_sha256=sha(__file__),field_source_sha256=sha(motor_rows.__file__),
        settings=[dict(name='main',rtol=2e-11,atol=2e-13,max_step=.04),
                  dict(name='strict',rtol=2e-12,atol=2e-14,max_step=.02)],
        qualification='Floating conversions and ordinary numerical integration; no basin radius or global attraction proof.')
    write(folder/'protocol.json',protocol)
    (folder/'source.py').write_bytes(Path(__file__).read_bytes())
    (folder/'motor_rows.py').write_bytes(Path(motor_rows.__file__).read_bytes())
    phase_grid=np.linspace(0.,period,12001,endpoint=False)
    cycle=reference.sol(phase_grid).T
    tree=cKDTree(physical(cycle))
    spacing=period/len(phase_grid)
    def distance_to_orbit(states):
        target=physical(states)
        nearest=tree.query(target)[1]
        distances=[];matched=[]
        for z,k in zip(target,nearest):
            center=phase_grid[k]
            def objective(t):return float(np.sum((physical(reference.sol(t%period))-z)**2))
            fit=minimize_scalar(objective,bounds=(center-2*spacing,center+2*spacing),method='bounded',
                                options=dict(xatol=1e-12,maxiter=80))
            assert fit.success
            distances.append(math.sqrt(max(fit.fun,0.)));matched.append(float(fit.x%period))
        return np.array(distances),np.array(matched)
    histories=[];records=[]
    sample_indices=np.arange(0,len(times),10)
    for i,start in enumerate(starts):
        solutions=[]
        for setting in protocol['settings']:
            sol=solve_ivp(field,(0.,times[-1]),start,method='DOP853',t_eval=times,
                          **{k:setting[k] for k in ['rtol','atol','max_step']})
            if not sol.success:
                write(folder/f'failure-{i:02}-{setting["name"]}.json',dict(message=sol.message))
                raise RuntimeError(sol.message)
            solutions.append(sol.y.T)
        distance,phase=distance_to_orbit(solutions[1][sample_indices])
        error=float(np.max(abs(solutions[0]-solutions[1])))
        pair_separations=[float(np.linalg.norm(physical(start).reshape(3,8).T[g[0]]-physical(start).reshape(3,8).T[g[1]]))
                          for g in GROUPS if len(g)==2]
        record=dict(index=i,solver_max_state_difference=error,initial_distance_to_orbit=float(distance[0]),
                    final_distance_to_orbit=float(distance[-1]),minimum_distance_to_orbit=float(distance.min()),
                    maximum_distance_to_orbit=float(distance.max()),initial_pair_separations=pair_separations,
                    all_pairs_initially_split=all(v>1e-8 for v in pair_separations),
                    integration_control_passed=error<1e-7,
                    approached_reference_cycle=bool(distance[-1]<1e-5 and distance[-1]<distance[0]*1e-3))
        np.savez_compressed(folder/f'trajectory-{i:02}.npz',time=times,main_states=solutions[0],strict_states=solutions[1],
                            distance_time=times[sample_indices],full_state_orbit_distance=distance,matched_phase=phase,
                            initial_state=start,parameters=parameters)
        record['trajectory_sha256']=sha(folder/f'trajectory-{i:02}.npz')
        records.append(record);histories.append(solutions[1])
        write(folder/'progress.json',dict(completed=len(records),declared=len(starts),runs=records))
        print(json.dumps(record),flush=True)
    np.savez_compressed(folder/'reference.npz',time=phase_grid,states=cycle,parameters=parameters,period=period)
    result=dict(status='complete',runs=records,all_declared_starts_retained=True,
                all_integration_controls_passed=all(v['integration_control_passed'] for v in records),
                all_six_approach_same_reference_orbit=all(v['approached_reference_cycle'] for v in records),
                distance_definition='Euclidean norm of all24 physical(c,s,X) components, minimizing phase on the saved periodic orbit.',
                phase_matching='12,001-point full-state nearest-neighbor search followed by bounded dense-orbit phase minimization.',
                reference_sha256=sha(folder/'reference.npz'),protocol_sha256=sha(folder/'protocol.json'),
                limitation='Numerical illustration for six specified starts. Neither a validated basin boundary nor a same-phase synchronization claim.')
    write(folder/'analysis.json',result)
    if not result['all_integration_controls_passed']:raise ArithmeticError('Solver control did not pass')
    print('complete',result['all_six_approach_same_reference_orbit'],flush=True)

if __name__=='__main__':main()
