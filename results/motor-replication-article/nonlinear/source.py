"""Predeclared nonlinear N=16 perturbations for the replication manuscript.

All eight starts and both integration settings are retained. These ordinary
numerical trajectories illustrate the separate interval theorem.
"""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
import motor_rows

ROOT=Path(__file__).resolve().parents[1]
GROUPS=((0,1,2,3),(4,5,6,7),(8,9),(10,11,12,13),(14,15))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mid(value):
    return math.ldexp(float(value['midpoint'][0]),value['midpoint'][1])


def physical(states):
    return np.concatenate((states,-np.sum(states[...,32:],axis=-1,keepdims=True)),axis=-1)


def write(path,data):
    path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'results/motor-five-cycle-proof/recheck-01.json'
    block_source=ROOT/'results/motor-five-replication-proof/singleton-blocks.json'
    cert=json.loads(source.read_text())
    blocks=json.loads(block_source.read_text())
    assert cert['status']=='attracting_five_cluster_cycle_certified'
    q=np.array([mid(x) for x in cert['initial_state_box']])
    period=mid(cert['positive_period'])
    params=np.array([129/55,1+.95/np.sqrt(32001/500000),1.])
    weights=np.array([2,2,1,2,1])/8
    x=np.r_[q[10:],-np.dot(weights[:4],q[10:])/weights[-1]]
    rows=np.empty((16,3))
    for j,group in enumerate(GROUPS):
        rows[list(group)]=[q[j],q[5+j],x[j]]
    initial=np.r_[rows[:,0],rows[:,1],rows[:-1,2]]
    field=lambda t,z:motor_rows.rhs(t,z,params,16)
    declared=[]
    for cluster in (2,4):
        matrix=np.array([[mid(x) for x in row] for row in blocks['singleton_blocks'][str(cluster)]['full_return']])
        values,vectors=np.linalg.eig(matrix)
        j=np.argmax(abs(values));mu=float(values[j].real)
        assert abs(values[j].imag)<1e-12 and mu < -1
        vector=vectors[:,j].real
        vector/=np.linalg.norm(vector)
        if vector[np.argmax(abs(vector))]<0:
            vector=-vector
        for amplitude in (1e-6,1e-8):
            for sign in (-1,1):
                perturb=np.zeros((16,3))
                a,b=GROUPS[cluster]
                perturb[a]=sign*amplitude*vector/np.sqrt(2)
                perturb[b]=-perturb[a]
                delta=np.r_[perturb[:,0],perturb[:,1],perturb[:-1,2]]
                declared.append(dict(cluster=cluster,amplitude=amplitude,sign=sign,
                                     vector=vector.tolist(),mu=mu,start=(initial+delta).tolist()))
    settings=[dict(name='main',rtol=2e-11,atol=2e-13,max_step=.04),
              dict(name='strict',rtol=2e-12,atol=2e-14,max_step=.02)]
    protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  period=period,parameters=params.tolist(),groups=GROUPS,horizon_periods=4,
                  certificate_sha256=sha(source),matrix_sha256=sha(block_source),
                  source_sha256=sha(__file__),field_sha256=sha(motor_rows.__file__),
                  settings=settings,starts=declared,
                  numerical_scope='Full47D nonlinear integrations. Eigenvectors of interval-matrix midpoints propose initial directions only. The interval proof is separate.',
                  gates=dict(max_state_disagreement=1e-7,small_amplitude_stroboscopic_relative_error_first_three_returns=.01))
    write(out/'protocol.json',protocol)
    (out/'source.py').write_bytes(Path(__file__).read_bytes())
    (out/'motor_rows.py').write_bytes(Path(motor_rows.__file__).read_bytes())
    times=np.linspace(0,4*period,2401)
    records=[]
    for index,start in enumerate(declared):
        solutions=[]
        for setting in settings:
            kwargs={k:v for k,v in setting.items() if k!='name'}
            sol=solve_ivp(field,(0,times[-1]),start['start'],method='DOP853',t_eval=times,**kwargs)
            assert sol.success
            solutions.append(sol.y.T)
        strict=physical(solutions[1])
        a,b=GROUPS[start['cluster']]
        contrasts=(strict[:,[a,16+a,32+a]]-strict[:,[b,16+b,32+b]])/np.sqrt(2)
        signed=contrasts@np.array(start['vector'])
        contrast_norm=np.linalg.norm(contrasts,axis=1)
        indices=np.arange(0,len(times),600)
        prediction=start['sign']*start['amplitude']*start['mu']**np.arange(5)
        relative=signed[indices]/prediction-1
        disagreement=float(np.max(abs(solutions[0]-solutions[1])))
        path=out/f'trajectory-{index:02}.npz'
        np.savez_compressed(path,time=times,states=solutions[1],main_states=solutions[0],
                            signed_contrast=signed,contrast_norm=contrast_norm,
                            returns=times[indices],observed=signed[indices],prediction=prediction)
        record=dict(index=index,cluster=start['cluster'],amplitude=start['amplitude'],sign=start['sign'],
                    max_state_disagreement=disagreement,solver_control_passed=disagreement<1e-7,
                    relative_return_errors=relative.tolist(),
                    first_three_small_signal_returns_passed=bool(np.max(abs(relative[1:4]))<.01) if start['amplitude']==1e-8 else None,
                    final_contrast_norm=float(contrast_norm[-1]),trajectory_sha256=sha(path))
        records.append(record)
        write(out/'progress.json',dict(runs=records))
        print(index,start['cluster'],start['amplitude'],start['sign'],disagreement,relative[1:4],flush=True)
    # A different solver checks one complete nonlinear history, including the
    # larger perturbation and departure from its small-signal approximation.
    probe=declared[0]
    other=solve_ivp(field,(0,times[-1]),probe['start'],method='Radau',t_eval=times,
                    jac=lambda t,z:motor_rows.field_jac(z,params,16)[1],
                    rtol=2e-11,atol=2e-13,max_step=.04)
    assert other.success
    dop=np.load(out/'trajectory-00.npz')['states']
    independent_error=float(np.max(abs(other.y.T-dop)))
    np.savez_compressed(out/'radau-control.npz',time=times,states=other.y.T)
    passed=all(r['solver_control_passed'] for r in records) and all(
        r['first_three_small_signal_returns_passed'] for r in records if r['amplitude']==1e-8) and independent_error<1e-7
    write(out/'analysis.json',dict(status='passed' if passed else 'failed',runs=records,
                                  independent_radau_max_state_difference=independent_error,
                                  radau_sha256=sha(out/'radau-control.npz'),
                                  qualification='Numerical illustrations, not a validated nonlinear escape or basin theorem. All declared starts are retained.'))
    print('passed' if passed else 'failed','Radau',independent_error,flush=True)
    if not passed:
        raise SystemExit(2)


if __name__=='__main__':
    main()
