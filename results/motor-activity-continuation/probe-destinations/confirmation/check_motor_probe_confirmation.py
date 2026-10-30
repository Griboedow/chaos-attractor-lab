"""Independent full Radau matrix and phase/permutation controls for probe cycles."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment, minimize_scalar
from scipy.sparse import eye, kron, csc_matrix
import motor_rows

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/motor-activity-continuation/probe-destinations'
REFERENCE=ROOT/'results/motor-activity-continuation/lane49/periodic/cycle.npz'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(state): return np.stack((state[:8],state[8:16],np.r_[state[16:],-sum(state[16:])]),axis=1)


def groups(state):
    r=rows(state);todo=set(range(8));result=[]
    while todo:
        a=min(todo);g=[b for b in sorted(todo) if np.linalg.norm(r[a]-r[b])<1e-10]
        result.append(g);todo-=set(g)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('Choose a new output')
    args.output.mkdir(parents=True)
    sources=[Path(__file__),Path(motor_rows.__file__)]
    for path in sources:(args.output/path.name).write_bytes(path.read_bytes())
    protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_sha256={p.name:sha(p) for p in sources},scope='Numerical confirmation, not an interval certificate.',
        reference_path=REFERENCE.relative_to(ROOT).as_posix(),reference_sha256=sha(REFERENCE))
    (args.output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    reference=np.load(REFERENCE);p=reference['parameters'];period=float(reference['time'][-1])
    ref=solve_ivp(lambda t,y:motor_rows.rhs(t,y,p,8),(0,period),reference['full_states'][0],
        method='DOP853',rtol=5e-13,atol=5e-15,max_step=.01,dense_output=True)
    if not ref.success: raise RuntimeError(ref.message)
    comparisons=[]
    for case in (16,17,18,19):
        path=DATA/f'case-{case}/periodic/cycle.npz';target_file=np.load(path)
        target=rows(target_file['full_states'][0])
        def fit(t):
            state=ref.sol(float(t)%period);r=rows(state)
            cost=np.sum((target[:,None,:]-r[None,:,:])**2,axis=2)
            a,b=linear_sum_assignment(cost)
            return float(np.mean(cost[a,b])),state,b
        grid=np.linspace(0,period,1025)[:-1]
        center=float(grid[np.argmin([fit(t)[0] for t in grid])])
        solution=minimize_scalar(lambda t:fit(t)[0],bounds=(center-period/512,center+period/512),
                                 method='bounded',options=dict(xatol=1e-14,maxiter=200))
        phase=float(solution.x);history=[]
        for _ in range(12):
            value,state,assignment=fit(phase)
            # Gauss-Newton scalar phase correction avoids bounded optimizer's
            # relative sqrt(eps) stopping floor near a nonzero phase.
            displacement=rows(state)[assignment]-target
            velocity=rows(motor_rows.rhs(0,state,p,8))[assignment]
            correction=float(np.sum(displacement*velocity)/np.sum(velocity*velocity))
            history.append(dict(phase=phase,rms=float(np.sqrt(value)),correction=correction))
            trial=phase-correction
            if fit(trial)[0]>value: break
            phase=trial
            if abs(correction)<1e-14: break
        value,state,assignment=fit(phase)
        record=dict(case=case,cycle_sha256=sha(path),period=float(target_file['time'][-1]),
            inferred_root_groups=groups(target_file['full_states'][0]),reference_period=period,
            phase=phase%period,target_row_to_reference_row=assignment.tolist(),
            full_row_rms_distance=float(np.sqrt(value)),phase_refinement=history,
            limitation='Floating-point geometric comparison; the assignment is an observed symmetry relation, not a certified equivalence.')
        comparisons.append(record)
    (args.output/'orbit-comparisons.json').write_text(json.dumps(comparisons,indent=2)+'\n',encoding='utf-8')
    print('comparisons',[(v['case'],v['full_row_rms_distance']) for v in comparisons],flush=True)
    start=time.perf_counter();path=DATA/'case-16/periodic/cycle.npz';source=np.load(path)
    p=source['parameters'];period=float(source['time'][-1]);initial=source['full_states'][0]
    flow=solve_ivp(lambda t,y:motor_rows.rhs(t,y,p,8),(0,period),initial,
        jac=lambda t,y:motor_rows.field_jac(y,p,8)[1],method='Radau',rtol=2e-12,atol=2e-14,
        max_step=.015,dense_output=True)
    if not flow.success: raise RuntimeError(flow.message)
    def jac_at(t): return motor_rows.field_jac(flow.sol(t),p,8)[1]
    def rhs(t,z): return (jac_at(t)@z.reshape((23,23),order='F')).ravel(order='F')
    def linear_jac(t,z): return kron(eye(23,format='csc'),csc_matrix(jac_at(t)),format='csc')
    test=np.sin(np.arange(23*23)+.37);test_error=float(np.max(abs(rhs(.37,test)-linear_jac(.37,test)@test)))
    if test_error>1e-12: raise ArithmeticError('Column-major sparse Jacobian control failed')
    variation=solve_ivp(rhs,(0,period),np.eye(23).ravel(order='F'),jac=linear_jac,
        method='Radau',rtol=2e-12,atol=2e-14,max_step=.02)
    if not variation.success: raise RuntimeError(variation.message)
    matrix=variation.y[:,-1].reshape((23,23),order='F')
    eig=np.linalg.eigvals(matrix);old=np.linalg.eigvals(source['full_monodromy'])
    a,b=linear_sum_assignment(abs(eig[:,None]-old[None,:]))
    flow_index=int(np.argmin(abs(eig-1)));remaining=np.delete(eig,flow_index)
    field=motor_rows.rhs(0,initial,p,8)
    record=dict(status='full_radau_spectrum_computed',input_path=path.relative_to(ROOT).as_posix(),input_sha256=sha(path),
        period=period,state_method='Full 23D Radau with analytic Jacobian and dense output; no cluster projection.',
        variational_method='All 23 columns, Radau on the time-dependent full Jacobian along the independently computed Radau state.',
        settings=dict(rtol=2e-12,atol=2e-14,state_max_step=.015,matrix_max_step=.02),
        sparse_column_convention_error=test_error,state_closure_max=float(np.max(abs(flow.y[:,-1]-initial))),
        matrix_difference_from_DOP853=float(np.max(abs(matrix-source['full_monodromy']))),
        eigenvalue_matching_max_error=float(np.max(abs(eig[a]-old[b]))),
        flow_identity_max_error=float(np.max(abs(matrix@field-field))),
        flow_multiplier_real=float(eig[flow_index].real),flow_multiplier_imag=float(eig[flow_index].imag),
        maximum_nonflow_modulus=float(max(abs(remaining))),
        unstable_count_at_1e_minus_6=int(np.sum(abs(eig)>1+1e-6)),
        neutral_count_at_1e_minus_6=int(np.sum(abs(abs(eig)-1)<1e-6)),
        multipliers=[dict(real=float(z.real),imag=float(z.imag),modulus=float(abs(z))) for z in sorted(eig,key=abs,reverse=True)],
        wall_seconds=time.perf_counter()-start,
        limitations=['The state interpolant and matrix integration are numerical, not interval enclosures.',
                     'This supports a strongly attracting periodic candidate, not a computer-assisted stability theorem or historical novelty.'])
    np.savez_compressed(args.output/'case-16-radau.npz',monodromy=matrix,eigenvalues=eig,parameters=p,
                        state_time=flow.t,state=flow.y.T,period=period)
    (args.output/'case-16-radau.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('full Radau nonflow maximum',record['maximum_nonflow_modulus'],'matrix error',record['matrix_difference_from_DOP853'],flush=True)


if __name__=='__main__':main()
