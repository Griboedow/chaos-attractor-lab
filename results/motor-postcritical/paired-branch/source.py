"""Continue the numerically doubled four-pair motor cycle and its Floquet blocks."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag
from scipy.optimize import linear_sum_assignment

from motor_postcritical_cycles import field_jac, lift, unpack, integrate, shoot, audit, ROOT
from motor_rows import rhs as full_rhs, field_jac as full_field_jac

GROUPS=[[0,1],[2,3],[4,5],[6,7]]
WEIGHTS=np.ones(4)/4


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def symmetry(source, output):
    if output.exists():
        raise FileExistsError(output)
    a=np.load(source/'cycle.npz');p=a['parameters'];state=a['quotient_initial_state'];period=a['time'][-1]
    q=integrate(state, period/4, p, WEIGHTS, augmented=True, rtol=1e-12, max_step=.01)
    half=integrate(state, period/2, p, WEIGHTS, augmented=True, rtol=1e-12, max_step=.01)
    permutation=np.array([2,3,1,0])
    def transform(z):
        c,s,x=unpack(z, WEIGHTS)
        return np.r_[c[permutation],s[permutation],x[permutation][:-1]]
    pm=np.column_stack([transform(e) for e in np.eye(11)])
    quarter=q.y[11:,-1].reshape(11,11)
    twisted=np.linalg.solve(pm,quarter)
    mon=a['quotient_monodromy']
    full=lift(state,GROUPS)
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                parameters=p.tolist(),period=float(period),quarter_period=float(period/4),
                initial_state=state.tolist(),state_layout='c0..c3,s0..s3,X0..X2; X3=-sum(X0..X2)',
                phase='X0=0, upward',phase_velocity=float(field_jac(state,p,WEIGHTS)[0][8]),
                permutation_new_from_old=permutation.tolist(),permutation_matrix=pm.tolist(),
                quarter_endpoint_residual=float(np.max(abs(q.y[:11,-1]-transform(state)))),
                half_endpoint_residual=float(np.max(abs(half.y[:11,-1]-transform(transform(state))))),
                fourth_power_twisted_derivative_error=float(np.max(abs(np.linalg.matrix_power(twisted,4)-mon))),
                quarter_fundamental_matrix=quarter.tolist(),twisted_quarter_derivative=twisted.tolist(),
                source=str(source.relative_to(ROOT)).replace('\\','/'),
                source_sha256=hashlib.sha256((source/'cycle.npz').read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                caveat='Numerical controls supporting an exact equivariant formulation, not an interval certificate.')
    write(output,record)
    print(json.dumps(record,indent=2))


def augmented(state,period,p,dense=False):
    def fun(t,z):
        f,j=field_jac(z[:11],p,WEIGHTS)
        c,s,x=unpack(z[:11],WEIGHTS);v=-p[0]*x+p[1]*(s-np.mean(s))
        values=[f,(j@z[11:132].reshape(11,11)).ravel()]
        for k in range(4):
            mat=np.array([[-1,-v[k]-p[1]*s[k],p[0]*s[k]],
                          [v[k],-1+p[1]*c[k],-p[0]*c[k]],[0,p[1],-p[0]]])
            values.append((mat@z[132+9*k:141+9*k].reshape(3,3)).ravel())
        return np.concatenate(values)
    initial=np.r_[state,np.eye(11).ravel(),np.tile(np.eye(3).ravel(),4)]
    sol=solve_ivp(fun,(0,period),initial,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.025,dense_output=dense)
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def spectral(state,period,p,keep_cycle=True):
    sol=augmented(state,period,p,dense=keep_cycle)
    quotient=sol.y[11:132,-1].reshape(11,11)
    blocks=[sol.y[132+9*k:141+9*k,-1].reshape(3,3) for k in range(4)]
    eigc=np.linalg.eigvals(quotient)
    neutral=np.argmin(abs(eigc-1));nonneutral=np.delete(eigc,neutral)
    eigb=[np.linalg.eigvals(q) for q in blocks]
    serial=lambda es:[dict(real=float(e.real),imag=float(e.imag),modulus=float(abs(e))) for e in es[np.argsort(-abs(es))]]
    record=dict(quotient_multipliers=serial(eigc),transverse_multipliers=[serial(e) for e in eigb],
                quotient_neutral_error=float(abs(eigc[neutral]-1)),
                largest_quotient_nonflow_modulus=float(max(abs(nonneutral))),
                largest_transverse_modulus=float(max(max(abs(e)) for e in eigb)),
                unstable_quotient_count=int(np.sum(abs(nonneutral)>1+1e-6)),
                unstable_transverse_count=int(sum(np.sum(abs(e)>1+1e-6) for e in eigb)),
                closure_error=float(np.max(abs(sol.y[:11,-1]-state))),
                transverse_det_I_plus_M=[float(np.linalg.det(np.eye(3)+q)) for q in blocks])
    data=dict(quotient_monodromy=quotient,transverse_monodromy=np.array(blocks),parameters=p,
              quotient_initial_state=state)
    if keep_cycle:
        times=np.linspace(0,period,1001);sample=sol.sol(times)[:11].T
        diff=np.linalg.norm(sample[:,[0,4,8]]-sample[:,[1,5,9]],axis=1)
        record['split_pair_rms_distance']=float(np.sqrt(np.mean(diff**2)))
        record['split_pair_max_distance']=float(max(diff))
        record['max_disk_radius_squared']=float(np.max(sample[:,:4]**2+sample[:,4:8]**2))
        data.update(time=times,cluster_states=sample,full_states=np.array([lift(x,GROUPS) for x in sample]))
    return record,data


def run(source,output):
    source,output=source.resolve(),output.resolve()
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    old=np.load(source/'cycle.npz');p0=old['parameters'];s0=old['quotient_initial_state'];t0=old['time'][-1]
    star=3.3702183995790422
    grids=[('near_onset',star-np.array([.005,.002,.001,.0005,.0002,.0001])),
           ('toward_chaos',np.arange(p0[0],3.225,-.005))]
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                source=source.relative_to(ROOT).as_posix(),source_sha256=hashlib.sha256((source/'cycle.npz').read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                model='Exact N8 restriction to four clusters of size two',groups=GROUPS,
                rtol=2e-11,atol=2e-13,max_step=.025,cases=[],failures=[],
                caveat='Numerical branch continuation, not an interval bifurcation theorem.')
    began=time.perf_counter()
    for label,grid in grids:
        state,period=s0.copy(),float(t0)
        for index,rho in enumerate(grid):
            p=np.array([rho,p0[1],p0[2]])
            try:
                state,period,sol,history=shoot(state,period,p,WEIGHTS)
                spec,data=spectral(state,period,p)
                name=f'{label}-{index:02d}'
                np.savez_compressed(output/(name+'.npz'),**data)
                item=dict(label=label,index=index,rho=float(rho),delta=float((p[1]-1)**2/rho**2-1),
                           offset_below_first_crossing=float(star-rho),period=period,
                           initial_state=state.tolist(),shooting=history,data_file=name+'.npz',**spec)
                record['cases'].append(item)
                print(label,index,rho,period,spec['largest_quotient_nonflow_modulus'],spec['largest_transverse_modulus'],flush=True)
            except Exception as error:
                record['failures'].append(dict(label=label,index=index,rho=float(rho),error=repr(error)))
                write(output/'analysis.json',record)
                raise
            record['wall_seconds']=time.perf_counter()-began
            write(output/'analysis.json',record)
    (output/'source.py').write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['symmetry','continue'])
    parser.add_argument('--source',type=Path,default=ROOT/'results/motor-postcritical/paired-cycle-rho-3p3652184-retry')
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    if a.command=='symmetry':symmetry(a.source.resolve(),a.output.resolve())
    else:run(a.source,a.output)
