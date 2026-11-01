"""Numerical pseudo-arclength continuation of the nearby four-row cycle.

This locates candidate organizing bifurcations, not an interval certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from check_motor_asymmetry import augmented

FREE=np.array([i for i in range(11) if i!=8]);B_SCALE=100.


@njit(cache=True)
def extended(z,p):
    base=augmented(z[:133],p,4)
    probe=np.zeros(133);probe[:11]=z[:11];probe[11:132]=np.eye(11).ravel()
    jac=augmented(probe,p,4)[11:132].reshape(11,11)
    derivative=jac@z[133:144]
    for i in range(4,8):derivative[i]+=1.
    return np.r_[base,derivative]


def shoot(y,rho,r):
    q=np.zeros(11);q[FREE]=y[:10];period=y[10];b=y[11]/B_SCALE
    if not (2<period<30 and .8<b<2.2):raise ValueError('outside declared cycle continuation box')
    p=np.array([rho,r,b]);z=np.r_[q,np.eye(11).ravel(),0.,np.zeros(11)]
    sol=solve_ivp(lambda t,z:extended(z,p),(0,period),z,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.02)
    if not sol.success:raise RuntimeError(sol.message)
    end=sol.y[:11,-1];matrix=sol.y[11:132,-1].reshape(11,11)
    probe=np.r_[end,np.eye(11).ravel(),0.];velocity=augmented(probe,p,4)[:11]
    derivative=np.column_stack([(matrix-np.eye(11))[:,FREE],velocity,sol.y[133:144,-1]/B_SCALE])
    return end-q,derivative,matrix,q,period,p


def correct(y,rho,r,predictor=None,tangent=None,fixed_b=None):
    y=y.copy();history=[]
    for it in range(15):
        value,j,m,q,period,p=shoot(y,rho,r)
        constraint=(y[11]-fixed_b*B_SCALE) if fixed_b is not None else np.dot(y-predictor,tangent)
        row=np.r_[np.zeros(11),1.] if fixed_b is not None else tangent
        residual=np.r_[value,constraint];err=float(np.max(abs(residual)));history.append(err)
        if err<2e-10:
            assert np.linalg.norm(q[8:])>.1
            return y,dict(iterations=it,residual=err,residual_history=history,monodromy=m,seed=q,period=period,parameters=p)
        step=np.linalg.solve(np.vstack([j,row]),-residual)
        accepted=False
        for factor in [1.,.5,.25,.125,.0625]:
            trial=y+factor*step
            try:
                val=shoot(trial,rho,r)[0]
                con=(trial[11]-fixed_b*B_SCALE) if fixed_b is not None else np.dot(trial-predictor,tangent)
                if np.max(abs(np.r_[val,con]))<err:y=trial;accepted=True;break
            except (ValueError,RuntimeError):pass
        if not accepted:raise RuntimeError(f'Newton line search failed, residual {err}')
    raise RuntimeError('Newton iteration limit')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);began=time.perf_counter()
    d=np.load(args.source);p=d['parameters'];rho,r,_=p;q=d['seed'];period=float(d['period'])
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        method='Full11D shooting with initial X0=0; analytic state, period and b derivatives; pseudo-arclength continuation.',
        clock='original tau',rtol=2e-11,atol=2e-13,maxstep=.02,phase_index=8,b_scale=B_SCALE,
        initial_b=[1.556,1.555],steps=args.steps,initial_arclength_step=.01,minimum_arclength_step=1e-5,
        scope='Numerical cycle branch and multipliers only; a +1 or -1 near-crossing is not a validated bifurcation theorem.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['continue_motor_braid_cycles.py','check_motor_asymmetry.py','motor_rows.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[];solutions=[];cycle_data=[];attempts=[]
    def save(y,record):
        m=record.pop('monodromy');seed=record.pop('seed');parameters=record.pop('parameters');period=record['period']
        eigen=np.linalg.eigvals(m);flow=int(np.argmin(abs(eigen-1)));nonflow=np.delete(eigen,flow);nonflow=nonflow[np.argsort(abs(nonflow))[::-1]]
        record.update(index=len(records),b=float(parameters[2]),parameters=parameters.tolist(),
            flow_multiplier=[float(eigen[flow].real),float(eigen[flow].imag)],
            nonflow_multipliers=[[float(z.real),float(z.imag)] for z in nonflow],
            unstable_count=int((abs(nonflow)>1+1e-7).sum()),largest_nonflow_abs=float(abs(nonflow[0])))
        records.append(record);solutions.append(y.copy());cycle_data.append((seed.copy(),period,m.copy()))
        (out/'analysis.json').write_text(json.dumps(dict(completed=len(records),records=records,attempts=attempts,wall_seconds=time.perf_counter()-began),indent=2)+'\n')
        np.savez_compressed(out/'cycles.npz',coordinates=solutions,seeds=np.array([c[0] for c in cycle_data]),periods=np.array([c[1] for c in cycle_data]),monodromies=np.array([c[2] for c in cycle_data]))
    y=np.r_[q[FREE],period,1.556*B_SCALE]
    y,record=correct(y,rho,r,fixed_b=1.556);save(y,record)
    y[11]=1.555*B_SCALE;y,record=correct(y,rho,r,fixed_b=1.555);save(y,record)
    ds=.01
    for i in range(args.steps):
        tangent=solutions[-1]-solutions[-2];tangent/=np.linalg.norm(tangent)
        while True:
            predicted=solutions[-1]+ds*tangent
            try:
                y,record=correct(predicted,rho,r,predictor=predicted,tangent=tangent)
                record['arclength_step']=ds;save(y,record);break
            except (RuntimeError,ValueError,np.linalg.LinAlgError) as error:
                attempts.append(dict(step=i,ds=ds,error=str(error)));ds*=.5
                if ds<1e-5:
                    (out/'failure.json').write_text(json.dumps(dict(reason=str(error),step=i,ds=ds),indent=2)+'\n');return
        if (i+1)%10==0:print(i+1,record['b'],record['largest_nonflow_abs'],flush=True)
        ds=min(.01,ds*1.2)
    print('completed',len(records),'minb',min(r['b'] for r in records),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=160);main(p.parse_args())
