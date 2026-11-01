"""Leading tangent growth of the implicit BDF rod map, including its history."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
import spatial_accelerated as m


@njit(cache=True)
def curvature_force_variation(q,dk,beta):
    vectors,mids=m.strains_midframes(q)
    grad=np.zeros((len(q),3))
    for i in range(len(q)):
        moment=-beta*np.array([dk[i,0],dk[i,1],0.0])
        grad[i]+=m.inv_left_jac_action(vectors[i],moment,-1)
        if i>0:grad[i-1]-=m.inv_left_jac_action(vectors[i],moment,1)
    return grad


@njit(cache=True)
def drag_velocity_variation(q,dv,zeta,xi):
    frames,nodes=m.geometry(q);n=len(q);h=1.0/n
    grad=np.zeros((n,3));tail=np.zeros(3)
    for i in range(n-1,-1,-1):
        r=frames[i];drag=(r*zeta)@r.T*xi
        force=tail+h*drag@(dv[i]/6+dv[i+1]/3)
        body=r.T@force
        grad[i]=h*np.array([-body[1],body[0],0.0])
        tail+=h*drag@((dv[i]+dv[i+1])/2)
    return grad


def node_variation(q,theta):
    frames=m.geometry(q)[0]
    local=np.zeros_like(theta);local[:,0]=theta[:,1];local[:,1]=-theta[:,0]
    delta=np.einsum('nij,njp->nip',frames,local)/len(q)
    return np.concatenate([np.zeros((1,3,theta.shape[-1])),np.cumsum(delta,axis=0)])


def propagate(qnew,kold,kprev,nold,nprev,dkold,dkprev,dnold,dnprev,
              beta,zeta,xi,chi,eta,dt,order,fd_step=1e-6):
    n=len(qnew);count=dkold.shape[-1]
    alpha=1.0 if order==1 else 1.5
    denominator=(alpha+dt*np.array([1.,eta**-2]))[None,:,None]
    dkfree=dkold/denominator if order==1 else (2*dkold-.5*dkprev)/denominator
    dv=-dnold/dt if order==1 else (-2*dnold+.5*dnprev)/dt
    rhs=np.empty((3*n,count))
    for j in range(count):
        force=curvature_force_variation(qnew,dkfree[:,:,j],beta)+drag_velocity_variation(qnew,dv[:,:,j],zeta,xi)
        rhs[:,j]=-force.ravel()/beta[0]
    jac=m.mechanical_jacobian(qnew,kold,kprev,nold,nprev,beta,zeta,xi,chi,eta,dt,order)
    theta=np.linalg.solve(jac,rhs).reshape(n,3,count)
    dk=dkfree.copy()
    for j in range(count):
        direction=theta[:,:,j]
        eps=fd_step/max(1.,float(np.max(abs(direction))))
        plus=m.target_curvature(m.rotate_increment(qnew,eps*direction),chi,eta)
        minus=m.target_curvature(m.rotate_increment(qnew,-eps*direction),chi,eta)
        dk[:,:,j]+=dt*(plus-minus)/(2*eps)/denominator[0,:,0]
    return theta,dk,node_variation(qnew,theta)


def run(args):
    seed_report=json.loads((args.seed/'report.json').read_text(encoding='utf-8'))
    parameters=seed_report['parameters']
    assert seed_report['success']
    with np.load(args.seed/'checkpoint.npz') as c:
        q=c['q'].copy();k=c['k'].copy();nodes=c['nodes'].copy()
        kp=c['k_previous'].copy();npv=c['nodes_previous'].copy();predictor=c['predictor'].copy()
        dt=float(c['dt']);origin=float(c['time'])
        saved_tangent={name:c[name].copy() for name in ['dk','dk_previous','dnodes','dnodes_previous','log_growth']} if args.resume_tangent else None
    beta,zeta,xi,chi,eta=m.parameters(lam=parameters['slenderness'],eta=parameters['eta'],
        xi=parameters['xi'],chi=parameters['chi'],nu=parameters['nu'])
    assert abs(chi)>0
    n=len(q);count=3;rng=np.random.default_rng(args.random_seed)
    # Generic curvature-history directions. Their initial algebraic mismatch
    # is resolved by the implicit step; the declared burn interval is excluded.
    dk=rng.normal(size=(n,2,count))*abs(chi)
    dkprev=dk.copy();dn=np.zeros((n+1,3,count));dnprev=dn.copy()
    initial_observed=0.0
    if args.resume_tangent:
        assert args.burn==0
        dk=saved_tangent['dk'];dkprev=saved_tangent['dk_previous']
        dn=saved_tangent['dnodes'];dnprev=saved_tangent['dnodes_previous']
        initial_observed=seed_report['observed_duration']
    steps=round(args.duration/dt);burn=round(args.burn/dt)
    assert abs(steps*dt-args.duration)<1e-9 and steps>burn>=0
    totals=saved_tangent['log_growth'].copy() if args.resume_tangent else np.zeros(count)
    observations=[];success=True;failure=None
    start=time.perf_counter();worst=0.;max_hinge=0.;completed=0
    sample_stride=max(1,round(.012/dt))
    times=[origin];qs=[q.copy()];ks=[k.copy()];tips=[nodes[-1].copy()]
    for step in range(1,steps+1):
        qn,kn,inc,its,ok,res,hinge=m.implicit_step(q,k,nodes,kp,npv,beta,zeta,xi,chi,eta,dt,2,predictor,args.tolerance)
        worst=max(worst,res);max_hinge=max(max_hinge,hinge)
        if not ok:
            success=False;failure=dict(step=step,time=origin+step*dt,residual=res);break
        theta,dkn,dnn=propagate(qn,k,kp,nodes,npv,dk,dkprev,dn,dnprev,beta,zeta,xi,chi,eta,dt,2,args.fd_step)
        vector=np.concatenate([theta.reshape(3*n,count),dkn.reshape(2*n,count)/abs(chi)],axis=0)/np.sqrt(n)
        basis,upper=np.linalg.qr(vector,mode='reduced')
        diagonal=np.abs(np.diag(upper))
        assert np.all(diagonal>1e-100) and np.all(np.isfinite(diagonal))
        inverse=np.linalg.solve(upper,np.eye(count))
        dkprev=np.einsum('nip,pq->niq',dk,inverse)
        dnprev=np.einsum('nip,pq->niq',dn,inverse)
        dk=np.einsum('nip,pq->niq',dkn,inverse)
        dn=np.einsum('nip,pq->niq',dnn,inverse)
        kp,npv=k.copy(),nodes.copy();q,k=qn,kn;nodes=m.geometry(q)[1];predictor=inc;completed=step
        if step>burn:totals+=np.log(diagonal)
        if step%sample_stride==0 or step==steps:
            times.append(origin+step*dt);qs.append(q.copy());ks.append(k.copy());tips.append(nodes[-1].copy())
        if step%1000==0 or step==steps:
            row=dict(time=origin+step*dt,observed_duration=initial_observed+max(0,(step-burn)*dt),
                exponents=(totals/(initial_observed+(step-burn)*dt)).tolist() if step>burn else None,
                elapsed=time.perf_counter()-start,residual=res)
            observations.append(row);print(json.dumps(row),flush=True)
    np.savez_compressed(args.output/'trajectory.npz',t=np.array(times),q=np.array(qs),k=np.array(ks),tip=np.array(tips))
    np.savez_compressed(args.output/'checkpoint.npz',q=q,k=k,nodes=nodes,k_previous=kp,
        nodes_previous=npv,predictor=predictor,time=origin+completed*dt,dt=dt,
        dk=dk,dk_previous=dkprev,dnodes=dn,dnodes_previous=dnprev,log_growth=totals)
    result=dict(parameters=parameters,segments=n,dt=dt,time_origin=origin,requested_duration=args.duration,
        completed_time=float(times[-1]),final_accepted_time=origin+completed*dt,success=success,
        failure=failure,plane='spatial',order=2,burn_duration=burn*dt,
        initial_observed_duration=initial_observed,new_observed_duration=max(0,(completed-burn)*dt),
        observed_duration=initial_observed+max(0,(completed-burn)*dt),
        exponents=(totals/(initial_observed+(completed-burn)*dt)).tolist() if completed>burn else None,
        max_scaled_residual=worst,max_hinge_angle=max_hinge,newton_tolerance=args.tolerance,
        drive_directional_difference_step=args.fd_step,random_seed=args.random_seed,
        initial_tangents='Restored current and previous tangent histories and cumulative log growth' if args.resume_tangent else 'Three independent Gaussian curvature directions; algebraic startup excluded by burn',
        qr_norm='Mean squared body orientation variation plus squared curvature variation divided by chi squared; the same QR transform is applied to both BDF history states.',
        elapsed_seconds=time.perf_counter()-start,observations=observations,
        seed_checkpoint=dict(external_path=str(args.seed/'checkpoint.npz'),sha256=hashlib.sha256((args.seed/'checkpoint.npz').read_bytes()).hexdigest()),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        solver_sha256=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest(),
        scope='Three finite-time exponents of the constant-step implicit BDF2 map. Derivatives of old-state forcing are analytic; the curvature-drive directional derivative is centered numerically. No continuous-time, asymptotic or continuum chaos certificate.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ['parameters','segments','dt','success','observed_duration','exponents','elapsed_seconds']}),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--duration',type=float,default=120)
    ap.add_argument('--burn',type=float,default=24)
    ap.add_argument('--resume-tangent',action='store_true')
    ap.add_argument('--random-seed',type=int,default=260409490)
    ap.add_argument('--fd-step',type=float,default=1e-6)
    ap.add_argument('--tolerance',type=float,default=1e-8)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);run(args)
