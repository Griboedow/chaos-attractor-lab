"""Independent endpoint-initialized rod IVPs, with exact BDF restart records."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import spatial_fast as m


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def remesh(q,k,n):
    old=len(q)
    if n==old:return q.copy(),k.copy()
    source_s=np.r_[0.0,(np.arange(old)+.5)/old]
    source_q=np.vstack([np.array([1.,0,0,0]),q])
    target_s=(np.arange(n)+.5)/n
    out=np.empty((n,4))
    # Piecewise geodesic interpolation; the last interval also supplies the
    # short extrapolation to a finer mesh's last orientation center.
    for j,s in enumerate(target_s):
        left=min(len(source_s)-2,max(0,np.searchsorted(source_s,s)-1))
        fraction=(s-source_s[left])/(source_s[left+1]-source_s[left])
        qa=source_q[left];qb=source_q[left+1]
        delta=m.logq(m.mul(np.r_[qa[0],-qa[1:]],qb))
        out[j]=m.mul(qa,m.expq(fraction*delta))
    # Internal curvatures live at bending-interval midpoints, not at
    # orientation centers. Endpoint extension is constant for this initial
    # interpolation only; it is not a boundary condition during evolution.
    source_k=np.r_[.25/old,np.arange(1,old)/old]
    target_k=np.r_[.25/n,np.arange(1,n)/n]
    kout=np.column_stack([np.interp(target_k,source_k,k[:,j]) for j in range(2)])
    return out,kout


def consistent_twist(q,k,beta,tolerance=1e-10):
    n=len(q);spin=np.zeros((n,3));initial=None
    for iteration in range(12):
        current=m.rotate_increment(q,spin)
        res=m.energy_gradient(current,k,beta)[1][:,2]/beta[0]
        norm=float(np.max(abs(res)))
        if initial is None:initial=norm
        if norm<tolerance:
            return current,dict(initial_scaled_torque=initial,final_scaled_torque=norm,
                iterations=iteration,maximum_material_spin=float(np.max(abs(spin[:,2]))))
        jac=np.empty((n,n))
        for j in range(n):
            delta=np.zeros((n,3));delta[j,2]=1e-6
            plus=m.energy_gradient(m.rotate_increment(current,delta),k,beta)[1][:,2]/beta[0]
            minus=m.energy_gradient(m.rotate_increment(current,-delta),k,beta)[1][:,2]/beta[0]
            jac[:,j]=(plus-minus)/2e-6
        update=np.linalg.solve(jac,-res);accepted=False
        for ls in range(12):
            trial=spin.copy();trial[:,2]+=update*2.0**(-ls)
            newres=m.energy_gradient(m.rotate_increment(q,trial),k,beta)[1][:,2]/beta[0]
            if np.max(abs(newres))<norm:
                spin=trial;accepted=True;break
        if not accepted:break
    raise RuntimeError('Initial axial moment balance did not converge')


def run(args):
    report=json.loads((args.seed/'report.json').read_text(encoding='utf-8'))
    payload=args.seed/'trajectory.npz'
    assert report['success'],'Only a completed seed is accepted'
    with np.load(payload) as d:
        q0=d['q'][-1];k0=d['k'][-1];origin=float(d['t'][-1])
    pars=report['parameters']
    beta,zeta,xi,chi,eta=m.parameters(lam=pars['slenderness'],eta=pars['eta'],
        xi=pars['xi'],chi=pars['chi'],nu=pars['nu'])
    q,k=remesh(q0,k0,args.segments)
    q,twist=consistent_twist(q,k,beta)
    nodes=m.geometry(q)[1];kp=k.copy();npv=nodes.copy()
    predictor=np.zeros((args.segments,3))
    steps=round(args.duration/args.dt)
    assert abs(steps*args.dt-args.duration)<1e-10
    stride=max(1,round(.01/args.dt))
    times=[origin];qs=[q.copy()];ks=[k.copy()];tips=[nodes[-1].copy()]
    worst=0.;hinge_max=0.;max_its=0;total_its=0
    start=time.perf_counter();success=True;completed_steps=0;failure=None
    print(json.dumps(dict(starting=True,segments=args.segments,dt=args.dt,
        time_origin=origin,requested_duration=args.duration,twist_projection=twist)),flush=True)
    for step in range(1,steps+1):
        order=1 if step==1 else 2
        qn,kn,inc,its,ok,res,hinge=m.implicit_step(q,k,nodes,kp,npv,beta,zeta,xi,chi,eta,args.dt,order,predictor,args.tolerance)
        worst=max(worst,res);hinge_max=max(hinge_max,hinge);max_its=max(max_its,its);total_its+=its
        if not ok:
            success=False;failure=dict(attempted_step=step,time=origin+step*args.dt,residual=res)
            print(json.dumps(failure),flush=True);break
        kp,npv=k.copy(),nodes.copy()
        q,k=qn,kn;nodes=m.geometry(q)[1];predictor=inc;completed_steps=step
        if step%stride==0 or step==steps:
            times.append(origin+step*args.dt);qs.append(q.copy());ks.append(k.copy());tips.append(nodes[-1].copy())
        if step%1000==0:
            print(json.dumps(dict(time=origin+step*args.dt,elapsed=time.perf_counter()-start,
                residual=res,tip=nodes[-1].tolist())),flush=True)
    arrays=dict(t=np.array(times),q=np.array(qs),k=np.array(ks),tip=np.array(tips))
    np.savez_compressed(args.output/'trajectory.npz',**arrays)
    np.savez_compressed(args.output/'checkpoint.npz',q=q,k=k,nodes=nodes,k_previous=kp,
        nodes_previous=npv,predictor=predictor,time=origin+completed_steps*args.dt,
        dt=args.dt,completed_steps=completed_steps)
    result=dict(parameters=pars,segments=args.segments,dt=args.dt,time_origin=origin,
        requested_duration=args.duration,completed_time=float(times[-1]),
        final_accepted_time=origin+completed_steps*args.dt,success=success,failure=failure,
        plane='spatial',order=2,startup='One backward-Euler step for a new endpoint-initialized IVP',
        newton_tolerance=args.tolerance,max_scaled_residual=worst,max_hinge_angle=hinge_max,
        max_newton_iterations=max_its,mean_newton_iterations=total_its/step,
        initial_twist=twist,elapsed_seconds=time.perf_counter()-start,
        seed=dict(external_path=str(payload),sha256=digest(payload),report_sha256=digest(args.seed/'report.json'),
            segments=len(q0),time=origin),
        remeshing='None' if args.segments==len(q0) else 'Geodesic orientation interpolation/extrapolation and linear curvature interpolation at hinge midpoints, followed by axial moment projection at fixed centerline',
        code_sha256=digest(Path(__file__)),solver_sha256=digest(Path(m.__file__)),
        python=platform.python_version(),numpy=np.__version__,
        scope='A new initial-value run from a saved endpoint. It is not a bitwise continuation of the earlier BDF history. Spatially refined initial states are explicitly interpolated and twist-projected.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--segments',type=int,required=True)
    ap.add_argument('--dt',type=float,required=True)
    ap.add_argument('--duration',type=float,required=True)
    ap.add_argument('--tolerance',type=float,default=1e-8)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    run(args)
