"""Overdamped spatial morphoelastic rod: geometric segment approximation.

Constant-step implicit BDF1/BDF2 on position and spontaneous curvature.
Directors use unit quaternions. Axial torque is algebraic: no rotational drag
or inertia is introduced. Large histories are written only to the chosen path.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np
from numba import njit


@njit(cache=True)
def mul(a,b):
    return np.array([a[0]*b[0]-np.dot(a[1:],b[1:]),
        a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
        a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
        a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]])


@njit(cache=True)
def expq(v):
    angle=np.sqrt(np.dot(v,v))
    f=np.sin(angle/2)/angle if angle>1e-8 else 0.5-angle*angle/48
    return np.array([np.cos(angle/2),f*v[0],f*v[1],f*v[2]])


@njit(cache=True)
def logq(q):
    q=q.copy()
    if q[0]<0: q=-q
    size=np.sqrt(np.dot(q[1:],q[1:]))
    f=2*np.arctan2(size,q[0])/size if size>1e-10 else 2/q[0]
    return f*q[1:]


@njit(cache=True)
def matrix(q):
    w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],
                     [2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],
                     [2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]])


@njit(cache=True)
def rotate_increment(q,inc):
    out=np.empty_like(q)
    for i in range(len(q)):
        out[i]=mul(q[i],expq(inc[i]))
        out[i]/=np.sqrt(np.dot(out[i],out[i]))
    return out


@njit(cache=True)
def geometry(q):
    n=len(q)
    frames=np.empty((n,3,3))
    nodes=np.zeros((n+1,3))
    for i in range(n):
        frames[i]=matrix(q[i])
        nodes[i+1]=nodes[i]+frames[i,:,2]/n
    return frames,nodes


@njit(cache=True)
def strains_midframes(q):
    n=len(q)
    vectors=np.empty((n,3))
    mids=np.empty((n,3,3))
    prev=np.array([1.0,0.0,0.0,0.0])
    for i in range(n):
        conj=prev.copy();conj[1:]*=-1
        vectors[i]=logq(mul(conj,q[i]))
        mids[i]=matrix(mul(prev,expq(vectors[i]/2)))
        prev=q[i]
    return vectors,mids


@njit(cache=True)
def inv_left_jac_action(v,w,sign):
    angle=np.sqrt(np.dot(v,v))
    coef=(1-0.5*angle/np.tan(angle/2))/(angle*angle) if angle>1e-4 else 1/12+angle*angle/720
    cross=np.cross(v,w)
    return w+sign*0.5*cross+coef*np.cross(v,cross)


@njit(cache=True)
def energy_gradient(q,kappa,beta):
    n=len(q)
    vectors,mids=strains_midframes(q)
    grad=np.zeros((n,3))
    energy=0.0
    max_hinge=0.0
    for i in range(n):
        length=(0.5 if i==0 else 1.0)/n
        u=vectors[i]/length
        strain=u-np.array([kappa[i,0],kappa[i,1],0.0])
        moment=beta*strain
        energy+=0.5*length*np.dot(strain,moment)
        grad[i]+=inv_left_jac_action(vectors[i],moment,-1)
        if i>0: grad[i-1]-=inv_left_jac_action(vectors[i],moment,1)
        max_hinge=max(max_hinge,np.sqrt(np.dot(vectors[i],vectors[i])))
    return energy,grad,max_hinge


@njit(cache=True)
def target_curvature(q,chi,eta):
    vectors,mids=strains_midframes(q)
    out=np.empty((len(q),2))
    out[:,0]=-chi*mids[:,2,1]
    out[:,1]=chi*eta**-3*mids[:,2,0]
    return out


@njit(cache=True)
def drag_gradient(q,nodes_old,nodes_prev,zeta,xi,dt,order):
    frames,nodes=geometry(q)
    n=len(q);h=1.0/n
    if order==1:
        velocities=(nodes-nodes_old)/dt
    else:
        velocities=(1.5*nodes-2*nodes_old+0.5*nodes_prev)/dt
    grad=np.zeros((n,3));tail=np.zeros(3)
    for i in range(n-1,-1,-1):
        r=frames[i]
        drag=(r*zeta)@r.T*xi
        v0,v1=velocities[i],velocities[i+1]
        weighted=drag@(v0/6+v1/3)*h
        force=tail+weighted
        body=r.T@force
        grad[i]=h*np.array([-body[1],body[0],0.0])
        tail+=h*drag@((v0+v1)/2)
    return grad


@njit(cache=True)
def skew(v):
    return np.array([[0.0,-v[2],v[1]],[v[2],0.0,-v[0]],[-v[1],v[0],0.0]])


@njit(cache=True)
def drag_jacobian(q,nodes_old,nodes_prev,zeta,xi,dt,order):
    # The translational mass is symmetric. Precontract each tail with its
    # right tangent action, then fill only its two active rotational axes.
    # Drag-coefficient variations depend on the varied segment, so their
    # force vectors are formed once before propagation to upstream segments.
    frames,nodes=geometry(q)
    n=len(q);h=1.0/n;factor=(1.0 if order==1 else 1.5)/dt
    velocity=(nodes-nodes_old)/dt if order==1 else (1.5*nodes-2*nodes_old+0.5*nodes_prev)/dt
    ez=np.array([0.0,0.0,1.0]);s=-skew(ez)
    drags=np.empty((n,3,3));actions=np.empty((n,3,3));tails=np.zeros((n+1,3,3));forces=np.zeros((n+1,3))
    for i in range(n-1,-1,-1):
        drags[i]=(frames[i]*zeta)@frames[i].T*xi
        actions[i]=frames[i]@s
        tails[i]=tails[i+1]+drags[i]
        forces[i]=forces[i+1]+h*drags[i]@((velocity[i]+velocity[i+1])/2)
    jac=np.zeros((3*n,3*n))
    for j in range(n):
        right=(tails[j+1]+drags[j]/2)@actions[j]
        diagonal=(tails[j+1]+drags[j]/3)@actions[j]
        for i in range(j+1):
            block=diagonal if i==j else right
            for a in range(2):
                for b in range(2):
                    value=0.0
                    for component in range(3):
                        value+=actions[i,component,a]*block[component,b]
                    value*=factor*h**3
                    jac[3*i+a,3*j+b]=value
                    jac[3*j+b,3*i+a]=value
        weighted=velocity[j]/6+velocity[j+1]/3
        body=frames[j].T@(forces[j+1]+h*drags[j]@weighted)
        jac[3*j,3*j]-=h*body[2]
        jac[3*j+1,3*j+1]-=h*body[2]
        jac[3*j,3*j+2]+=h*body[0]
        jac[3*j+1,3*j+2]+=h*body[1]
        for axis in range(3):
            basis=np.zeros(3);basis[axis]=1.0
            cross=skew(basis)
            local=cross*np.expand_dims(zeta,0)-np.expand_dims(zeta,1)*cross
            delta=xi*frames[j]@local@frames[j].T
            avg_force=delta@((velocity[j]+velocity[j+1])/2)
            weighted_force=delta@weighted
            for i in range(j+1):
                force=weighted_force if i==j else avg_force
                for a in range(2):
                    value=0.0
                    for component in range(3):
                        value+=actions[i,component,a]*force[component]
                    jac[3*i+a,3*j+axis]+=h*h*value
    return jac


@njit(cache=True)
def hinge_eval(prev,current,k_old,k_prev,beta,chi,eta,dt,order,length):
    conj=prev.copy();conj[1:]*=-1
    vector=logq(mul(conj,current))
    middle=matrix(mul(prev,expq(vector/2)))
    drive=np.array([-chi*middle[2,1],chi*eta**-3*middle[2,0]])
    relax=np.array([1.0,eta**-2])
    k=(k_old+dt*drive)/(1+dt*relax) if order==1 else (2*k_old-0.5*k_prev+dt*drive)/(1.5+dt*relax)
    moment=beta*(vector/length-np.array([k[0],k[1],0.0]))
    return -inv_left_jac_action(vector,moment,1),inv_left_jac_action(vector,moment,-1)


@njit(cache=True)
def mechanical_jacobian(q,kold,kprev,nodes_old,nodes_prev,beta,zeta,xi,chi,eta,dt,order):
    jac=drag_jacobian(q,nodes_old,nodes_prev,zeta,xi,dt,order)
    n=len(q);identity=np.array([1.0,0.0,0.0,0.0])
    eps=1e-6
    for j in range(n):
        for axis in range(3):
            direction=np.zeros(3);direction[axis]=eps
            plus=mul(q[j],expq(direction));minus=mul(q[j],expq(-direction))
            prev=identity if j==0 else q[j-1]
            length=(0.5 if j==0 else 1.0)/n
            lp,rp=hinge_eval(prev,plus,kold[j],kprev[j],beta,chi,eta,dt,order,length)
            lm,rm=hinge_eval(prev,minus,kold[j],kprev[j],beta,chi,eta,dt,order,length)
            if j>0: jac[3*(j-1):3*j,3*j+axis]+=(lp-lm)/(2*eps)
            jac[3*j:3*j+3,3*j+axis]+=(rp-rm)/(2*eps)
            if j<n-1:
                lp,rp=hinge_eval(plus,q[j+1],kold[j+1],kprev[j+1],beta,chi,eta,dt,order,1.0/n)
                lm,rm=hinge_eval(minus,q[j+1],kold[j+1],kprev[j+1],beta,chi,eta,dt,order,1.0/n)
                jac[3*j:3*j+3,3*j+axis]+=(lp-lm)/(2*eps)
                jac[3*(j+1):3*(j+2),3*j+axis]+=(rp-rm)/(2*eps)
    return jac/beta[0]


@njit(cache=True)
def residual_q(q,k_old,nodes_old,k_prev,nodes_prev,beta,zeta,xi,chi,eta,dt,order):
    drive=target_curvature(q,chi,eta)
    relax=np.array([1.0,eta**-2])
    if order==1:
        k=(k_old+dt*drive)/(1+dt*relax)
    else:
        k=(2*k_old-0.5*k_prev+dt*drive)/(1.5+dt*relax)
    energy,elastic,max_hinge=energy_gradient(q,k,beta)
    drag=drag_gradient(q,nodes_old,nodes_prev,zeta,xi,dt,order)
    return (elastic+drag)/beta[0],q,k,max_hinge


@njit(cache=True)
def residual(inc,q_old,k_old,nodes_old,k_prev,nodes_prev,beta,zeta,xi,chi,eta,dt,order):
    q=rotate_increment(q_old,inc)
    return residual_q(q,k_old,nodes_old,k_prev,nodes_prev,beta,zeta,xi,chi,eta,dt,order)


@njit(cache=True)
def implicit_step(q_old,k_old,nodes_old,k_prev,nodes_prev,beta,zeta,xi,chi,eta,dt,order,predictor,tol):
    q=rotate_increment(q_old,predictor);n=len(q_old);inc=predictor.copy()
    iterations=0
    max_residual=1e300
    for it in range(15):
        r,q,k,max_hinge=residual_q(q,k_old,nodes_old,k_prev,nodes_prev,beta,zeta,xi,chi,eta,dt,order)
        norm=np.max(np.abs(r))
        max_residual=norm
        if norm<tol:
            for i in range(n):
                conj=q_old[i].copy();conj[1:]*=-1
                inc[i]=logq(mul(conj,q[i]))
            return q,k,inc,it,True,norm,max_hinge
        jac=mechanical_jacobian(q,k_old,k_prev,nodes_old,nodes_prev,beta,zeta,xi,chi,eta,dt,order)
        change=np.linalg.solve(jac,-r.ravel()).reshape((n,3))
        factor=1.0
        accepted=False
        for ls in range(12):
            trial=rotate_increment(q,factor*change)
            rn=residual_q(trial,k_old,nodes_old,k_prev,nodes_prev,beta,zeta,xi,chi,eta,dt,order)[0]
            if np.max(np.abs(rn)) < norm:
                q=trial;accepted=True;break
            factor*=0.5
        iterations=it+1
        if not accepted: break
    r,q,k,max_hinge=residual_q(q,k_old,nodes_old,k_prev,nodes_prev,beta,zeta,xi,chi,eta,dt,order)
    return q,k,inc,iterations,False,float(np.max(np.abs(r))),max_hinge


def parameters(lam=50,eta=0.9,xi=5e-6,chi=100,nu=0.5):
    log=np.log(2*lam/(1+eta))
    beta=np.array([np.pi*eta/(4*lam**4),np.pi*eta**3/(4*lam**4),
                   np.pi*eta**3/(2*lam**4*(1+eta**2)*(1+nu))])
    zeta=np.array([4*np.pi/(log+eta/(1+eta)),4*np.pi/(log+1/(1+eta)),2*np.pi/(log-0.5)])
    return beta,zeta,xi,chi,eta


def consistent_initial(n,beta,plane):
    s=(np.arange(n)+0.5)/n
    q=np.zeros((n,4))
    for i in range(n):
        x=0.01*s[i]*(2-s[i])
        y=0.007*s[i]**2*(3-2*s[i])
        q[i]=expq(np.array([-x,0.0,0.0]) if plane=='yz' else
                  np.array([0.0,x,0.0]) if plane=='xz' else np.array([-x,y,0.0]))
    k=np.zeros((n,2));spin=np.zeros((n,3))
    for it in range(10):
        current=rotate_increment(q,spin)
        res=energy_gradient(current,k,beta)[1][:,2]/beta[0]
        if max(abs(res))<1e-11:
            return current,k,float(max(abs(res)))
        jac=np.empty((n,n))
        for j in range(n):
            shifted=spin.copy();shifted[j,2]+=1e-6
            plus=energy_gradient(rotate_increment(q,shifted),k,beta)[1][:,2]/beta[0]
            shifted[j,2]-=2e-6
            minus=energy_gradient(rotate_increment(q,shifted),k,beta)[1][:,2]/beta[0]
            jac[:,j]=(plus-minus)/2e-6
        spin[:,2]+=np.linalg.solve(jac,-res)
    raise RuntimeError('Initial twist equilibrium did not converge')


def run(n=16,dt=0.002,end=30,lam=50,eta=0.9,xi=5e-6,chi=100,nu=0.5,plane='spatial',order=2,tol=1e-8,output=None):
    beta,zeta,xi,chi,eta=parameters(lam,eta,xi,chi,nu)
    q,k,initial_torque=consistent_initial(n,beta,plane)
    nodes=geometry(q)[1]
    kprev,nprev=k.copy(),nodes.copy()
    predictor=np.zeros((n,3))
    stride=max(1,round(0.01/dt));steps=round(end/dt)
    times=[0.0];qs=[q.copy()];ks=[k.copy()];tips=[nodes[-1].copy()]
    worst=0.0;hinge_max=0.0;max_iterations=0;total_iterations=0
    start=time.perf_counter();success=True
    for step in range(1,steps+1):
        actual_order=1 if step==1 else order
        qnew,knew,inc,its,ok,res,hinge=implicit_step(q,k,nodes,kprev,nprev,beta,zeta,xi,chi,eta,dt,actual_order,predictor,tol)
        worst=max(worst,res);hinge_max=max(hinge_max,hinge)
        max_iterations=max(max_iterations,its);total_iterations+=its
        if not ok:
            success=False
            print(json.dumps(dict(failed_step=step,time=step*dt,residual=res)),flush=True)
            break
        kprev,nprev=k.copy(),nodes.copy()
        q,k=qnew,knew;nodes=geometry(q)[1]
        predictor=inc
        if step%stride==0:
            times.append(step*dt);qs.append(q.copy());ks.append(k.copy());tips.append(nodes[-1].copy())
        if step%1000==0:
            print(json.dumps(dict(time=step*dt,elapsed=time.perf_counter()-start,residual=res,tip=nodes[-1].tolist())),flush=True)
    times=np.array(times);qs=np.array(qs);ks=np.array(ks);tips=np.array(tips)
    late=times>=times[-1]/2
    centered=tips[late]-tips[late].mean(axis=0)
    singular=np.linalg.svd(centered,compute_uv=False)
    report=dict(parameters=dict(slenderness=lam,eta=eta,xi=xi,chi=chi,nu=nu),segments=n,dt=dt,
                requested_duration=end,completed_time=float(times[-1]),plane=plane,order=order,newton_tolerance=tol,
                success=success,max_scaled_residual=worst,max_hinge_angle=hinge_max,max_newton_iterations=max_iterations,
                mean_newton_iterations=total_iterations/step,initial_axial_torque_scaled=initial_torque,
                elapsed_seconds=time.perf_counter()-start,late_tip_singular_values=singular.tolist(),
                late_tip_ranges=np.c_[tips[late].min(axis=0),tips[late].max(axis=0)].tolist(),
                quaternion_norm_error=float(np.max(np.abs(np.sum(qs*qs,axis=-1)-1))),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=platform.python_version(),numpy=np.__version__)
    if output is not None:
        np.savez_compressed(output/'trajectory.npz',t=times,q=qs,k=ks,tip=tips)
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--segments',type=int,default=16)
    ap.add_argument('--dt',type=float,default=0.002)
    ap.add_argument('--duration',type=float,default=30)
    ap.add_argument('--slenderness',type=float,default=50)
    ap.add_argument('--eta',type=float,default=0.9)
    ap.add_argument('--xi',type=float,default=5e-6)
    ap.add_argument('--chi',type=float,default=100)
    ap.add_argument('--nu',type=float,default=0.5)
    ap.add_argument('--plane',choices=['yz','xz','spatial'],default='spatial')
    ap.add_argument('--order',type=int,choices=[1,2],default=2)
    ap.add_argument('--tolerance',type=float,default=1e-8)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    run(a.segments,a.dt,a.duration,a.slenderness,a.eta,a.xi,a.chi,a.nu,a.plane,a.order,a.tolerance,a.output)


if __name__=='__main__': main()
