"""Straight and chaotic geodesics of a three-cosine surface in a flat phase torus."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp

TAU = 2 * np.pi


@njit(cache=True)
def field(state, n, k):
    x,y,z,vx,vy,vz = state
    w = x-y+n*z
    u = k*np.sin(w)
    gx,gy,gz = -np.sin(x)-u,-np.sin(y)+u,-n*u
    gg = gx*gx+gy*gy+gz*gz
    f = np.cos(x)+np.cos(y)+k*np.cos(w)
    fd = gx*vx+gy*vy+gz*vz
    h = -np.cos(x)*vx*vx-np.cos(y)*vy*vy-k*np.cos(w)*(vx-vy+n*vz)**2
    normal = -(h+3*fd+2*f)/gg
    tx,ty,tz = vx-fd*gx/gg,vy-fd*gy/gg,vz-fd*gz/gg
    speed = 1-(tx*tx+ty*ty+tz*tz)
    return np.array([vx,vy,vz,normal*gx+speed*tx,normal*gy+speed*ty,normal*gz+speed*tz])


@njit(cache=True)
def step(s, n, k, dt):
    a = field(s,n,k)
    b = field(s+dt*a/2,n,k)
    c = field(s+dt*b/2,n,k)
    d = field(s+dt*c,n,k)
    return s+(dt/6)*(a+2*b+2*c+d)


@njit(cache=True)
def sample(s,n,k,dt,stride,count):
    points=np.empty((count,6))
    for i in range(count):
        for _ in range(stride):
            s=step(s,n,k,dt)
        points[i]=s
    return points


def start(n,k,angle):
    q=np.array([np.pi/2,np.pi/2,np.pi/(2*n)])
    a=np.array([1.,-1.,-2/n])
    a/=np.linalg.norm(a)
    normal=np.array([-1-k,-1+k,-n*k])
    b=np.cross(normal,a)
    b/=np.linalg.norm(b)
    return np.r_[q,np.cos(angle)*a+np.sin(angle)*b]


def independent(t,s,n,k):
    L=np.array([[1.,0.,0.],[0.,1.,0.],[1.,-1.,float(n)]])
    amp=np.array([1.,1.,k])
    q,v=s[:3],s[3:]
    p=L@q
    f=amp@np.cos(p)
    g=-L.T@(amp*np.sin(p))
    quadratic=-(amp*np.cos(p))@((L@v)**2)
    gd=g@v
    tangent=v-g*(gd/(g@g))
    return np.r_[v,-g*((quadratic+3*gd+2*f)/(g@g))+(1-tangent@tangent)*tangent]


def floquet(n,k):
    speed=np.sqrt(2.+4./n**2)
    # Curvature repeats after pi in x. Use this short segment to avoid a poorly conditioned large monodromy.
    duration=np.pi*speed
    def variational(t,y):
        x=np.pi/2+t/speed
        minus_K=(n*k*np.cos(x)/(2*np.sin(x)**2+(n*n+2)*k*k))**2
        return [y[2],y[3],minus_K*y[0],minus_K*y[1]]
    sol=solve_ivp(variational,(0,duration),[1.,0.,0.,1.],method="DOP853",rtol=2e-12,atol=2e-14)
    assert sol.success
    m=sol.y[:,-1].reshape(2,2)
    assert abs(np.linalg.det(m)-1)<1e-9
    exponent=np.arccosh(np.trace(m)/2)/duration
    cycles=n//np.gcd(n,2)
    period=TAU*cycles*speed
    return dict(n=n,k=k,genus=2*n+1,ray_period=float(period),ray_winding=[int(cycles),-int(cycles),int(-2*cycles/n)],
                transverse_rate=float(exponent),transverse_multiplier=float(np.exp(period*exponent)),short_monodromy=m.tolist())


def errors(s,n,k):
    q,v=s[:,:3],s[:,3:]
    w=q[:,0]-q[:,1]+n*q[:,2]
    f=np.cos(q[:,0])+np.cos(q[:,1])+k*np.cos(w)
    fd=-np.sin(q[:,0])*v[:,0]-np.sin(q[:,1])*v[:,1]-k*np.sin(w)*(v[:,0]-v[:,1]+n*v[:,2])
    return dict(surface=float(abs(f).max()),tangency=float(abs(fd).max()),speed=float(abs(np.linalg.norm(v,axis=1)-1).max()))


def run(out,report):
    out.mkdir(parents=True,exist_ok=True)
    report.parent.mkdir(parents=True,exist_ok=True)
    n,k=7,0.375
    rows=[floquet(m,c) for m in (1,2,3,4,5,7,9,13,21) for c in (.125,.25,.375,.5,.75,1.,1.5,1.875)]
    ray=floquet(n,k)
    perturb=[]
    for angle in (1e-4,1e-6,1e-8):
        s=start(n,k,angle)
        def leave(t,state):
            delta=state[:3]-s[:3]-t*start(n,k,0.)[3:]
            return np.linalg.norm(delta)-.2
        leave.terminal=True
        leave.direction=1
        estimates=[]
        for method,rtol,atol in (("DOP853",2e-12,2e-14),("RK45",2e-11,2e-13)):
            sol=solve_ivp(independent,(0,180),s,args=(n,k),method=method,rtol=rtol,atol=atol,
                          events=leave_wrapped(leave),max_step=.15)
            assert sol.success and len(sol.t_events[0])
            estimates.append(sol.t_events[0][0])
        perturb.append(dict(angle=angle,dop853=estimates[0],rk45=estimates[1]))
    s=start(n,k,.321)
    controls={}
    for angle in (.321,-.51,1.08):
        y=start(n,k,angle)
        tr=sample(y,n,k,.0125,4,160)
        sol=solve_ivp(independent,(0,8),y,args=(n,k),method="DOP853",t_eval=np.arange(1,161)*.05,rtol=2e-12,atol=2e-14)
        controls[str(angle)]=dict(short_state_error=float(np.max(abs(sol.y.T-tr))),**errors(tr,n,k))
    h=.0125
    raw=sample(s,n,k,h,8,280000)
    q=(raw[:,:3]+np.pi)%TAU-np.pi
    vol,_=np.histogramdd(q,bins=36,range=[[-np.pi,np.pi]]*3)
    # The global sheet index is defined with x,y in the reference cell.
    p=q[:,0]-q[:,1]+n*q[:,2]
    sheets=np.floor(np.mod(p,TAU*n)/np.pi).astype(int)
    checks=dict(errors=errors(raw,n,k),occupied_bins=int(np.count_nonzero(vol)),sheet_counts=np.bincount(sheets,minlength=2*n).tolist())
    file=out/'lattice-cloud.npz'
    np.savez_compressed(file,points=q)
    # Display positions from one and the same trajectory; the exact ray uses the analytic formula.
    meta=dict(selected=ray,sweep=rows,perturbations=perturb,controller_checks=controls,long_run=checks,
              samples=len(raw),step=h,stride=8,initial=s.tolist(),angle=.321,file=file.name,
              sha256=hashlib.sha256(file.read_bytes()).hexdigest())
    report.write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(ray=ray,perturbations=perturb,metrics=checks)),flush=True)


def leave_wrapped(function):
    def with_parameters(t,s,n,k):
        return function(t,s)
    with_parameters.terminal=True
    with_parameters.direction=1
    return with_parameters


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    run(args.out,args.report)
