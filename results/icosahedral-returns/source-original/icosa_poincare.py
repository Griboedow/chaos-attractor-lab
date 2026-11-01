"""Adaptive local section return with the event-time derivative correction."""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

base=Path(__file__).parent;d=np.load(base/'icosa-v1/input.npz');frame=np.load(base/'icosa-return-scout/data.npz')
c=d['centers'][177];u=d['axes'][177];s=d['sigma'][177];g=d['gamma'][177];x0=frame['point'];n=frame['normal'];B=frame['basis'];inv=frame['inverse_basis']
skew=np.zeros((60,3,3));skew[:,0,1]=-u[:,2];skew[:,0,2]=u[:,1];skew[:,1,0]=u[:,2];skew[:,1,2]=-u[:,0];skew[:,2,0]=-u[:,1];skew[:,2,1]=u[:,0]
def field(t,x):
    r=x[:3]-c;v=np.cross(u,r);w=np.exp(-.5*(r*r).sum(1)/s**2);f=-g*x[:3]+(w[:,None]*v).sum(0)
    if len(x)==3:return f
    J=-g*np.eye(3)+(w[:,None,None]*(skew-v[:,:,None]*r[:,None,:]/s**2)).sum(0)
    return np.r_[f,(J@x[3:].reshape(3,2)).ravel()]
def section(t,x):return (x[:3]-x0)@n
section.direction=1
def evaluate(uv,rtol=2e-12,atol=2e-14,max_step=.04,horizon=25.):
    start=x0+B@uv;sol=solve_ivp(field,[0,horizon],np.r_[start,B.ravel()],method='DOP853',rtol=rtol,atol=atol,max_step=max_step,events=section,dense_output=True)
    assert sol.success,sol.message
    for t,x in zip(sol.t_events[0],sol.y_events[0]):
        q=inv@(x[:3]-x0);f=field(0,x[:3]);speed=n@f
        if t>1. and abs(q[0])<.15 and abs(q[1])<.01 and speed>.05:
            J=inv@(np.eye(3)-np.outer(f,n)/speed)@x[3:].reshape(3,2)
            return q,float(t),J,float(speed)
    raise RuntimeError(f'No local return by {horizon}: {uv}')
