"""Numerical shooting of the regular competitor found by fresh basin starts."""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root as solve_root
from scipy.spatial import cKDTree

base=Path(__file__).parent;d=np.load(base/'icosa-v1/input.npz');b=np.load(base/'icosa-basins/data.npz');out=base/'icosa-cycle';out.mkdir(exist_ok=False)
ident=183;c=d['centers'][ident];u=d['axes'][ident];s=d['sigma'][ident];g=d['gamma'][ident]
skew=np.zeros((60,3,3));skew[:,0,1]=-u[:,2];skew[:,0,2]=u[:,1];skew[:,1,0]=u[:,2];skew[:,1,2]=-u[:,0];skew[:,2,0]=-u[:,1];skew[:,2,1]=u[:,0]
def field(t,x):
    r=x[:3]-c;v=np.cross(u,r);w=np.exp(-.5*(r*r).sum(1)/s**2);f=-g*x[:3]+(w[:,None]*v).sum(0)
    if len(x)==3:return f
    j=-g*np.eye(3)+(w[:,None,None]*(skew-v[:,:,None]*r[:,None,:]/s**2)).sum(0)
    return np.r_[f,(j@x[3:].reshape(3,3)).ravel()]
xref=b['end'][7];normal=field(0,xref);normal/=np.linalg.norm(normal)
def event(t,x):return (x-xref)@normal
event.direction=1
sol=solve_ivp(field,[0,100],xref,method='DOP853',rtol=1e-11,atol=1e-13,max_step=.05,events=event)
times=sol.t_events[0];states=sol.y_events[0];mask=(times>1)&(np.linalg.norm(states-xref,axis=1)<1e-5)
assert mask.any();guess=float(times[mask][0]);cached={}
def evaluate(v):
    if 'v' in cached and np.array_equal(v,cached['v']):return cached['F'],cached['J']
    assert v[3]>0
    q=solve_ivp(field,[0,v[3]],np.r_[v[:3],np.eye(3).ravel()],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.05)
    end=q.y[:3,-1];M=q.y[3:,-1].reshape(3,3);F=np.r_[end-v[:3],(v[:3]-xref)@normal]
    J=np.block([[M-np.eye(3),field(0,end)[:,None]],[normal[None,:],np.zeros((1,1))]])
    cached.update(v=v.copy(),F=F,J=J,M=M);return F,J
begin=time.time();result=solve_root(lambda v:evaluate(v)[0],np.r_[xref,guess],jac=lambda v:evaluate(v)[1],tol=1e-10)
F,J=evaluate(result.x);assert np.linalg.norm(F)<1e-10
ev=np.linalg.eigvals(cached['M']);nontrivial=np.delete(ev,np.argmin(np.abs(ev-1)));assert np.all(np.abs(nontrivial)<1)
period=float(result.x[3]);orbit=solve_ivp(field,[0,period],result.x[:3],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.02,t_eval=np.linspace(0,period,8000,endpoint=False)).y.T
others=[]
for row in [7,9,10,12,14]:
    tail=b['tail'][row,-2000:];errors=[float(np.sqrt(np.mean(cKDTree(orbit@a.T).query(tail)[0]**2))) for a in d['group']];others.append({'row':row,'best_rotation_rms':min(errors)})
report={'id':183,'basin_row':7,'initial_period_guess':guess,'period':period,'point':result.x[:3].tolist(),'shooting_residual':float(np.linalg.norm(F)),'multipliers':[[float(v.real),float(v.imag)] for v in ev],'nontrivial_multiplier_moduli':np.abs(nontrivial).tolist(),'log_abs_determinant_over_period':float(np.log(abs(np.linalg.det(cached['M'])))/period),'trace_target':float(-3*g),'other_regular_starts':others,'seconds':time.time()-begin,'limits':'Numerical periodic-orbit shooting and Floquet stability, not an interval existence proof. Least period is not certified. Symmetry matches compare finite samples.'}
np.savez_compressed(out/'orbit.npz',points=orbit,monodromy=cached['M']);(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
