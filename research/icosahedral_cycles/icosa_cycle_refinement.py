"""Full-cycle multiple shooting and independent RK4 convergence for two candidates."""
from pathlib import Path
import sys
from cycle_paths import ROOT,RESEARCH,INPUT,STRUCTURE,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from numba import njit
from tetra_spectra import field as numba_field

base=OUTPUT;d=np.load(INPUT/'icosa-v1/input.npz');previous=json.loads((base/'icosa-relative-cycles/report.json').read_text());out=base/'icosa-cycle-refinement';out.mkdir(exist_ok=False)
c=d['centers'][177];u=d['axes'][177];s=d['sigma'][177];gamma=d['gamma'][177];group=d['group']
skew=np.zeros((60,3,3));skew[:,0,1]=-u[:,2];skew[:,0,2]=u[:,1];skew[:,1,0]=u[:,2];skew[:,1,2]=-u[:,0];skew[:,2,0]=-u[:,1];skew[:,2,1]=u[:,0]
def field(t,x):
    r=x[:3]-c;v=np.cross(u,r);w=np.exp(-.5*(r*r).sum(1)/s**2);f=-gamma*x[:3]+(w[:,None]*v).sum(0)
    if len(x)==3:return f
    j=-gamma*np.eye(3)+(w[:,None,None]*(skew-v[:,:,None]*r[:,None,:]/s**2)).sum(0)
    return np.r_[f,(j@x[3:].reshape(3,3)).ravel()]
@njit
def rk4(x,T,nominal_h):
    n=int(np.ceil(T/nominal_h));h=T/n;z=x.copy()
    for _ in range(n):
        a=numba_field(z,c,u,s,gamma)[0];b=numba_field(z+.5*h*a,c,u,s,gamma)[0];e=numba_field(z+.5*h*b,c,u,s,gamma)[0];f=numba_field(z+h*e,c,u,s,gamma)[0]
        z+=h/6*(a+2*b+2*e+f)
    return z,n,h
rows=[]
for index in [0,20]:
    prev=previous['rows'][index];m=prev['rotation_order'];G=group[prev['guess']['rotation']];xref=np.array(prev['point']);normal=field(0,xref);normal/=np.linalg.norm(normal);cached={}
    initial=np.r_[np.concatenate([np.linalg.matrix_power(G,k)@xref for k in range(m)]),prev['relative_period']]
    def evaluate(v):
        if 'v' in cached and np.array_equal(v,cached['v']):return cached['F'],cached['J']
        nodes=v[:-1].reshape(m,3);tau=v[-1];F=np.empty(3*m+1);J=np.zeros((3*m+1,3*m+1));Ms=[]
        for k in range(m):
            sol=solve_ivp(field,[0,tau],np.r_[nodes[k],np.eye(3).ravel()],method='DOP853',rtol=3e-14,atol=3e-16,max_step=.025)
            end=sol.y[:3,-1];M=sol.y[3:,-1].reshape(3,3);Ms.append(M)
            F[3*k:3*k+3]=end-nodes[(k+1)%m];J[3*k:3*k+3,3*k:3*k+3]+=M;J[3*k:3*k+3,3*((k+1)%m):3*((k+1)%m)+3]-=np.eye(3);J[3*k:3*k+3,-1]=field(0,end)
        F[-1]=(nodes[0]-xref)@normal;J[-1,:3]=normal
        cached.update(v=v.copy(),F=F,J=J,Ms=Ms);return F,J
    begin=time.time();fit=root(lambda v:evaluate(v)[0],initial,jac=lambda v:evaluate(v)[1],tol=1e-10);F,J=evaluate(fit.x);assert np.linalg.norm(F)<1e-10
    nodes=fit.x[:-1].reshape(m,3);tau=float(fit.x[-1]);T=m*tau;point=nodes[0]
    M=np.eye(3)
    for a in cached['Ms']:M=a@M
    ev=np.linalg.eigvals(M);unstable=float(ev[np.argmax(abs(ev))].real);stable_inferred=float(np.exp(-3*gamma*T)/abs(unstable));checks=[]
    sol=solve_ivp(field,[0,T],point,method='DOP853',rtol=3e-14,atol=3e-16,max_step=.01,t_eval=np.linspace(0,T,40001))
    independent_error=float(np.linalg.norm(sol.y[:,-1]-point))
    for h in [.004,.002,.001,.0005]:
        end,n,actual_h=rk4(point,T,h);checks.append({'nominal_h':h,'actual_h':actual_h,'steps':n,'endpoint':end.tolist(),'closure_error':float(np.linalg.norm(end-point))})
    np.savez_compressed(out/f'cycle-{index}.npz',points=sol.y[:,:-1].T,nodes=nodes,segment_monodromies=np.array(cached['Ms']),full_monodromy=M)
    row={'index':index,'segments':m,'point':point.tolist(),'segment_time':tau,'period':T,'shooting_residual':float(np.linalg.norm(F)),'adaptive_full_closure_error':independent_error,'multipliers_direct':[[float(a.real),float(a.imag)] for a in ev],'unstable_multiplier':unstable,'stable_modulus_inferred_from_divergence':stable_inferred,'rotation':prev['guess']['rotation'],'segment_symmetry_residual':float(max(np.linalg.norm(nodes[k]-np.linalg.matrix_power(G,k)@nodes[0]) for k in range(m))),'rk4':checks,'seconds':time.time()-begin}
    rows.append(row);(out/'report.json').write_text(json.dumps({'rows':rows,'limits':'Full periodic-orbit multiple shooting uses the actual rounded coefficient field. Direct tiny Floquet multipliers can lose precision; the stable magnitude is also inferred from exact divergence and the neutral multiplier. Existence and least period are not interval certified.'},indent=2));print(json.dumps(row),flush=True)
