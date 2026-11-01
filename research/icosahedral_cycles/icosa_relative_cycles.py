"""Adaptive shooting of relative periodic-orbit candidates, retaining failures."""
from pathlib import Path
import sys
from cycle_paths import ROOT,RESEARCH,INPUT,STRUCTURE,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

base=OUTPUT;d=np.load(INPUT/'icosa-v1/input.npz');rec=json.loads((base/'icosa-recurrences/candidates.json').read_text());out=base/'icosa-relative-cycles';out.mkdir(exist_ok=False)
ident=177;c=d['centers'][ident];u=d['axes'][ident];s=d['sigma'][ident];gamma=d['gamma'][ident];group=d['group']
skew=np.zeros((60,3,3));skew[:,0,1]=-u[:,2];skew[:,0,2]=u[:,1];skew[:,1,0]=u[:,2];skew[:,1,2]=-u[:,0];skew[:,2,0]=-u[:,1];skew[:,2,1]=u[:,0]
def field(t,x):
    r=x[:3]-c;v=np.cross(u,r);w=np.exp(-.5*(r*r).sum(1)/s**2);f=-gamma*x[:3]+(w[:,None]*v).sum(0)
    if len(x)==3:return f
    j=-gamma*np.eye(3)+(w[:,None,None]*(skew-v[:,:,None]*r[:,None,:]/s**2)).sum(0)
    return np.r_[f,(j@x[3:].reshape(3,3)).ravel()]
rows=[];orbits={};start_time=time.time()
for index,guess in enumerate(rec['selected'][:24]):
    G=group[guess['rotation']];xref=np.array(guess['start']);normal=field(0,xref);normal/=np.linalg.norm(normal);cache={}
    def evaluate(v):
        if 'v' in cache and np.array_equal(cache['v'],v):return cache['F'],cache['J']
        sol=solve_ivp(field,[0,v[3]],np.r_[v[:3],np.eye(3).ravel()],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.05)
        end=sol.y[:3,-1];M=sol.y[3:,-1].reshape(3,3);F=np.r_[end-G@v[:3],(v[:3]-xref)@normal]
        J=np.block([[M-G,field(0,end)[:,None]],[normal[None,:],np.zeros((1,1))]])
        cache.update(v=v.copy(),F=F,J=J,M=M);return F,J
    t=guess['time_guess'];begin=time.time()
    fit=least_squares(lambda v:evaluate(v)[0],np.r_[xref,t],jac=lambda v:evaluate(v)[1],bounds=(np.r_[xref-.3,max(.2,t*.6)],np.r_[xref+.3,t*1.4]),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=50)
    F,J=evaluate(fit.x);residual=float(np.linalg.norm(F));speed=float(np.linalg.norm(field(0,fit.x[:3])));accepted=residual<1e-8 and speed>1e-4
    row={'index':index,'guess':guess,'accepted_relative_root':accepted,'point':fit.x[:3].tolist(),'relative_period':float(fit.x[3]),'residual':residual,'speed':speed,'evaluations':fit.nfev,'seconds':time.time()-begin}
    if accepted:
        ev=np.linalg.eigvals(G.T@cache['M']);order=next(m for m in range(1,7) if np.linalg.norm(np.linalg.matrix_power(G,m)-np.eye(3))<1e-10)
        neutral=int(np.argmin(abs(ev-1)));other=np.delete(ev,neutral)
        row.update(rotation_order=order,full_period_candidate=float(order*fit.x[3]),relative_multipliers=[[float(a.real),float(a.imag)] for a in ev],nontrivial_moduli=np.abs(other).tolist(),trace_error=float(abs(np.log(abs(np.linalg.det(cache['M'])))/fit.x[3]+3*gamma)))
        sol=solve_ivp(field,[0,fit.x[3]],fit.x[:3],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.02,t_eval=np.linspace(0,fit.x[3],2000,endpoint=False));orbits[str(index)]=sol.y.T
    rows.append(row);(out/'report.json').write_text(json.dumps({'rows':rows,'seconds':time.time()-start_time,'limits':'Numerical relative shooting with rounded group matrices. Full-period closure and independent precision checks remain separate. No interval periodic-orbit or chaos proof.'},indent=2));np.savez_compressed(out/'orbits.npz',**orbits)
    print(json.dumps({k:v for k,v in row.items() if k!='guess'}),flush=True)
