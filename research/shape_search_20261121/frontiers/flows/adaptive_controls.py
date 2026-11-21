"""Independent fifth returns from a trigonometric polar RHS and two adaptive solvers."""
import json
import sys
import numpy as np
from scipy.integrate import solve_ivp
from native_section import ROOT,RAW,RESULT,P,packed
from return_tube import graph,CENTER,ZLO,ZHI,WIDTH,grid
sys.path.insert(0,str(ROOT/'research/shape_search_20261121/flows'))
from adaptive import polar_ring

def direct(x,z,method,rtol):
    def f(t,u):return polar_ring(u,P)
    def event(t,u):return u[1]
    event.direction=1
    sol=solve_ivp(f,(0,49),[x,0,z],method=method,events=event,rtol=rtol,atol=rtol*.01)
    ts=sol.t_events[0];es=sol.y_events[0]
    ts,es=ts[ts>1e-5],es[ts>1e-5]
    if not sol.success or len(ts)!=5:raise ValueError((sol.message,ts))
    normal=np.array([f(0,u)[1] for u in es])
    return np.r_[ts[-1],es[-1,[0,2]],normal.min(),np.diff(np.r_[0,ts]).min()]

def main():
    g=graph();rng=np.random.default_rng(91844)
    z=np.r_[np.linspace(ZLO,ZHI,30),rng.uniform(ZLO,ZHI,50)]
    d=np.r_[np.resize([-WIDTH,WIDTH,0.],30),rng.uniform(-WIDTH,WIDTH,50)]
    x=g(z-CENTER)+d;results={}
    for method,rtol in [('DOP853',2e-12),('RK45',3e-12)]:
        n=len(x) if method=='DOP853' else 18
        a=np.array([direct(x[i],z[i],method,rtol) for i in range(n)])
        np.savez_compressed(RAW/f'adaptive-{method}.npz',initial=np.c_[x[:n],z[:n]],result=a)
        rk=grid(x[:n],z[:n],.00625,packed())
        b=a[:,2];dd=a[:,1]-g(b-CENTER)
        results[method]=dict(rtol=rtol,n=n,step_reference=.00625,max_RK4_discrepancy=np.max(abs(rk[:,:3]-a[:,:3]),axis=0).tolist(),
            z_image=[b.min(),b.max()],delta_image=[dd.min(),dd.max()],min_positive_normal_speed=a[:,3].min(),
            minimum_first_and_adjacent_return_gap=a[:,4].min())
        if method=='RK45':
            other=np.load(RAW/'adaptive-DOP853.npz')['result'][:n,:3]
            results['adaptive_pair_max_discrepancy_t_x_z']=np.max(abs(other-a[:,:3]),axis=0).tolist()
        (RAW/'adaptive.json').write_text(json.dumps(results,indent=2)+'\n')
        print(method,results[method],flush=True)

if __name__=='__main__':main()
