"""Independent adaptive NumPy implementation for the leading loop candidate."""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp

root=Path(__file__).parent;d=np.load(root/'icosa-v1/input.npz')
out=root/'icosa-dop853';out.mkdir(exist_ok=False)
ident=183;c=d['centers'][ident];u=d['axes'][ident];sigma=d['sigma'][ident];gamma=d['gamma'][ident]
skew=np.zeros((len(c),3,3));skew[:,0,1]=-u[:,2];skew[:,0,2]=u[:,1];skew[:,1,0]=u[:,2];skew[:,1,2]=-u[:,0];skew[:,2,0]=-u[:,1];skew[:,2,1]=u[:,0]
def field(t,state):
    z=state[:3];r=z-c;cross=np.cross(u,r);g=np.exp(-.5*(r*r).sum(1)/sigma**2)
    f=-gamma*z+(g[:,None]*cross).sum(0)
    jac=-gamma*np.eye(3)+(g[:,None,None]*(skew-cross[:,:,None]*r[:,None,:]/sigma**2)).sum(0)
    return np.r_[f,(jac@state[3:].reshape(3,3)).ravel()]
z=d['starts'][ident].copy();q=np.eye(3);acc=np.zeros(3);windows=[];points=[];evals=0;begin=time.time()
for block in range(4000):
    sol=solve_ivp(field,[0,2],np.r_[z,q.ravel()],method='DOP853',rtol=1e-10,atol=1e-12,max_step=.25)
    assert sol.success
    z=sol.y[:3,-1];q,r=np.linalg.qr(sol.y[3:,-1].reshape(3,3));logs=np.log(np.abs(np.diag(r)));evals+=sol.nfev
    if block>=1000:
        acc+=logs;points.append(z)
        if (block-1000+1)%500==0:
            windows.append((acc/1000).tolist());acc[:]=0
            print(json.dumps({'block':block+1,'last_window':windows[-1],'seconds':time.time()-begin}),flush=True)
row={'id':ident,'warm_time':2000,'observed_time':6000,'qr_interval':2,'rtol':1e-10,'atol':1e-12,'max_step':.25,'windows':windows,'mean':np.mean(windows,axis=0).tolist(),'trace_error':float(np.max(np.abs(np.sum(windows,axis=1)+3*gamma))),'nfev':evals,'seconds':time.time()-begin}
np.savez_compressed(out/'183.npz',points=np.array(points))
(out/'summary.json').write_text(json.dumps(row,indent=2));print(json.dumps(row),flush=True)
