"""Independent adaptive full-spectrum check of a saved curved-chain pilot."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from follower_link import solve3
from search_curved_follower import field


@njit(cache=True)
def field_jac(x,p):
    length=np.array([1.,p[5],p[6]])
    tail=np.array([1.+p[3]+p[4],p[3]+p[4],p[4]])
    rest=np.array([0.,p[7],p[7]+p[8]])
    elastic=np.array([[1.+p[1],-p[1],0.],[-p[1],p[1]+p[2],-p[2]],[0.,-p[2],p[2]]])
    m=np.empty((3,3));s=np.empty((3,3))
    for i in range(3):
        for j in range(3):
            base=tail[max(i,j)]*length[i]*length[j]
            m[i,j]=base*np.cos(x[i]-x[j]);s[i,j]=base*np.sin(x[i]-x[j])
    f=solve3(m,p[0]*length*np.sin(x-x[2])-elastic@(x-rest))
    jac=np.empty((3,3))
    for j in range(3):
        g=np.empty(3)
        for i in range(3):
            g[i]=p[0]*length[i]*np.cos(x[i]-x[2])*((i==j)-(2==j))-elastic[i,j]
            for a in range(3):
                g[i]+=s[i,a]*((i==j)-(a==j))*f[a]
        jac[:,j]=solve3(m,g)
    return f,jac


@njit(cache=True)
def augmented(z,p):
    f,j=field_jac(z[:3],p)
    ans=np.empty(13);ans[:3]=f;ans[3:12]=(j@z[3:12].reshape(3,3)).ravel();ans[12]=np.trace(j)
    return ans


def write(path,record):
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf8')


def crossings(states,p):
    # Strict crossings of nonadjacent zero-thickness segments; no near-contact test.
    direction=np.stack((np.cos(states),np.sin(states)),axis=2)*np.array([1.,p[5],p[6]])[None,:,None]
    positions=np.concatenate((np.zeros((len(states),1,2)),np.cumsum(direction,axis=1)),axis=1)
    def cross(a,b):return a[:,0]*b[:,1]-a[:,1]*b[:,0]
    a,b,c,d=positions[:,0],positions[:,1],positions[:,2],positions[:,3]
    intersect=(cross(b-a,c-a)*cross(b-a,d-a)<0)&(cross(d-c,a-c)*cross(d-c,b-c)<0)
    return dict(sampled_strict_intersections=int(intersect.sum()),samples=len(states),
                sampled_intersection_fraction=float(intersect.mean()),
                warning='Sampled zero-thickness segments only; absence does not certify positive clearance.')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);p=source['parameters'][args.id];x=source['states'][args.id,:3].copy()
    initial=x.copy();f,j=field_jac(x,p)
    reference=np.column_stack([field(x+1e-25j*np.eye(3)[i],p).imag/1e-25 for i in range(3)])
    error=float(np.max(abs(j-reference)));assert error<1e-9
    delta2=p[7]-p[0]*p[5]*np.sin(p[8])/p[1]
    eq1=p[0]*(-np.sin(delta2+p[8])-p[5]*np.sin(p[8]))
    equilibrium=np.array([eq1,eq1+delta2,eq1+delta2+p[8]])
    eq_residual=float(np.max(abs(field(equilibrium,p))));assert eq_residual<1e-9
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        id=args.id,parameters=p.tolist(),initial=initial.tolist(),warm=args.warm,duration=args.duration,
        method='DOP853 full 3D variational QR; independent general-matrix and complex-step controls',
        rtol=1e-10,atol=1e-12,max_step=args.step,qr=.1,sampling=.01,algebra_error=error,
        equilibrium=equilibrium.tolist(),equilibrium_residual=eq_residual,
        selection='High mean pilot rates examined even when the four-window retention test failed.'))
    for name in ['check_curved_follower.py','search_curved_follower.py','follower_link.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    basis=np.eye(3);logs=np.zeros(3);trace=0.;elapsed=0.;times=[];states=[];windows=[]
    previous=np.zeros(3);prevtrace=0.;section_times=[];section_states=[];began=time.perf_counter()
    def section(t,z):return z[0]-equilibrium[0]
    section.direction=1
    for iteration in range(int(round((args.warm+args.duration)/.1))):
        end=(iteration+1)*.1
        grid=np.linspace(elapsed,end,11)
        sol=solve_ivp(lambda t,z:augmented(z,p),(elapsed,end),np.r_[x,basis.ravel(),0.],
                      method='DOP853',rtol=1e-10,atol=1e-12,max_step=args.step,t_eval=grid,events=section)
        if not sol.success:raise RuntimeError(sol.message)
        x=sol.y[:3,-1];basis,r=np.linalg.qr(sol.y[3:12,-1].reshape(3,3))
        if elapsed>=args.warm-1e-9:
            logs+=np.log(abs(np.diag(r)));trace+=sol.y[-1,-1]
            times.extend(grid[1:]);states.extend(sol.y[:3,1:].T)
            section_times.extend(sol.t_events[0]);section_states.extend(np.asarray(sol.y_events[0]).reshape(-1,13)[:,:3])
        elapsed=end
        if elapsed>args.warm and abs((elapsed-args.warm)%(args.duration/4))<1e-8:
            window=dict(end=elapsed,spectrum=((logs-previous)/(args.duration/4)).tolist(),
                        mean_divergence=(trace-prevtrace)/(args.duration/4))
            windows.append(window);previous=logs.copy();prevtrace=trace
            print(args.id,window,flush=True)
    states=np.array(states)
    np.savez_compressed(out/'trace.npz',time=times,states=states,parameters=p,final_state=x,final_basis=basis,
                        log_sums=logs,integrated_trace=trace,section_times=section_times,section_states=np.array(section_states).reshape(-1,3))
    result=dict(spectrum=(logs/args.duration).tolist(),windows=windows,trace_residual_per_tau=float(abs(logs.sum()-trace)/args.duration),
                min_state=states.min(axis=0).tolist(),max_state=states.max(axis=0).tolist(),sections=len(section_times),
                crossing_check=crossings(states,p),wall_seconds=time.perf_counter()-began,
                status='Finite-time numerical evidence only; no historical or topological novelty claim.')
    write(out/'analysis.json',result)
    assert result['trace_residual_per_tau']<1e-6


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--id',type=int,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--warm',type=float,default=100.);p.add_argument('--duration',type=float,default=2000.)
    p.add_argument('--step',type=float,default=.01);main(p.parse_args())
