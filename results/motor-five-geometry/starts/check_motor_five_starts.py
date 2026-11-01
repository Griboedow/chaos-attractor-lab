"""Fresh five-row histories and independent short source-equation controls."""
import argparse,hashlib,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from motor_rows import source_force_constraint_field


@njit(cache=True)
def rhs(t,x,p):
    n=5;rho,r,b=p;v=np.empty(5);answer=np.empty(14);mean=np.sum(x[5:10])/5
    for i in range(5):
        pos=x[10+i] if i<4 else -np.sum(x[10:])
        v[i]=-rho*pos+r*(x[5+i]-mean)
        answer[i]=1-x[i]-v[i]*x[5+i];answer[5+i]=b-x[5+i]+v[i]*x[i]
        if i<4:answer[10+i]=v[i]
    return answer


def run(job):
    p=np.array(job['parameters']);x=np.array(job['initial']);grid=np.linspace(500,2000,30001)
    sol=solve_ivp(rhs,(0,2000),x,args=(p,),method='DOP853',rtol=2e-11,atol=2e-13,max_step=.01,t_eval=grid)
    assert sol.success
    path=Path(job['output'])/f"case-{job['id']:02}.npz"
    np.savez_compressed(path,time=grid,states=sol.y.T,parameters=p,initial=x)
    result=dict(id=job['id'],file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        maximum_disk_radius_squared=float(np.max(sol.y[:5]**2+sol.y[5:10]**2)))
    assert result['maximum_disk_radius_squared']<=1+p[2]**2+1e-8
    return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=np.load(args.source);p=d['parameters'];rng=np.random.default_rng(2026110117);starts=[]
    # Eight displacement-only starts and eight independent full disk starts.
    for i in range(16):
        X=rng.normal(size=5);X-=X.mean();X*=.1 if i<8 else 1.
        if i<8:c=np.ones(5);s=np.full(5,p[2])
        else:
            angle=rng.uniform(0,2*np.pi,5);rad=np.sqrt(1+p[2]**2)*np.sqrt(rng.uniform(0,1,5))
            c=rad*np.cos(angle);s=rad*np.sin(angle)
        starts.append(np.r_[c,s,X[:4]])
    controls=d['states'][np.linspace(0,len(d['states'])-1,4).astype(int)]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),seed=2026110117,
        parameters=p.tolist(),initial=np.array(starts).tolist(),source_control_initial=controls.tolist(),
        controls='Four uniformly indexed reference phases, independent force-constraint Radau vs compiled DOP853, both2e-12/2e-14,maxstep.002,time0..2.',
        histories='16 fresh histories, DOP8532e-11/2e-13,maxstep.01; warm500 then1500tau,dt.05.',
        scope='These states test attraction numerically, not a basin measure or asymptotic chaotic classification. All runs and original row labels retained.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_motor_five_starts.py','motor_rows.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    errors=[];paired=[];grid=np.linspace(0,2,201)
    for initial in controls:
        a=solve_ivp(lambda t,x:source_force_constraint_field(x,p,5),(0,2),initial,method='Radau',rtol=2e-12,atol=2e-14,max_step=.002,t_eval=grid)
        b=solve_ivp(rhs,(0,2),initial,args=(p,),method='DOP853',rtol=2e-12,atol=2e-14,max_step=.002,t_eval=grid)
        assert a.success and b.success
        err=float(np.max(abs(a.y-b.y)/(1+abs(a.y))));assert err<1e-8
        errors.append(err);paired.append([a.y.T,b.y.T])
    np.savez_compressed(out/'independent-controls.npz',time=grid,states=paired,initial=controls,parameters=p)
    (out/'independent-controls.json').write_text(json.dumps(dict(scaled_errors=errors,maximum=max(errors)),indent=2)+'\n')
    records=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs=[dict(id=i,parameters=p.tolist(),initial=x.tolist(),output=str(out)) for i,x in enumerate(starts)]
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();records.append(r)
            (out/'analysis.json').write_text(json.dumps(dict(declared=16,completed=len(records),runs=sorted(records,key=lambda x:x['id'])),indent=2)+'\n')
            print('completed',r['id'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
