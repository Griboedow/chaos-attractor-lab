"""Independent fresh-start histories and full-state permutation comparisons.

No inference of common attractors from leading rates alone. Every declared
anchor start is integrated; all 24 row relabelings are compared with one saved
reference trajectory in physically scaled full state coordinates.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
from itertools import permutations
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from motor_rows import source_force_constraint_field
from audit_motor_permutation_geometry import full


@njit(cache=True)
def rhs(t, x, p):
    rho,r,b=p
    velocity=-rho*np.array([x[8],x[9],x[10],-sum(x[8:])])+r*(x[4:8]-sum(x[4:8])/4)
    return np.concatenate((1-x[:4]-velocity*x[4:8],b-x[4:8]+velocity*x[:4],velocity[:3]))


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')


def run(job):
    began=time.perf_counter();p=np.array(job['parameters']);initial=np.array(job['initial'])
    reference=source_force_constraint_field(initial,p,4)
    error=float(np.max(abs(rhs(0.,initial,p)-reference)/(1+abs(reference))))
    assert error<1e-11
    grid=np.linspace(500.,3500.,60001)
    sol=solve_ivp(rhs,(0.,3500.),initial,args=(p,),method='DOP853',rtol=2e-10,atol=2e-12,max_step=.03,t_eval=grid)
    assert sol.success
    x=sol.y.T;positions=full(x,4)[:,8:];velocity=-p[0]*positions+p[1]*(x[:,4:8]-x[:,4:8].mean(axis=1)[:,None])
    means=positions.mean(axis=0);distinguished=int(np.argmin(means));others=[i for i in range(4) if i!=distinguished]
    i,j,k=others
    circulation=float(np.mean((positions[:,i]-positions[:,k])*(velocity[:,j]-velocity[:,k])-(positions[:,j]-positions[:,k])*(velocity[:,i]-velocity[:,k])))
    crossings=[]
    for row in range(4):
        ids=np.flatnonzero((positions[:-1,row]<0)&(positions[1:,row]>=0))
        event=grid[ids]-.05*positions[ids,row]/(positions[ids+1,row]-positions[ids,row])
        crossings.append(dict(row=row,events=len(ids),mean_interval=float(np.diff(event).mean()) if len(event)>1 else None))
    path=Path(job['output'])/f"case-{job['id']:02}.npz"
    np.savez_compressed(path,time=grid,states=x,parameters=p,initial=initial)
    result=dict(id=job['id'],parameters=p.tolist(),initial=initial.tolist(),mean_position=means.tolist(),
        quarter_mean_positions=[part.mean(axis=0).tolist() for part in np.array_split(positions,4)],
        distinguished_by_smallest_mean_position=distinguished,circulation_of_sorted_other_rows=circulation,
        crossings=crossings,algebra_relative_error=error,maximum_disk_radius_squared=float(np.max(x[:,:4]**2+x[:,4:8]**2)),
        file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-began)
    assert result['maximum_disk_radius_squared']<1+p[2]**2+1e-8
    write(path.with_suffix('.json'),result)
    return result


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    inputs=np.load(args.starts);reference=np.load(args.reference)
    jobs=[dict(id=i,parameters=inputs['parameters'][i].tolist(),initial=inputs['initial'][i].tolist(),output=str(out)) for i in range(32)]
    write(out/'protocol.json',dict(starts=str(args.starts),starts_sha256=hashlib.sha256(args.starts.read_bytes()).hexdigest(),
        reference=str(args.reference),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),jobs=jobs,
        method='DOP853 state-only,rtol2e-10,atol2e-12,maxstep.03; warm500 then3000tau,uniform.05samples.',
        comparison='Reference full12ktau stride2; every61st fresh-start sample; all24 row permutations; same full-state scaling as permutation audit.',
        limitation='Finite numerical histories; no basin probability, attractor identity proof, or certified isotropy group.'))
    for name in ['check_motor_four_starts.py','motor_rows.py','audit_motor_permutation_geometry.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();records.append(r)
            write(out/'analysis.json',dict(declared=32,completed=len(records),runs=sorted(records,key=lambda x:x['id'])))
            print('completed',r['id'],'distinguished',r['distinguished_by_smallest_mean_position'],'circulation',round(r['circulation_of_sorted_other_rows'],4),flush=True)
    p=reference['parameters'];weights=np.r_[np.ones(8),np.full(4,p[0]/p[1])]/np.sqrt(1+p[2]**2)
    tree=cKDTree(full(reference['states'][::2],4)*weights)
    for r in sorted(records,key=lambda x:x['id']):
        x=np.load(out/r['file'])['states'][::61];query=full(x,4)*weights;distances=[]
        for perm in permutations(range(4)):
            ids=np.r_[perm,np.array(perm)+4,np.array(perm)+8]
            dd=tree.query(query[:,ids],workers=1)[0]
            distances.append(dict(permutation=list(perm),median=float(np.median(dd)),p95=float(np.quantile(dd,.95)),maximum=float(dd.max())))
        r['comparisons']=distances;r['best_by_median']=min(distances,key=lambda a:a['median'])
        write(out/'comparisons.json',dict(completed=sum('comparisons' in r for r in records),declared=32,runs=sorted(records,key=lambda a:a['id'])))
        print('compared',r['id'],r['best_by_median'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--starts',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
