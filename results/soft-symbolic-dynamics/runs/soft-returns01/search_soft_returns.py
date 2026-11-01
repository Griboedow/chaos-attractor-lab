"""Finite tangent-segment return proposals near a numerical saddle.

Neither the initial segment nor terminal plane is a validated invariant
manifold. Only later complete window coverings can establish symbolic dynamics.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.optimize import brentq
from find_soft_cycles import iterate


@njit(cache=True,parallel=True)
def trace(center,direction,grid,beta,tau,points,count):
    result=np.empty((count+1,len(grid),3));status=np.zeros(len(grid),np.int64)
    for i in prange(len(grid)):
        x=center+grid[i]*direction;result[0,i]=x
        for k in range(count):
            x,J,r,flag=iterate(x,beta,tau,points,1);result[k+1,i]=x
            if flag:status[i]=flag;break
    return result,status


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    with np.load(args.source) as f:beta,tau=f['params'];points=f['points'];center=f['orbit'][0];values=f['multipliers'];basis=f['eigenvectors'].real
    assert np.max(abs(values.imag))<1e-12 and abs(values[0])>1 and np.max(abs(values[1:]))<1
    for j in range(3):
        basis[:,j]/=np.linalg.norm(basis[:,j])
        if basis[np.argmax(abs(basis[:,j])),j]<0:basis[:,j]*=-1
    inverse=np.linalg.inv(basis);grid=np.linspace(-args.radius,args.radius,args.samples)
    states,status=trace(center,basis[:,0],grid,beta,tau,points,args.steps);assert not status.any();coordinates=(states-center)@inverse.T
    np.savez_compressed(args.out/'trace.npz',params=[beta,tau],points=points,center=center,basis=basis,multipliers=values.real,grid=grid,states=states,coordinates=coordinates,status=status)
    attempts=[]
    for k in range(3,args.steps+1):
        u=coordinates[k,:,0]
        for i in np.where(u[:-1]*u[1:]<0)[0]:
            if max(abs(grid[i]),abs(grid[i+1]))<1e-9:continue
            if min(np.linalg.norm(coordinates[k,i,1:]),np.linalg.norm(coordinates[k,i+1,1:]))>args.max_distance:continue
            if np.linalg.norm(states[k,i+1]-states[k,i])>.1:continue
            try:
                fun=lambda s:float((inverse@(iterate(center+s*basis[:,0],beta,tau,points,k)[0]-center))[0])
                s=brentq(fun,grid[i],grid[i+1],xtol=1e-18,rtol=1e-14)
                x=center+s*basis[:,0];orbit=[x]
                for j in range(k):x=iterate(x,beta,tau,points,1)[0];orbit.append(x)
                orbit=np.array(orbit);end=inverse@(x-center);excursion=float(np.max(np.linalg.norm(orbit-center,axis=1)))
                if excursion<.04:continue
                attempts.append(dict(steps=k,grid_index=int(i),parameter=float(s),end_coordinates=end.tolist(),end_distance=float(np.linalg.norm(x-center)),max_excursion=excursion,root_residual=float(abs(end[0]))))
            except (ValueError,RuntimeError) as exc:attempts.append(dict(steps=k,grid_index=int(i),error=str(exc)))
        group=[r for r in attempts if r['steps']==k and 'end_distance' in r]
        if group:print(json.dumps(dict(steps=k,roots=len(group),best=min(group,key=lambda r:r['end_distance']))),flush=True)
    ordered=sorted([r for r in attempts if 'end_distance' in r],key=lambda r:r['end_distance'])
    result=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),radius=args.radius,samples=args.samples,steps=args.steps,attempts=attempts,closest=ordered[:30],seconds=time.time()-begin,scope='Numerical tangent-segment returns to a linear plane only, not a homoclinic orbit or entropy certificate.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'find_soft_cycles.py','search_soft_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--radius',type=float,default=.002);p.add_argument('--samples',type=int,default=16385);p.add_argument('--steps',type=int,default=18);p.add_argument('--max-distance',type=float,default=.02);main(p.parse_args())
