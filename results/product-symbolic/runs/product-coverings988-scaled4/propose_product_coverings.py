"""Sampled one-step branching-window proposals for the fourth iterate of the product sphere map."""
import argparse,itertools,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from numba import njit,prange,set_num_threads
from product_return_map import iterate


@njit(cache=True,parallel=True)
def maps(states,beta,tau,points):
    result=np.empty_like(states)
    for i in prange(len(states)):result[i]=iterate(states[i],beta,tau,points,1)[0]
    return result


def propose(data,attempt,max_iterations):
    beta,tau=data['params'];points=data['points'];c=data['center'];B=data['basis'];L=attempt['steps'];s=attempt['parameter'];C=np.linalg.inv(B)
    orbit=[c+s*B[:,0]] if 'orbit' not in attempt else np.array(attempt['orbit']);jac=[];directions=[B[:,0]];u=B[:,0]
    for k in range(L):
        x,J,r,flag=iterate(orbit[k],beta,tau,points,1);assert not flag
        if 'orbit' not in attempt:orbit.append(x)
        jac.append(J);u=J@u;directions.append(u/np.linalg.norm(u))
    normals=[None]*(L+1);normals[L]=C[0]/np.linalg.norm(C[0])
    for k in range(L-1,0,-1):v=jac[k].T@normals[k+1];normals[k]=v/np.linalg.norm(v)
    centers=np.array([c]+list(orbit[1:L]));bases=[B]
    for k in range(1,L):
        n=normals[k];slow=B[:,1]-np.dot(B[:,1],n)*n;slow/=np.linalg.norm(slow);fast=np.cross(n,slow);bases.append(np.column_stack([directions[k],slow,fast]))
    bases=np.array(bases);inverses=np.linalg.inv(bases);condition=float(max(np.linalg.cond(b) for b in bases))
    edges=np.array([(0,0),(0,1)]+[(i,(i+1)%L) for i in range(1,L)]);linear=[];offset=[]
    for i,j in edges:
        y,J,r,flag=iterate(centers[i],beta,tau,points,1);assert not flag;linear.append(inverses[j]@J@bases[i]);offset.append(inverses[j]@(y-centers[j]))
    linear=np.array(linear);offset=np.array(offset);gap=abs((C@(centers[1]-c))[0]);upper=np.full((L,3),.03);upper[0,0]=.7*gap;upper[1,0]=.2*gap
    grid=np.array(list(itertools.product(np.linspace(-1,1,9),repeat=3)));errors=np.full((len(edges),3),1e-10);history=[];success=False;radii=None
    for iteration in range(max_iterations):
        Arows=[];rhs=[]
        for e,((i,j),A,d) in enumerate(zip(edges,linear,offset)):
            error=1.5*errors[e]+1e-10
            row=np.zeros((L,3));row[i,0]-=abs(A[0,0]);row[i,1:]+=abs(A[0,1:]);row[j,0]+=1.12;Arows.append(row.ravel());rhs.append(-abs(d[0])-error[0])
            for coordinate in [1,2]:
                row=np.zeros((L,3));row[i]+=abs(A[coordinate]);row[j,coordinate]-=.9;Arows.append(row.ravel());rhs.append(-abs(d[coordinate])-error[coordinate])
        row=np.zeros((L,3));row[0,0]=1;row[1]=abs(C[0]@bases[1]);Arows.append(row.ravel());rhs.append(.97*gap)
        # HiGHS' default absolute feasibility tolerance exceeds several window
        # widths. Solve in microunits and independently inspect every constraint.
        scale=1e6
        fit=linprog(1/upper.ravel(),A_ub=Arows,b_ub=np.array(rhs)*scale,bounds=[(1e-9*scale,h*scale) for h in upper.ravel()],method='highs',options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
        if not fit.success:history.append(dict(iteration=iteration,linear_feasible=False,message=fit.message));break
        radii=fit.x.reshape(L,3)/scale;checks=[]
        violation=float(np.max(np.array(Arows)@radii.ravel()-rhs))
        if violation>1e-14:history.append(dict(iteration=iteration,linear_feasible=False,maximum_constraint_violation=violation));break
        for e,((i,j),A,d) in enumerate(zip(edges,linear,offset)):
            q=grid*radii[i];image=(maps(centers[i]+q@bases[i].T,beta,tau,points)-centers[j])@inverses[j].T
            errors[e]=np.max(abs(image-d-q@A.T),axis=0);normalized=image/radii[j];orientation=np.sign(A[0,0]);left=normalized[grid[:,0]==-1,0];right=normalized[grid[:,0]==1,0]
            face=float(min(np.min(-orientation*left),np.min(orientation*right))-1);stable=float(1-np.max(abs(normalized[:,1:])));checks.append(dict(source=int(i),target=int(j),face_margin=face,stable_margin=stable))
        margin=min(min(q['face_margin'],q['stable_margin']) for q in checks);history.append(dict(iteration=iteration,linear_feasible=True,maximum_constraint_violation=violation,minimum_sampled_margin=margin,radii=radii.tolist(),sampled_remainders=errors.tolist(),checks=checks))
        if margin>.025:success=True;break
    record=dict(steps=L,parameter=s,grid_index=attempt['grid_index'],end_distance=attempt['end_distance'],maximum_basis_condition=condition,history=history,sampled_passed=success)
    arrays=None if radii is None else dict(params=data['params'],points=points,centers=centers,bases=bases,radii=radii,edges=edges,orientation=np.sign(linear[:,0,0]),linear=linear,offset=offset,orbit=np.array(orbit))
    return record,arrays


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time();data=dict(np.load(args.source/'trace.npz'));report=json.loads((args.source/'report.json').read_text(encoding='utf-8'));selected=[]
    for name in [Path(__file__).name,'product_return_map.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for n in sorted(set(r['steps'] for r in report['attempts'])):
        group=[r for r in report['attempts'] if r['steps']==n and r.get('root_residual',1)<1e-8 and r.get('accepted_numerical',True)]
        selected+=sorted(group,key=lambda r:r['end_distance'])[:args.per_length]
    rows=[]
    for ident,attempt in enumerate(selected):
        result,arrays=propose(data,attempt,args.iterations);result['id']=ident;rows.append(result)
        if arrays is not None:np.savez_compressed(args.out/f'windows-{ident}.npz',**arrays)
        print(json.dumps({k:v for k,v in result.items() if k!='history'}),flush=True)
        (args.out/'report.json').write_text(json.dumps(dict(source=str(args.source),source_sha256=hashlib.sha256((args.source/'trace.npz').read_bytes()).hexdigest(),selection='Smallest terminal distance, up to the declared count at each return length; every selected outcome retained.',per_length=args.per_length,grid=9,proposals=rows,seconds=time.time()-begin,scope='Linear programs and sampled nonlinear remainders only. No interval covering relation, invariant manifold, entropy or historical priority claim.'),indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'product_return_map.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--per-length',type=int,default=3);p.add_argument('--iterations',type=int,default=20);main(p.parse_args())
