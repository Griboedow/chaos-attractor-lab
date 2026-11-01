"""Numerical window proposals for a branching return graph; no interval proof."""
import argparse,itertools,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog,brentq
from numba import njit,prange,set_num_threads
from continue_raar_cycles import iterate
from search_raar_returns import evaluate
from trace_raar_cycle_seeds import value

@njit(cache=True,parallel=True)
def maps(points,p):
    result=np.empty_like(points)
    for i in prange(len(points)):
        x=points[i]
        for _ in range(7):x,flag,axis=value(x,p)
        result[i]=x
    return result

def proposal(data,pre,post,s):
    p=data['params'];c=data['center'];B=data['basis'];lam=data['multipliers'][0];coeff=data['coefficients'];C=np.linalg.inv(B);L=pre+post
    seed=evaluate(coeff,s/lam**pre);points=[seed];jac=[];u=np.array([sum(k*coeff[i,k]*(s/lam**pre)**(k-1) for k in range(1,coeff.shape[1])) for i in range(3)]);directions=[u/np.linalg.norm(u)]
    for _ in range(L):
        x,J,_,_,status=iterate(points[-1],p,7);assert status==0;points.append(x);jac.append(J);u=J@u;directions.append(u/np.linalg.norm(u))
    normals=[None]*(L+1);normals[L]=C[0]/np.linalg.norm(C[0])
    for i in range(L-1,0,-1):v=jac[i].T@normals[i+1];normals[i]=v/np.linalg.norm(v)
    centers=[c]+points[1:L];bases=[B]
    for i in range(1,L):
        normal=normals[i];slow=B[:,1]-np.dot(B[:,1],normal)*normal;slow/=np.linalg.norm(slow);fast=np.cross(normal,slow);bases.append(np.column_stack([directions[i],slow,fast]))
    bases=np.array(bases);centers=np.array(centers);inverses=np.linalg.inv(bases);edges=[(0,0),(0,1)]+[(i,(i+1)%L) for i in range(1,L)];linear=[];offset=[]
    for i,j in edges:
        y,J,_,_,_=iterate(centers[i],p,7);linear.append(inverses[j]@J@bases[i]);offset.append(inverses[j]@(y-centers[j]))
    first_gap=abs((C@(centers[1]-c))[0]);upper=np.full((L,3),.035);upper[0,0]=.70*first_gap;upper[1,0]=.20*first_gap
    # Linear sufficient inequalities with slack, later checked on actual grids.
    rows=[];rhs=[]
    for (i,j),A,d in zip(edges,linear,offset):
        row=np.zeros((L,3));row[i,0]-=abs(A[0,0]);row[i,1:]+=abs(A[0,1:]);row[j,0]+=1.08;rows.append(row.ravel());rhs.append(-abs(d[0])-1e-8)
        for k in [1,2]:
            row=np.zeros((L,3));row[i]+=abs(A[k]);row[j,k]-=.94;rows.append(row.ravel());rhs.append(-abs(d[k])-1e-8)
    # Separate the fixed-point window and first excursion window along its unstable covector.
    row=np.zeros((L,3));row[0,0]=1;row[1]=abs(C[0]@bases[1]);rows.append(row.ravel());rhs.append(.97*first_gap)
    result=linprog(1/upper.ravel(),A_ub=np.array(rows),b_ub=rhs,bounds=[(1e-9,h) for h in upper.ravel()],method='highs')
    info=dict(pre=pre,post=post,length=L,parameter=s,linear_feasible=bool(result.success),optimizer_message=result.message,max_basis_condition=float(max(np.linalg.cond(v) for v in bases)))
    if not result.success:return info,None
    radii=result.x.reshape(L,3);grid=np.array(list(itertools.product(np.linspace(-1,1,9),repeat=3)));checks=[]
    for (i,j),A in zip(edges,linear):
        states=centers[i]+(grid*radii[i])@bases[i].T;image=(maps(states,p)-centers[j])@inverses[j].T/radii[j]
        left=image[grid[:,0]==-1,0];right=image[grid[:,0]==1,0];orientation=np.sign(A[0,0]);face=min(float(np.min(-orientation*left)),float(np.min(orientation*right)));stable=float(np.max(abs(image[:,1:])))
        checks.append(dict(source=i,target=j,face_margin=face-1,stable_margin=1-stable,passed=bool(face>1 and stable<1)))
    info.update(checks=checks,all_sampled_passed=all(r['passed'] for r in checks),min_sampled_margin=min(min(r['face_margin'],r['stable_margin']) for r in checks),radii=radii.tolist(),terminal_coordinates=(C@(points[-1]-c)).tolist())
    arrays=dict(centers=centers,bases=bases,radii=radii,params=p,edges=np.array(edges),orientation=np.sign([A[0,0] for A in linear]),linear=np.array(linear),offset=np.array(offset),orbit=np.array(points),seed=seed)
    return info,arrays

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(8);begin=time.time();data=dict(np.load(args.source/'trace.npz'));report=json.loads((args.source/'report.json').read_text());initial=[r for r in report['attempts'] if r.get('grid_index')==867 and r.get('returns')<=11 and r.get('root_residual',1)<1e-8];rows=[]
    for pre in [3,4,5,6]:
        for row in initial:
            info,arrays=proposal(data,pre,row['returns'],row['parameter']);rows.append(info)
            if arrays is not None:np.savez_compressed(args.out/f"pre-{pre}-post-{row['returns']}.npz",**arrays)
            print(json.dumps({k:v for k,v in info.items() if k not in ['checks','radii']}),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(cases=rows,seconds=time.time()-begin,scope='Linear-program and finite-grid proposals only. No covering relation or homoclinic/chaos theorem has been certified.'),indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
