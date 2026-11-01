"""Refine numerical window sizes using sampled nonlinear remainders."""
import argparse,itertools,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from numba import set_num_threads
from propose_raar_coverings import maps

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);data=dict(np.load(args.input));centers=data['centers'];bases=data['bases'];inverse=np.linalg.inv(bases);p=data['params'];edges=data['edges'];linear=data['linear'];offset=data['offset'];radii=data['radii'];L=len(centers);set_num_threads(8);begin=time.time();history=[]
    gap=abs((inverse[0]@(centers[1]-centers[0]))[0]);upper=np.full((L,3),.05);upper[0,0]=.70*gap;upper[1,0]=.20*gap
    grid=np.array(list(itertools.product(np.linspace(-1,1,17),repeat=3)))
    for iteration in range(30):
        remainders=[];checks=[]
        for (i,j),A,d in zip(edges,linear,offset):
            q=grid*radii[i];states=centers[i]+q@bases[i].T;image=(maps(states,p)-centers[j])@inverse[j].T
            remainders.append(np.max(abs(image-d-q@A.T),axis=0));normalized=image/radii[j];orientation=np.sign(A[0,0]);left=normalized[grid[:,0]==-1,0];right=normalized[grid[:,0]==1,0];face=min(np.min(-orientation*left),np.min(orientation*right))-1;stable=1-np.max(abs(normalized[:,1:]));checks.append(dict(source=int(i),target=int(j),face_margin=float(face),stable_margin=float(stable)))
        remainder=np.array(remainders);minimum=min(min(r['face_margin'],r['stable_margin']) for r in checks);history.append(dict(iteration=iteration,radii=radii.tolist(),sampled_remainders=remainder.tolist(),checks=checks,min_margin=minimum));print(json.dumps(dict(iteration=iteration,min_margin=minimum,max_radius=float(radii.max()))),flush=True)
        if minimum>.025:
            data['radii']=radii;np.savez_compressed(args.out/'windows.npz',**data);break
        rows=[];rhs=[]
        for k,((i,j),A,d) in enumerate(zip(edges,linear,offset)):
            error=1.5*remainder[k]+1e-9
            row=np.zeros((L,3));row[i,0]-=abs(A[0,0]);row[i,1:]+=abs(A[0,1:]);row[j,0]+=1.12;rows.append(row.ravel());rhs.append(-abs(d[0])-error[0])
            for coord in [1,2]:
                row=np.zeros((L,3));row[i]+=abs(A[coord]);row[j,coord]-=.9;rows.append(row.ravel());rhs.append(-abs(d[coord])-error[coord])
        row=np.zeros((L,3));row[0,0]=1;row[1]=abs(inverse[0,0]@bases[1]);rows.append(row.ravel());rhs.append(.97*gap)
        fit=linprog(1/upper.ravel(),A_ub=rows,b_ub=rhs,bounds=[(1e-8,h) for h in upper.ravel()],method='highs')
        if not fit.success:history[-1]['optimizer_failure']=fit.message;break
        radii=fit.x.reshape(L,3)
    (args.out/'report.json').write_text(json.dumps(dict(history=history,seconds=time.time()-begin,scope='Sampled remainders and grids propose windows. They are not interval bounds or covering certificates.'),indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
