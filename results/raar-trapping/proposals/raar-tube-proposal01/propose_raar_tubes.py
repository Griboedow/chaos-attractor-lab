"""Numerical local-frame tube diagnostics; no invariant-region certificate."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.spatial import cKDTree
from search_raar_circles import step
from trace_raar_cycle_seeds import value

@njit(cache=True)
def frames_and_radii(z,p,tangent,axis_fraction):
    n=len(z);frames=np.empty((n,3,3));h=np.empty(n);curvature=np.empty(n);normal_rate=np.empty(n);transport=np.empty(n);weight=np.empty(n);Q=np.eye(3);s=.001;w=1.
    for k in range(n):
        frames[k]=Q;transport[k]=s;weight[k]=w;y,J,axis,status=step(z[k],p);assert status==0;hk=min(tangent,axis_fraction*axis);h[k]=hk;V=J@Q;q1=V[:,0]/np.linalg.norm(V[:,0]);v2=V[:,1]-np.dot(V[:,1],q1)*q1;q2=v2/np.linalg.norm(v2);q3=np.cross(q1,q2);Qnext=np.column_stack((q1,q2,q3));rate=abs(np.dot(q3,V[:,2]));normal_rate[k]=rate;eta=0.
        for a in [-1.,0.,1.]:
            for b in [-1.,0.,1.]:
                delta=hk*(a*Q[:,0]+b*Q[:,1]);v,flag,_=value(z[k]+delta,p);assert flag==0;eta=max(eta,abs(np.dot(q3,v-y-J@delta)))
        curvature[k]=eta;s=rate*s+1.2*eta;w=rate*w+1;Q=Qnext
    return frames,h,curvature,normal_rate,transport,weight

@njit(cache=True,parallel=True)
def images(z,frames,h,normal,ids,grid,p):
    count=len(ids)*len(grid);points=np.empty((count,3));mapped=np.empty_like(points);flags=np.zeros(count,np.int64);axes=np.empty(count)
    for i in prange(len(ids)):
        index=ids[i];r=np.array([h[index],h[index],normal[index]])
        for j in range(len(grid)):
            x=z[index]+frames[index]@(r*grid[j]);y,flag,axis=value(x,p);k=i*len(grid)+j;points[k]=x;mapped[k]=y;flags[k]=flag;axes[k]=axis
    return points,mapped,flags,axes

@njit(cache=True,parallel=True)
def scores(mapped,neighbors,z,frames,h,normal):
    score=np.full(len(mapped),np.inf);choice=np.full(len(mapped),-1,np.int64)
    for i in prange(len(mapped)):
        for j in neighbors[i]:
            q=frames[j].T@(mapped[i]-z[j]);v=max(abs(q[0])/h[j],abs(q[1])/h[j],abs(q[2])/normal[j])
            if v<score[i]:score[i]=v;choice[i]=j
    return score,choice

def main(args):
    repo=Path(__file__).resolve().parents[1];args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6);begin=time.time();source=repo/'results/raar-circles/render/case-1051-animation.npz';d=np.load(source);z=d['states'];p=d['params'];frames,h,curvature,rate,transport,weight=frames_and_radii(z,p,args.tangent,args.axis_fraction)
    np.savez_compressed(args.out/'frames.npz',frames=frames,h=h,curvature=curvature,normal_rate=rate,transport=transport,weight=weight);valid=np.arange(10000,len(z));tree=cKDTree(z[valid]);ids=np.linspace(20000,len(z)-2,args.samples,dtype=int);grid=np.array(list(itertools.product([-1.,0.,1.],repeat=3)));records=[]
    print(json.dumps(dict(stage='frames',seconds=time.time()-begin,normal_rate_quantiles=np.quantile(rate[10000:],[0,.5,.99,1]).tolist(),transport_quantiles=np.quantile(transport[10000:],[0,.5,.99,1]).tolist())),flush=True)
    for epsilon in args.epsilon:
        normal=transport+epsilon*weight;points,mapped,flags,axes=images(z,frames,h,normal,ids,grid,p);_,nearest=tree.query(mapped,k=args.neighbors,workers=4);neighbors=valid[nearest];score,choice=scores(mapped,neighbors,z,frames,h,normal);key=format(epsilon,'.0e');np.savez_compressed(args.out/f'epsilon-{key}.npz',source_indices=ids,grid=grid,points=points,images=mapped,flags=flags,axes=axes,scores=score,targets=choice)
        row=dict(epsilon=epsilon,samples=args.samples,points=len(points),undefined=int(np.count_nonzero(flags)),inside=int(np.sum(score<1)),inside_fraction=float(np.mean(score<1)),all_vertices_inside_boxes=int(np.sum(np.max(score.reshape(len(ids),-1),axis=1)<1)),score_quantiles=np.quantile(score,[0,.5,.9,.99,1]).tolist(),minimum_sampled_axis=float(axes.min()),normal_radius_quantiles=np.quantile(normal[10000:],[0,.5,.99,1]).tolist(),seconds=time.time()-begin);records.append(row);print(json.dumps(row),flush=True)
    result=dict(tangent=args.tangent,axis_fraction=args.axis_fraction,neighbors=args.neighbors,rows=records,reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Finite local-frame and vertex diagnostics only. Tangent frames, transported radii, nonlinear remainder estimates and membership searches are numerical. No enclosure of whole boxes or proof of trapping/attraction.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--tangent',type=float,default=.01);p.add_argument('--axis-fraction',type=float,default=.08);p.add_argument('--samples',type=int,default=5000);p.add_argument('--neighbors',type=int,default=32);p.add_argument('--epsilon',type=float,nargs='+',default=[1e-6,1e-5,1e-4]);main(p.parse_args())
