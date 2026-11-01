"""Numerical critical surfaces and postcritical images for the binary map.

This is a finite mesh and forward-orbit diagnostic, not a surface enclosure,
attractor boundary theorem, invariant-membership proof or novelty certificate.
"""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
from scipy.cluster.vq import kmeans2, vq
from search_product_alphabet import step


@njit(cache=True)
def critical_points(beta,tau,levels,lower,upper,size):
    points=np.empty((size*size*3*16,3));tags=np.empty(len(points),np.int64);count=0
    for axis in range(3):
        others=np.array([j for j in range(3) if j!=axis]);c=(levels[axis,0]+levels[axis,1])/2;a=(levels[axis,1]-levels[axis,0])/2
        for ia in range(size):
            for ib in range(size):
                x=np.zeros(3);x[others[0]]=lower[others[0]]+(upper[others[0]]-lower[others[0]])*ia/(size-1)
                x[others[1]]=lower[others[1]]+(upper[others[1]]-lower[others[1]])*ib/(size-1)
                previous=c-10*tau/(2*a);x[axis]=previous
                _,J,_,_,_,flag=step(x,beta,tau,levels,2);old=np.linalg.det(J)
                for scan in range(1,65):
                    current=c+(-10+20*scan/64)*tau/(2*a);x[axis]=current
                    _,J,_,_,_,flag=step(x,beta,tau,levels,2);new=np.linalg.det(J)
                    if old*new<0:
                        left,right=previous,current;fl=old
                        for _ in range(50):
                            middle=(left+right)/2;x[axis]=middle
                            _,J,_,_,_,flag=step(x,beta,tau,levels,2);fm=np.linalg.det(J)
                            if fl*fm<=0:right=middle
                            else:left=middle;fl=fm
                        x[axis]=(left+right)/2;points[count]=x;tags[count]=axis;count+=1
                    previous=current;old=new
    return points[:count],tags[:count]


@njit(cache=True)
def trace(points,beta,tau,levels,checkpoints):
    snapshots=np.empty((len(checkpoints),len(points),3));minimum=np.full(len(points),1e100)
    errors=np.empty(len(points));radii=np.empty(len(points))
    for k in range(len(points)):
        x=points[k].copy();idx=0
        _,J,r,_,_,flag=step(x,beta,tau,levels,2);errors[k]=abs(np.linalg.det(J));radii[k]=r
        for n in range(checkpoints[-1]+1):
            if n==checkpoints[idx]:
                snapshots[idx,k]=x;idx+=1
                if idx==len(checkpoints):break
            x,_,r,_,_,flag=step(x,beta,tau,levels,2);assert not flag
            minimum[k]=min(minimum[k],r)
    return snapshots,minimum,errors,radii


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in [Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:states=f['states'];beta,tau=f['params'];levels=f['levels'][:,:2]
    with np.load(args.cycle) as f:orbit=f['orbit']
    lower,upper=states.min(0)-args.pad,states.max(0)+args.pad
    points,tags=critical_points(beta,tau,levels,lower,upper,args.size)
    checkpoints=np.array([0,1,2,4,8,16,32,64,256,1024])
    snapshots,minima,errors,radii=trace(points,beta,tau,levels,checkpoints)
    tree=cKDTree(states);distances=np.array([tree.query(s,workers=1)[0] for s in snapshots])
    near=distances[0]<=args.near
    centroids,labels=kmeans2(states[::20],orbit,iter=30,minit='matrix')
    labels=vq(states,centroids)[0];transitions=np.bincount(4*labels[:-1]+labels[1:],minlength=16).reshape(4,4)
    with np.load(args.diagnostics) as f:diag=f['values']
    r=diag[:,0];dmax=diag[:,3:6].max(1)
    first=beta*np.maximum(1,abs(1-1/r))
    second=np.maximum(abs(1-2*beta),abs(1-2*beta+2*beta/r))
    norm_bound=first+second*dmax;contract=norm_bound<1
    np.savez_compressed(args.out/'critical-images.npz',points=points,tags=tags,checkpoints=checkpoints,snapshots=snapshots,minimum_radii=minima,source_residuals=errors,source_radii=radii,reference_distances=distances,params=[beta,tau],levels=levels)
    np.savez_compressed(args.out/'cluster-labels.npz',labels=labels,centroids=centroids,transition_counts=transitions)
    rows=[]
    for k,n in enumerate(checkpoints):
        rows.append(dict(iteration=int(n),distance_quantiles=np.quantile(distances[k],[0,.1,.5,.9,1]).tolist(),
                         near_source_distance_quantiles=np.quantile(distances[k,near],[0,.1,.5,.9,1]).tolist(),
                         coordinate_min=snapshots[k].min(0).tolist(),coordinate_max=snapshots[k].max(0).tolist()))
    report=dict(source_sha256=sha(args.source),diagnostics_sha256=sha(args.diagnostics),cycle_sha256=sha(args.cycle),
                grid_size=args.size,pad=args.pad,lower=lower.tolist(),upper=upper.tolist(),root_scan_h=[-10,10],root_scan_intervals=64,bisection_steps=50,
                critical_points=len(points),by_axis=np.bincount(tags,minlength=3).tolist(),maximum_determinant_residual=float(errors.max()),
                source_radius_minimum=float(radii.min()),minimum_followed_radius=float(minima.min()),
                near_threshold=args.near,near_sources=int(near.sum()),images=rows,
                numerical_clusters=dict(counts=np.bincount(labels,minlength=4).tolist(),centroids=centroids.tolist(),transitions=transitions.tolist(),
                                         next_phase_fraction=float(np.mean(labels[1:]==((labels[:-1]+1)%4)))),
                contraction_bound=dict(formula='beta*max(1,abs(1-1/r))+max(abs(1-2*beta),abs(1-2*beta+2*beta/r))*max_j(d_j)',
                                       strictly_below_one=int(contract.sum()),total=len(contract),
                                       leading_log_sum_inside=float(diag[contract,2].sum()),leading_log_sum_outside=float(diag[~contract,2].sum())),
                seconds=time.time()-begin,scope=__doc__)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report[k] for k in ['critical_points','maximum_determinant_residual','near_sources','minimum_followed_radius','numerical_clusters','contraction_bound','seconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','cycle','diagnostics','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--size',type=int,default=65);p.add_argument('--pad',type=float,default=.04);p.add_argument('--near',type=float,default=.005)
    main(p.parse_args())
