"""Map identities and finite-orbit geometric diagnostics for RAAR circles."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from search_raar_circles import step,project


@njit(cache=True)
def inspect(params,start,warm,steps):
    beta,r,tilt=params[:3];c=params[3:];n=np.array([0.,-np.sin(tilt),np.cos(tilt)])
    x=start.copy();minaxis=np.inf;maxnorm=0.;negative=0;max_det_error=0.;max_eigen_error=0.;max_eigen_relative=0.
    stats=np.empty((steps//100,6));index=0
    for k in range(warm+steps):
        new,j,dist,status=step(x,params)
        if status:raise ValueError('Projection axis')
        if k>=warm:
            pb,_,lb,_=project(x,c,n,r);pa,_,la,_=project(2*pb-x,np.zeros(3),np.array([0.,0.,1.]),1.)
            u=np.array([-pa[1],pa[0],0.]);v=np.cross(n,(pb-c)/r);dot=np.dot(u,v)
            a=1/la;b=r/lb;det=beta*beta*((1-a)*(beta+(1-2*beta)*b)+a*b*dot*dot)
            actual=np.linalg.det(j);max_det_error=max(max_det_error,abs(det-actual)/max(1.,abs(actual)))
            # A left eigenvector: w^T J = beta*w^T. Right eigenvectors need not equal w.
            w=np.cross(u,v);error=np.max(np.abs(w@j-beta*w));max_eigen_error=max(max_eigen_error,error);max_eigen_relative=max(max_eigen_relative,error/max(1.,np.max(np.abs(j))))
            negative+=int(det<0);minaxis=min(minaxis,dist);maxnorm=max(maxnorm,np.linalg.norm(x))
            if (k-warm)%100==0:
                stats[index]=np.array([la,lb,det,np.linalg.norm(new-x),np.linalg.norm(pa-pb),np.linalg.norm(x)])
                index+=1
        x=new
    return stats,minaxis,maxnorm,negative,max_det_error,max_eigen_error,max_eigen_relative


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    with np.load(args.source) as f:params=f['params']
    rows=[]
    for ident in [75,820,1051,1311,1483,1679]:
        begin=time.time();p=params[ident];s=np.array([.317,-.211,.419]);d=inspect(p,s,100000,1000000)
        stats,axis,norm,negative,deterr,eigerr,eigrel=d;assert deterr<1e-8 and eigrel<1e-9
        beta,r,tilt,cx,cy,cz=p;C=beta+abs(1-2*beta)*(np.linalg.norm(p[3:])+r);bound=C/(1-beta)
        assert norm<=max(np.linalg.norm(s),bound)*(1+1e-12)
        np.savez_compressed(args.out/f'case-{ident}.npz',params=p,statistics=stats)
        rows.append(dict(id=ident,params=p.tolist(),warm=100000,observed=1000000,statistics_stride=100,
            statistic_columns=['axis_distance_A_at_RB(x)','axis_distance_B_at_x','determinant','step_norm','projection_gap','state_norm'],
            min_axis_distance_all_observed=axis,max_state_norm=norm,affine_bound_C=C,asymptotic_radius_bound=bound,
            negative_determinant_count=int(negative),negative_determinant_fraction=negative/1000000,
            max_relative_determinant_formula_residual=deterr,max_left_eigenvector_formula_residual=eigerr,max_relative_left_eigenvector_residual=eigrel,
            statistics_quantiles=np.quantile(stats,[0,.01,.5,.99,1],axis=0).tolist(),seconds=time.time()-begin))
        print(json.dumps(rows[-1]),flush=True)
    for n in ['audit_raar_geometry.py','search_raar_circles.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    (args.out/'report.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),runs=rows,scope='Analytic identities checked numerically. All orbit bounds and axis gaps are finite observations; the triangle-inequality absorbing-radius estimate is an independent elementary argument.'),indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
