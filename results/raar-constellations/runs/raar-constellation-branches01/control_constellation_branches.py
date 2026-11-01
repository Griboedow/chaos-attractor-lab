"""Branch reduction, fixed-point diagnostics and a predeclared beta sweep."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.spatial import cKDTree
from search_raar_constellations import screen
from render_raar_constellations import orbit


def fixed_points(beta,points):
    rows=[]
    for i,b in enumerate(points):
        rho=np.linalg.norm(b[:2])
        row=dict(branch=i,radial_point_norm=float(rho),radial_minus_beta=float(rho-beta))
        if rho>beta:
            u=np.r_[b[:2]/rho,0.];x=(beta*u+(1-2*beta)*b)/(1-beta)
            distances=np.sum((points-x)**2,axis=1);other=int(np.argmin(distances))
            row.update(candidate=x.tolist(),nearest=other,branch_distance_excess=float(distances[i]-distances[other]),admissible=other==i)
        else:row['admissible']=False
        rows.append(row)
    return rows


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    protocol=dict(ids=args.ids,selected_after='Dense portraits, five-million-step checks and 33-start cohorts.',sweep=dict(beta_grid=[.85,.9995,300],warm=5000,observed=15000,windows=10,initial=[.317,-.211,.419]),reduction='Remove sites unvisited in the dense post-warmup sample; fresh run from the original initial state. Numerical comparison only; no global inactive-site proof.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    for name in [Path(__file__).name,'search_raar_constellations.py','render_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[];begin=time.time()
    for ident in args.ids:
        path=args.source/f'case-{ident}-states.npz'
        with np.load(path) as f:points=f['points'];params=f['params'];reference=f['states'];active=np.unique(f['branches'])
        beta=float(params[0]);reduced=points[active]
        z,branches,minima=orbit(beta,reduced,100000,1000000)
        tree=cKDTree(reference);scale=np.linalg.norm(np.ptp(reference,axis=0))
        distances=np.quantile(tree.query(z[::50],workers=1)[0]/scale,[.5,.99,1.])
        sweep_params=np.tile(params,(300,1));sweep_params[:,0]=np.linspace(.85,.9995,300)
        repeated=np.repeat(points[None],300,axis=0);counts=np.full(300,len(points))
        rates,tails,ends,status,mins,ratios=screen(sweep_params,repeated,counts,5000,15000)
        np.savez_compressed(args.out/f'case-{ident}-sweep.npz',params=sweep_params,points=points,rates=rates,tails=tails,ends=ends,status=status,minima=mins,boundratios=ratios)
        np.savez_compressed(args.out/f'case-{ident}-reduced.npz',states=z,branches=branches,points=reduced,active=active,params=params,minima=minima)
        row=dict(id=ident,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),active_sites=active.tolist(),reduced_distance_quantiles=distances.tolist(),fixed_point_diagnostics=fixed_points(beta,points),sweep_complete=int((status[:,0]==0).sum()),sweep_positive_all_windows=int(((status[:,0]==0)&(rates.min(1)>.005)).sum()))
        rows.append(row);(args.out/f'case-{ident}.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(row,seconds=time.time()-begin)),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,seconds=time.time()-begin),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())
