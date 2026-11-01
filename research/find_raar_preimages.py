"""Numerical multistart search for colliding preimages near the observed R1051 cloud."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.spatial import cKDTree
from scipy.stats import qmc
from search_raar_circles import step


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    with np.load(args.source) as d:z=d['states'];p=d['params']
    tree=cKDTree(z[1:]);reference=cKDTree(z);target_indices=np.linspace(10000,len(z)-10001,25,dtype=int);records=[];certificates=[]
    protocol=dict(source=str(args.source),source_sha256=sha(args.source),target_indices=target_indices.tolist(),nearest_output_candidates=256,global_seeds=args.global_seeds,global_seed=2026110196,global_box='target/beta plus cube of halfwidth 1+abs((1-2*beta)/beta)*(norm(center)+radius); includes the elementary preimage norm bound.',seed_separation=.005,root_tolerance=1e-11,accepted_residual=1e-9,distinct_preimage_separation=1e-6,scope='Numerical proposal search only. Targets and accepted preimages require independent interval validation; finite proximity to a reference cloud is not attractor membership. Global seeds, when requested, are an adaptive follow-up to the near-cloud-only search.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    for name in ['find_raar_preimages.py','search_raar_circles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for index in target_indices:
        target=z[index];neighbors=tree.query(target,k=256)[1];starts=[]
        for n in neighbors:
            point=z[n]
            if not starts or min(np.linalg.norm(point-s) for s in starts)>.005:starts.append(point)
        if args.global_seeds:
            assert args.global_seeds&(args.global_seeds-1)==0
            u=qmc.Sobol(3,scramble=True,seed=2026110196).random_base2(int(np.log2(args.global_seeds)));half=1+abs((1-2*p[0])/p[0])*(np.linalg.norm(p[3:])+p[1]);starts.extend(target/p[0]+half*(2*u-1))
        solutions=[];attempts=[]
        def fun(x):
            y,J,_,_=step(x,p);return y-target,J
        for start in starts:
            solution=root(fun,start,jac=True,tol=1e-11);y,J,axis,status=step(solution.x,p);res=float(np.max(np.abs(y-target)));good=not status and res<1e-9 and np.isfinite(solution.x).all()
            attempts.append(dict(start=start.tolist(),solution=solution.x.tolist(),residual=res,good=bool(good),solver_success=bool(solution.success)))
            if good and (not solutions or min(np.linalg.norm(solution.x-s['state']) for s in solutions)>1e-6):solutions.append(dict(state=solution.x.tolist(),residual=res,determinant=float(np.linalg.det(J)),reference_distance=float(reference.query(solution.x)[0]),axis_distance=float(axis)))
        record=dict(index=int(index),target=target.tolist(),params=p.tolist(),solutions=solutions,attempts=attempts);records.append(record)
        if len(solutions)>=2:certificates.append(dict(label=f'target-{index}',params=p.tolist(),target=target.tolist(),preimages=[s['state'] for s in solutions]))
        print(json.dumps(dict(index=int(index),starts=len(starts),solutions=len(solutions),reference_distances=[s['reference_distance'] for s in solutions])),flush=True)
    report=dict(protocol,targets=records,seconds=time.time()-begin);(args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(args.out/'certificate-input.json').write_text(json.dumps(dict(cases=certificates),indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=Path('results/raar-circles/render/case-1051-animation.npz'));p.add_argument('--global-seeds',type=int,default=0);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
