"""Backward beam search for singular preimages near the observed cloud/core.

Only proposals: Newton convergence does not prove exhaustiveness, enclosure,
membership in an attractor, or density of the singular preimage set.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.spatial import cKDTree
from scipy.stats import qmc
from search_product_alphabet import step


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    for name in [Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:reference=f['states'];beta,tau=f['params'];levels=f['levels']
    prior=json.loads(args.mechanism.read_text());s=np.array(prior['central_singularity_numerical'])
    cycle=np.load(args.cycle)['orbit'];tree=cKDTree(reference);nodes=[];attempts=[];history=[]
    def add(x,depth,parent,residual):
        ident=len(nodes);nodes.append(dict(id=ident,depth=depth,parent=parent,state=x.tolist(),residual=float(residual),cloud_distance=float(tree.query(x)[0]),core_distance=float(np.linalg.norm(cycle-x,axis=1).min())))
        return ident
    add(s,0,None,0.);frontier=[0]
    common=1.6*qmc.Sobol(3,scramble=True,seed=2026110202).random_base2(args.seed_power)-.8
    # Every accepted original proposal is included, not just the two certified ones.
    first=[np.array(row['state']) for row in prior['preimages']]
    save=lambda name,v:(args.out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    for depth in range(1,args.depth+1):
        accepted=[]
        for parent in frontier:
            target=np.array(nodes[parent]['state']);_,nearest=tree.query(target,k=args.neighbors)
            starts=np.vstack([reference[np.maximum(nearest-1,0)],common]+([np.array(first)] if depth==1 else []))
            local=[]
            for initial in starts:
                fit=root(lambda x:step(x,beta,tau,levels,2)[0]-target,initial,jac=lambda x:step(x,beta,tau,levels,2)[1],tol=1e-11)
                result=step(fit.x,beta,tau,levels,2);error=float(np.max(abs(result[0]-target)));ident=-1
                if np.isfinite(error) and error<1e-10 and np.max(abs(fit.x))<2.2 and not result[-1]:
                    for old in local:
                        if np.linalg.norm(fit.x-np.array(nodes[old]['state']))<1e-7:ident=old;break
                    if ident<0:ident=add(fit.x,depth,parent,error);local.append(ident);accepted.append(ident)
                attempts.append([depth,parent,*initial,float(fit.success),error if np.isfinite(error) else 1e300,ident])
        # Preserve two objectives: proximity to the cloud and to the core phase.
        a=sorted(accepted,key=lambda k:nodes[k]['cloud_distance']);b=sorted(accepted,key=lambda k:nodes[k]['core_distance'])
        frontier=list(dict.fromkeys(a[:args.width//2]+b[:args.width//2]))
        row=dict(depth=depth,accepted=len(accepted),expanded_next=frontier,
                 nearest_cloud=nodes[a[0]] if a else None,nearest_core=nodes[b[0]] if b else None,seconds=time.time()-begin)
        history.append(row);save('nodes.json',nodes);save('progress.json',history)
        print(json.dumps(dict(depth=depth,accepted=len(accepted),minimum_cloud=row['nearest_cloud']['cloud_distance'] if a else None,minimum_core=row['nearest_core']['core_distance'] if b else None,seconds=row['seconds'])),flush=True)
        if not frontier:break
    np.savez_compressed(args.out/'attempts.npz',attempts=np.array(attempts),common_starts=common,params=[beta,tau],levels=levels)
    save('report.json',dict(status='completed',sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.source,args.mechanism,args.cycle]},depth=args.depth,width=args.width,neighbors=args.neighbors,seed_power=args.seed_power,nodes=len(nodes),attempts=len(attempts),history=history,seconds=time.time()-begin,scope=__doc__))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','mechanism','cycle','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--depth',type=int,default=12);p.add_argument('--width',type=int,default=12)
    p.add_argument('--neighbors',type=int,default=24);p.add_argument('--seed-power',type=int,default=6)
    main(p.parse_args())
