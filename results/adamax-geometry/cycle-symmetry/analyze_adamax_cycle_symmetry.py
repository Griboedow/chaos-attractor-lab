"""Separate sign-related68-cycle destinations without counting copies as novelty."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from adamax_map import step_tangent


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=np.load(args.trajectories);closed=np.load(args.cycle);orbit=closed['orbit'];p=closed['parameters'];reflection=orbit*np.array([-1,-1,1])
    tails=d['tails'];unique,inverse=np.unique(tails.reshape(-1,3),axis=0,return_inverse=True)
    a=cKDTree(orbit).query(unique)[0][inverse].reshape(tails.shape[:2]).max(1)
    b=cKDTree(reflection).query(unique)[0][inverse].reshape(tails.shape[:2]).max(1)
    labels=(b<a).astype(np.int8);error=float(np.minimum(a,b).max());assert error<1e-10
    z=reflection[0].copy();product=np.eye(3)
    for _ in range(len(reflection)):z,J,_=step_tangent(z,np.eye(3),p);product=J@product
    R=np.diag([-1.,-1.,1.]);matrix_error=float(np.max(abs(product-R@closed['monodromy']@R)))
    closure=float(np.max(abs(z-reflection[0])));assert closure<1e-10 and matrix_error<1e-9
    np.savez_compressed(out/'symmetry.npz',original_cycle=orbit,reflected_cycle=reflection,labels=labels,
        maximum_tail_distance_to_original=a,maximum_tail_distance_to_reflected=b,reflected_monodromy=product)
    result=dict(period=len(orbit),records=len(labels),original_count=int((labels==0).sum()),reflected_count=int((labels==1).sum()),
        maximum_distance_to_assigned_cycle=error,reflected_closure=closure,monodromy_symmetry_error=matrix_error,
        distance_between_cycle_sets=float(cKDTree(orbit).query(reflection)[0].max()),
        trajectory_sha256=hashlib.sha256(args.trajectories.read_bytes()).hexdigest(),cycle_sha256=hashlib.sha256(args.cycle.read_bytes()).hexdigest(),
        scope='Two distinct sign-related68-cycles. Their symmetry equivalence is explicit; they are not counted as separate novel mechanisms. Counts concern continuation records, not independent initial conditions or basin probabilities.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in ['analyze_adamax_cycle_symmetry.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--trajectories',type=Path,required=True);p.add_argument('--cycle',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
