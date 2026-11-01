"""Select the closest opposite-determinant sample pair; no invariant-set claim."""
import argparse,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.spatial import cKDTree
from scipy.optimize import brentq
from search_raar_circles import step
from verify_raar_folding import sha

@njit(cache=True,parallel=True)
def determinants(z,p):
    out=np.empty(len(z));flags=np.empty(len(z),np.int64)
    for i in prange(len(z)):
        _,j,_,flag=step(z[i],p);out[i]=np.linalg.det(j);flags[i]=flag
    return out,flags
def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(8)
    with np.load(args.source) as d:z=d['states'];p=d['params']
    det,flag=determinants(z,p);assert not flag.any();positive=np.flatnonzero(det>0);negative=np.flatnonzero(det<0);distance,nearest=cKDTree(z[positive]).query(z[negative],workers=8);which=np.argmin(distance);ia=int(negative[which]);ib=int(positive[nearest[which]]);a=z[ia];b=z[ib];t=brentq(lambda t:np.linalg.det(step(a+t*(b-a),p)[1]),0,1,xtol=1e-14)
    row=dict(source=str(args.source),source_sha256=sha(args.source),params=p.tolist(),indices=[ia,ib],segment_endpoints=[a.tolist(),b.tolist()],endpoint_determinants=[float(det[ia]),float(det[ib])],pair_distance=float(distance[which]),t_proposal=float(t),negative_count=len(negative),positive_count=len(positive),selection='Closest opposite determinant-sign pair among all one million fixed reference states, Euclidean metric. Adaptive structural investigation after inverse-image counts.',scope='Proposal only. Proximity to a finite reference sample is not attractor membership or a proof that the fold creates the observed geometry.')
    np.savez_compressed(args.out/'determinants.npz',determinants=det,flags=flag,nearest_opposite_distances=distance,negative_indices=negative,nearest_positive_indices=positive[nearest]);(args.out/'input.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
    for name in ['propose_raar_fold.py','search_raar_circles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(row),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=Path('results/raar-circles/render/case-1051-animation.npz'));p.add_argument('--out',type=Path,required=True);main(p.parse_args())
