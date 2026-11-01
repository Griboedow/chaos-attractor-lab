"""Compare every fresh five-row history under all 120 fixed row relabelings."""
import argparse,hashlib,json
from itertools import permutations
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from audit_motor_permutation_geometry import full


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    ref=np.load(args.reference);p=ref['parameters'];w=np.r_[np.ones(10),np.full(5,p[0]/p[1])]/np.sqrt(1+p[2]**2)
    tree=cKDTree(full(ref['states'][::4],5)*w);perms=list(permutations(range(5)))
    files=sorted(args.source.glob('case-*.npz'));assert len(files)==16
    protocol=dict(reference=str(args.reference),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        source=str(args.source),histories=16,reference_stride=4,query_stride=301,permutations=120,weights=w.tolist(),
        metric='Full15D physical(c,s,X) Euclidean point-cloud distance, fixed row permutation for each comparison.',
        scope='One-sided sampled support comparison, not attractor identity, basin probability or proof of asymptotic chaos.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['compare_motor_five_starts.py','audit_motor_permutation_geometry.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    for path in files:
        d=np.load(path);q=full(d['states'][::301],5)*w;rows=[]
        for perm in perms:
            ids=np.r_[perm,np.array(perm)+5,np.array(perm)+10];dist=tree.query(q[:,ids],workers=1)[0]
            rows.append(dict(permutation=list(perm),median=float(np.median(dist)),p95=float(np.quantile(dist,.95)),maximum=float(dist.max())))
        best=min(rows,key=lambda r:r['median'])
        records.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),queries=len(q),best=best,comparisons=rows))
        (out/'analysis.json').write_text(json.dumps(dict(completed=len(records),declared=16,runs=records),indent=2)+'\n')
        print(path.name,'best',best,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
