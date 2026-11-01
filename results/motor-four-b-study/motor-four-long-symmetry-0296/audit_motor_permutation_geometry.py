"""Separate equation permutation symmetry from the sampled orbit's symmetry."""
import argparse
import hashlib
from itertools import permutations
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from motor_rows import field_jac


def full(x,n):return np.concatenate([x[...,:2*n],x[...,2*n:],-x[...,2*n:].sum(axis=-1,keepdims=True)],axis=-1)


def act(x,p,n):
    y=full(x,n);indices=np.r_[p,np.array(p)+n,np.array(p)+2*n]
    y=y[...,indices];return y[...,:3*n-1]


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    data=np.load(args.source);x=data['states'];params=data['parameters'];n=(x.shape[1]+1)//3;d=3*n-1
    assert x.shape[1]==d
    allperms=list(permutations(range(n)))
    weights=np.r_[np.ones(2*n),np.full(n,params[0]/params[1])]/np.sqrt(1+params[2]**2)
    a=full(x[:len(x)//2:2],n)*weights;b=full(x[len(x)//2::2],n)*weights
    trees=[cKDTree(a),cKDTree(b)];queries=[a[::37],b[::37]];records=[];max_equivariance=0.
    for p in allperms:
        # Check field/Jacobian covariance independently of cloud distances.
        matrix=np.column_stack([act(np.eye(d)[i],p,n) for i in range(d)])
        error=0.
        for state in x[[0,len(x)//2,-1]]:
            f,j=field_jac(state,params,n);fp,jp=field_jac(act(state,p,n),params,n)
            error=max(error,float(np.max(abs(fp-matrix@f)/(1+abs(fp)))),float(np.max(abs(jp@matrix-matrix@j)/(1+abs(matrix@j)))))
        assert error<1e-11;max_equivariance=max(max_equivariance,error)
        indices=np.r_[p,np.array(p)+n,np.array(p)+2*n]
        forward=trees[1].query(queries[0][:,indices],workers=1)[0]
        reverse=trees[0].query(queries[1][:,indices],workers=1)[0]
        combined=np.r_[forward,reverse]
        records.append(dict(permutation=list(p),median=float(np.median(combined)),p95=float(np.quantile(combined,.95)),
            maximum=float(combined.max()),equivariance_relative_error=error))
    baseline=records[0]
    for r in records:
        r['median_ratio_to_identity']=r['median']/baseline['median']
        r['p95_ratio_to_identity']=r['p95']/baseline['p95']
        r['within_three_times_identity']=bool(r['median']<=3*baseline['median'] and r['p95']<=3*baseline['p95'])
    record=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),N=n,dimension=d,
        metric='Euclidean distance of full(c,s,X) after c,s divided by sqrt(1+b^2), X additionally multiplied by rho/r. No learned projection.',
        reference='Each half of the full observation, stride2; queries every37th reference sample; map first-half queries to second-half cloud and conversely.',
        database_counts=[len(a),len(b)],query_counts=[len(q) for q in queries],scale=weights.tolist(),
        equation_equivariance_relative_error=max_equivariance,permutations=records,
        empirical_rule='Both median and95th-percentile distances <=3times the identity cross-half baseline. Heuristic, not a certified symmetry decision.',
        limitation='Finite point clouds cannot prove an invariant set symmetry, attraction, topology, equivalence or historical novelty.')
    (out/'audit.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    for name in ['audit_motor_permutation_geometry.py','motor_rows.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('baseline',baseline,flush=True)
    print('empirically close permutations',[r['permutation'] for r in records if r['within_three_times_identity']],flush=True)
    print('ratios',[(r['permutation'],round(r['median_ratio_to_identity'],2),round(r['p95_ratio_to_identity'],2)) for r in records],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
