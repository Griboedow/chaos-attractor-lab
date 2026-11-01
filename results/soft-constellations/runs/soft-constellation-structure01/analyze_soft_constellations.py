"""Standalone mean-shift controls and a declared split of the Case 29 cloud."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit
from search_soft_constellations import posterior


@njit(cache=True)
def mean_only(x,points,tau):
    residual=np.inf
    for k in range(10000):
        y,D,peak,entropy=posterior(x,points,tau);residual=np.linalg.norm(y-x);x=y
        if residual<1e-12:return x,k+1,residual
    return x,10000,residual


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[]
    for name in [Path(__file__).name,'search_soft_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for ident in [29,156]:
        source=args.source/f'case-{ident}.npz';cohort_source=args.cohort/f'case-{ident}.npz'
        with np.load(source) as f:points=f['points'];beta,tau=f['params'];z=f['states']
        with np.load(cohort_source) as f:starts=f['starts']
        baseline=[mean_only(x.copy(),points,tau) for x in starts]
        row=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cohort_sha256=hashlib.sha256(cohort_source.read_bytes()).hexdigest(),standalone_mean_control=[dict(initial=x.tolist(),last=y.tolist(),iterations=int(n),residual=float(r)) for x,(y,n,r) in zip(starts,baseline)])
        if ident==29:
            center=z.mean(0);_,v=np.linalg.eigh(np.cov(z.T));direction=v[:,-1];t=(z-center)@direction;q=np.sort(t);k=np.argmax(np.diff(q));cut=float((q[k]+q[k+1])/2)
            dense_source=args.dense/f'case-{ident}-states.npz'
            with np.load(dense_source) as f:dense=f['states']
            t=(dense-center)@direction;labels=(t>cut).astype(np.uint8)
            transitions=np.bincount(2*labels[:-1]+labels[1:],minlength=4).reshape(2,2)
            row['split']=dict(method='Largest empty interval along the first PCA direction in the 100000-state terminal long-run sample; tested on every consecutive state in the separate two-million-state render run.',center=center.tolist(),direction=direction.tolist(),cut=cut,dense_source_sha256=hashlib.sha256(dense_source.read_bytes()).hexdigest(),count=len(dense),class_counts=np.bincount(labels,minlength=2).tolist(),transitions=transitions.tolist(),observed_gap=[float(t[labels==0].max()),float(t[labels==1].min())],two_step_residual_quantiles=np.quantile(np.linalg.norm(dense[2:]-dense[:-2],axis=1),[0,.5,1]).tolist(),scope='Sampled alternating lobes, not an interval-certified invariant decomposition or a period-two orbit.')
        rows.append(row)
    report=dict(ids=[29,156],scope='Standalone averaging is a known mean-shift map. Stopping residuals and alternating labels are numerical observations only.',cases=rows)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps([dict(id=r['id'],baseline_max_iterations=max(x['iterations'] for x in r['standalone_mean_control']),baseline_max_residual=max(x['residual'] for x in r['standalone_mean_control']),split=r.get('split')) for r in rows]))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--cohort',type=Path,required=True);p.add_argument('--dense',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
