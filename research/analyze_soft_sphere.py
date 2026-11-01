"""Standalone averaging, competing regular outcomes and a sampled band test."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from analyze_soft_constellations import mean_only


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True, exist_ok=False); rows=[]
    for name in [Path(__file__).name, 'analyze_soft_constellations.py', 'search_soft_constellations.py', 'search_raar_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for ident in [155, 196, 251, 389]:
        source=args.source/f'case-{ident}.npz'; cohort=args.cohort/f'case-{ident}.npz'
        with np.load(source) as f: points=f['points']; beta,tau=f['params']; z=f['states']
        with np.load(cohort) as f: starts=f['starts']; tails=f['tails']; rates=f['rates']; status=f['status']
        baseline=[mean_only(x.copy(),points,tau) for x in starts]; outcomes=[]
        for k, tail in enumerate(tails):
            assert status[k,0]==0
            scale=max(float(np.linalg.norm(np.ptp(tail,axis=0))),1.)
            residuals=[float(np.max(np.linalg.norm(tail[-2048+lag:]-tail[-2048:-lag],axis=1))) for lag in range(1,513)]
            periods=[j+1 for j,v in enumerate(residuals) if v<1e-8*scale]
            outcomes.append(dict(start_index=k, detected_period=min(periods) if periods else None, minimum_window=float(rates[k].min()), mean_rate=float(rates[k].mean()), return_residuals=residuals))
        row=dict(id=ident, source_sha256=sha(source), cohort_sha256=sha(cohort), standalone_mean_control=[dict(initial=x.tolist(), last=y.tolist(), iterations=int(n), residual=float(r)) for x,(y,n,r) in zip(starts,baseline)], cohort_outcomes=outcomes)
        if ident==196:
            center=z.mean(0); _,v=np.linalg.eigh(np.cov(z.T)); direction=v[:,-1]; t=(z-center)@direction; ordered=np.sort(t)
            indices=np.argsort(np.diff(ordered))[-3:]; cuts=np.sort((ordered[indices]+ordered[indices+1])/2)
            dense=args.dense/f'case-{ident}-states.npz'
            with np.load(dense) as f: actual=f['states']
            scores=(actual-center)@direction; labels=np.searchsorted(cuts,scores)
            transitions=np.bincount(4*labels[:-1]+labels[1:],minlength=16).reshape(4,4)
            row['four_band_test']=dict(method='Three largest gaps along first PCA direction of the 100000-state terminal long sample; thresholds tested on the separate two-million-state rendering run.', center=center.tolist(), direction=direction.tolist(), cuts=cuts.tolist(), source_sha256=sha(dense), counts=np.bincount(labels,minlength=4).tolist(), transitions=transitions.tolist(), score_ranges=[[float(scores[labels==k].min()),float(scores[labels==k].max())] for k in range(4)], scope='A finite partition test, not an invariant decomposition, attractor certificate or new mechanism.')
        rows.append(row)
        print(json.dumps(dict(id=ident, baseline_max_iterations=max(r[1] for r in baseline), periods={str(p):sum(o['detected_period']==p for o in outcomes) for p in sorted(set(o['detected_period'] for o in outcomes),key=lambda p:-1 if p is None else p)}, band_test=row.get('four_band_test'))),flush=True)
    report=dict(cases=rows, scope='Numerical stopping residuals, finite period detections and sampled bands only. Known standalone mean-shift map.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','cohort','dense','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
