"""Account for periodic positive-growth flags and test two fixed-point captures."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from search_soft_sphere import step, screen, checks
from check_soft_sphere_precision import iterate


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x): p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8',newline='\n')


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps=80
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    for name in ['check_soft_sphere_precision.py','search_soft_sphere.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:
        params=f['params'].copy();points=f['points'].copy();counts=f['counts'].copy()
        values=screen(params,points,counts,5000,15000)
        comparisons={key:bool(np.array_equal(f[key],v,equal_nan=True)) for key,v in zip(['rates','tails','status','minima','soft_stats','boundratios'],values)}
    assert all(comparisons.values())
    followed=json.loads((args.follow/'report.json').read_text(encoding='utf-8'))['cases']
    periodic=[v for v in followed if v.get('detected_period') is not None]
    positive=[v for v in followed if v['minimum_leading_window']>.005]
    both=[v for v in periodic if v['minimum_leading_window']>.005]
    rows=[]
    for ident in [35,300]:
        source=args.follow/f'case-{ident}.npz'
        with np.load(source) as f:
            x=f['states'][-1].copy();z=f['states'][-2048:].copy()
        beta,tau=[mp.mpf(float(v)) for v in params[ident]]
        p=[[mp.mpf(float(v)) for v in vtx] for vtx in points[ident,:counts[ident]]]
        y,J,r,*_=step(x,float(beta),points[ident,:counts[ident]],float(tau))
        eig,vec=np.linalg.eig(J);i=int(np.argmax(abs(eig)));direction=vec[:,i].real;direction/=np.linalg.norm(direction)
        baseline=[mp.mpf(float(v)) for v in x]
        perturbed=[v+mp.mpf('1e-40')*mp.mpf(float(w)) for v,w in zip(baseline,direction)]
        separation=[];first=None
        for k in range(257):
            distance=mp.sqrt(sum((a-b)**2 for a,b in zip(baseline,perturbed)))
            if first is None and distance>mp.mpf('.01'):first=k
            separation.append(str(distance))
            baseline,_=iterate(baseline,beta,tau,p,mp);perturbed,_=iterate(perturbed,beta,tau,p,mp)
        assert first is not None
        row=dict(id=ident,source_sha256=sha(source),terminal_state=x.tolist(),terminal_state_hex=[float(v).hex() for v in x],terminal_range=np.ptp(z,axis=0).tolist(),one_step_binary64_residual=float(np.max(abs(y-x))),jacobian=J.tolist(),sampled_multipliers=[float(v.real) for v in eig],reflection_radius=float(r),perturbation='1e-40 times the numerically computed dominant real eigenvector',direction=direction.tolist(),digits=80,first_separation_above_0_01=first,separations=separation,scope='Independent finite-time instability of a numerically captured fixed point; neither a root certificate nor an asymptotic basin claim.')
        rows.append(row);save(args.out/f'case-{ident}.json',row)
        print(json.dumps(dict(id=ident,first_separation_above_0_01=first,multipliers=row['sampled_multipliers'])),flush=True)
    report=dict(source_sha256=sha(args.source),follow_report_sha256=sha(args.follow/'report.json'),screen_bitwise_replay=comparisons,formula_checks=checks(params,points,counts),followed=len(followed),positive_windows=len(positive),finite_period_flags=len(periodic),periodic_and_positive_ids=[v['id'] for v in both],eligible_nonperiodic_positive_ids=[v['id'] for v in positive if v.get('detected_period') is None],capture_controls=rows,scope='A positive tangent exponent does not turn a repeating numerical orbit into a chaotic attractor. All period flags are rejected from visual promotion regardless of tangent growth.')
    save(args.out/'report.json',report)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--follow',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mpmath-path',type=Path,required=True);main(p.parse_args())
