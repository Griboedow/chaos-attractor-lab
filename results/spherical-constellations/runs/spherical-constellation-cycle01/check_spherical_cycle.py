"""Exact rational checks for a circle of period-two points in sphere RAAR.

The analytic reduction is documented in docs/spherical-constellation-search.md.
This does not certify that a particular binary64 trajectory enters its basin.
"""
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path
import numpy as np


def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def exact(x): return F(*float(x).as_integer_ratio())
def encoded(x): return dict(numerator=str(x.numerator),denominator=str(x.denominator))


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source) as f:
        beta=exact(f['params'][0]); points=[[exact(x) for x in p] for p in f['points']]
        z=f['states']; branches=f['branches']
    b=points[args.site]; s=dot(b,b); r=beta/(1+beta); a2=r*r-s/4
    assert 0<beta<1 and 0<s<4*r*r and a2>0
    # Min over x=3b/2+w, w perpendicular to b, |w|^2=a2:
    # |x-c|^2-|x-b|^2 = A-2*w.(c-b) >= A-2*sqrt(a2*q).
    # The two rational inequalities below prove this exceeds the target.
    threshold=F(245,1000); rows=[]
    for i,c in enumerate(points):
        if i==args.site: continue
        d=[x-y for x,y in zip(c,b)]; bd=dot(b,d)
        A=dot(c,c)-s-3*bd; q=dot(d,d)-bd*bd/s
        assert q>=0
        first=A-threshold; second=first*first-4*a2*q
        assert first>0 and second>0
        rows.append(dict(site=i,A=encoded(A),projected_norm_squared=encoded(q),positive_first=encoded(first),positive_squared_gap=encoded(second)))
    gamma=s/(2*r*r)-1; tau=(1+beta)**2*gamma**2-2*beta; delta=beta*beta
    # The non-tangent eigenvalues of DF^2 have polynomial t^2-tau*t+delta.
    schur=[1-delta,1-tau+delta,1+tau+delta]
    assert all(v>0 for v in schur)
    discriminant=tau*tau-4*delta; assert discriminant<0
    last=int(branches[-1]); earlier=np.flatnonzero(branches!=last)
    float_b=np.array([float(x) for x in b]); terminal=z[-2048:]
    residual=float(np.max(np.linalg.norm(terminal[2:]-terminal[:-2],axis=1)))
    samples=[]
    for k in range(1,513): samples.append(float(np.max(np.linalg.norm(terminal[k:]-terminal[:-k],axis=1))))
    report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),site=args.site,beta=encoded(beta),point=[encoded(v) for v in b],norm_squared=encoded(s),sphere_projection_radius=encoded(r),cycle_circle_radius_squared=encoded(a2),uniform_squared_distance_margin=encoded(threshold),site_checks=rows,gamma=encoded(gamma),normal_trace=encoded(tau),normal_determinant=encoded(delta),schur_positive=[encoded(v) for v in schur],negative_discriminant=encoded(discriminant),numerical_only=dict(last_site=last,last_single_site_suffix_start_after_warmup=int(earlier[-1]+1),warmup=100000,period_two_residual=residual,period_residuals=samples,pair_sum_max_error=float(np.max(abs(terminal[1:]+terminal[:-1]-3*float_b))),circle_equation_max_error=float(np.max(abs(np.linalg.norm(terminal-2*float_b,axis=1)-float(r))))),scope='Exact arithmetic verifies an admissible normally attracting circle of period-two points for the exact stored inputs. Recorded trajectory convergence remains numerical; historical novelty is not asserted.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps(dict(status='passed',site=args.site,uniform_squared_distance_margin=float(threshold),cycle_circle_radius=float(a2)**.5,normal_two_step_modulus=float(beta),numerical_period_two_residual=residual)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--site',type=int,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
