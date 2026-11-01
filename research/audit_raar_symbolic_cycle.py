"""Certify the located 49-cycle's itinerary and one-dimensional instability."""
import argparse,hashlib,json
from pathlib import Path
import mpmath as mp
from certify_raar_cycles import exactfloat,iterate,matvec,det
from check_raar_coverings_ad import exact_inverse

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);mp.mp.dps=200;mp.iv.dps=200;iv=mp.iv;raw=json.loads(args.windows.read_text());cycle=json.loads(args.cycle.read_text());assert cycle['period']==49 and cycle['minimal_period_certified'] and cycle['unstable_certified'];assert cycle['parameters']==raw['params']
    centers=[[exactfloat(v,iv) for v in row] for row in raw['centers']];matrices=[[[exactfloat(v,iv)*exactfloat(r[j],iv) for j,v in enumerate(row)] for row in B] for B,r in zip(raw['bases'],raw['radii'])];inverses=[exact_inverse(v) for v in matrices]
    x=[iv.mpf(v) for v in cycle['box']];p=[exactfloat(v,iv) for v in cycle['parameters']];coordinates=[]
    for k in range(7):
        q=matvec(inverses[k],[x[i]-centers[k][i] for i in range(3)]);assert all(v.a>-1 and v.b<1 for v in q);coordinates.append([str(v) for v in q]);x,_,_,_=iterate(x,p,7,iv)
    monodromy=[[iv.mpf(v) for v in row] for row in cycle['monodromy']];a=[iv.mpf(v) for v in cycle['characteristic_coefficients']]
    def polynomial(t):return ((t+a[0])*t+a[1])*t+a[2]
    lo=iv.mpf('579.0979');hi=iv.mpf('579.0981');assert polynomial(lo).b<0 and polynomial(hi).a>0;r=iv.mpf([lo.a,hi.b]);b=a[0]+r;c=-a[2]/r;margins=[1-abs(c),1+b+c,1-b+c];assert all(v.a>0 for v in margins)
    point=[mp.mpf(v) for v in cycle['center_decimal']];params=[exactfloat(v,mp.mp) for v in cycle['parameters']];_,M,orbit,axes=iterate(point,params,49,mp.mp);eigen=mp.eig(mp.matrix(M),left=False,right=False)
    result=dict(status='passed',period=49,least_period_certified=True,itinerary=list(range(7)),coordinate_boxes=coordinates,unstable_multiplier_enclosure=str(r),remaining_quadratic_schur_margins=[str(v) for v in margins],unstable_dimensions=1,stable_dimensions=2,numerical_multipliers=[mp.nstr(v,60) for v in eigen],orbit_decimal=[[mp.nstr(v,80) for v in row] for row in orbit],cycle_certificate_sha256=hashlib.sha256(args.cycle.read_bytes()).hexdigest(),windows_sha256=hashlib.sha256(args.windows.read_bytes()).hexdigest(),scope='Located least-period-49 saddle, with certified itinerary through all seven windows and exactly one unstable multiplier. Long binary64 monodromy products in the preliminary numerical orbit collection are ill-conditioned; their small eigenvalues are not stability evidence.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:result[k] for k in ['status','period','unstable_multiplier_enclosure','unstable_dimensions','stable_dimensions','numerical_multipliers']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--windows',type=Path,required=True);p.add_argument('--cycle',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
