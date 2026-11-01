"""Independent automatic-derivative Krawczyk and multiplier certificates."""
import argparse,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
import numpy as np

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def determinant(A):return sum(A[0][i]*A[1][(i+1)%3]*A[2][(i+2)%3]-A[0][i]*A[1][(i+2)%3]*A[2][(i+1)%3] for i in range(3))

def main(args):
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    import check_soft_symbolic_ad as ad
    ad.mp=mp;mp.mp.dps=200;mp.iv.dps=200;iv=mp.iv;args.out.mkdir(parents=True,exist_ok=False);begin=time.time();primary=json.loads((args.source/'report.json').read_text());assert primary['status']=='passed'
    for name,digest in primary['sources'].items():assert sha(args.source/name)==digest
    assert sha(args.trap)==primary['trap_sha256']
    with np.load(args.trap) as f:assert f['params'].tolist()==primary['params'] and f['points'].tolist()==primary['points'];origin=list(map(lambda x:Fraction(float(x)),f['origin']));step=Fraction(float(f['step']));shape=f['shape'];occupied=set(f['occupied'].tolist())
    model=ad.Model(primary);I=[[iv.mpf(int(i==j)) for j in range(3)] for i in range(3)];rows=[]
    def iterate(x,period):
        M=I;orbit=[];axes=[]
        for _ in range(period):
            orbit.append(x);x,J,axis=model.value_jacobian(x);M=ad.mm(J,M);axes.append(axis)
        return x,M,orbit,axes
    def frac(v):
        sign,mantissa,exponent,_=v._mpi_[0];assert v._mpi_[0]==v._mpi_[1];return Fraction((-1 if sign else 1)*mantissa)*Fraction(2)**exponent
    for old in primary['cases']:
        period=old['period'];center=[iv.mpf(v) for v in old['center_decimal']];rad=iv.mpf('1e-90');box=[iv.mpf([(v-rad).a,(v+rad).b]) for v in center];y,MI,orbit,axes=iterate(box,period);yc,MC,_,_=iterate(center,period)
        mid=mp.matrix([[mp.mpf(MC[i][j].mid)-int(i==j) for j in range(3)] for i in range(3)]);inverse=mid**-1;C=[[iv.mpf(mp.nstr(inverse[i,j],210)) for j in range(3)] for i in range(3)];H=[[MI[i][j]-I[i][j] for j in range(3)] for i in range(3)];product=ad.mm(C,H);R=[[I[i][j]-product[i][j] for j in range(3)] for i in range(3)];correction=ad.mv(C,[yc[i]-center[i] for i in range(3)]);remainder=ad.mv(R,[box[i]-center[i] for i in range(3)]);K=[center[i]-correction[i]+remainder[i] for i in range(3)];assert all(k.a>b.a and k.b<b.b for k,b in zip(K,box));contraction=max(sum(abs(v).b for v in row) for row in R);assert contraction<1
        # The smaller primary box lies in this independently unique-root box.
        assert all(b.a<iv.mpf(v).a and iv.mpf(v).b<b.b for b,v in zip(box,old['box']))
        exclusions=[]
        for d in range(1,period):
            if period%d:continue
            difference=[orbit[d][j]-box[j] for j in range(3)];assert any(v.a>0 or v.b<0 for v in difference);exclusions.append(dict(divisor=d,difference=[str(v) for v in difference]))
        brackets=[];unstable=0
        for estimate in old['numerical_multipliers']:
            low=iv.mpf(estimate)-iv.mpf('1e-25');high=iv.mpf(estimate)+iv.mpf('1e-25');lv=iv.mpf([low.a,low.a]);hv=iv.mpf([high.b,high.b]);left=determinant([[lv*I[i][j]-MI[i][j] for j in range(3)] for i in range(3)]);right=determinant([[hv*I[i][j]-MI[i][j] for j in range(3)] for i in range(3)]);assert (left.a>0 and right.b<0) or (left.b<0 and right.a>0);root=iv.mpf([lv.a,hv.b]);expanding=bool(root.a>1 or root.b<-1);stable=bool(root.a>-1 and root.b<1);assert expanding or stable;unstable+=expanding;brackets.append(dict(interval=str(root),at_left=str(left),at_right=str(right),expanding=expanding))
        assert all(iv.mpf(brackets[i]['interval']).b<iv.mpf(brackets[i+1]['interval']).a for i in range(2));assert unstable==old['unstable_dimensions'];membership=[]
        for phase,q in enumerate(orbit):
            lo=[int((frac(v.a)-origin[j])//step) for j,v in enumerate(q)];hi=[int((frac(v.b)-origin[j])//step) for j,v in enumerate(q)];assert all(0<=lo[j]<=hi[j]<shape[j] for j in range(3));keys=[]
            for i in range(lo[0],hi[0]+1):
                for j in range(lo[1],hi[1]+1):
                    for k in range(lo[2],hi[2]+1):
                        key=int((i*shape[1]+j)*shape[2]+k);assert key in occupied;keys.append(key)
            membership.append(dict(phase=phase,box=[str(v) for v in q],intersected_keys=keys))
        row=dict(id=old['id'],period=period,radius='1e-90',box=[str(v) for v in box],krawczyk=[str(v) for v in K],contraction_bound=str(contraction),primary_root_same=True,proper_divisor_exclusions=exclusions,monodromy=[[str(v) for v in row] for row in MI],multiplier_brackets=brackets,unstable_dimensions=int(unstable),trapping_membership=membership,minimum_radius=str(min(a.a for a in axes)));rows.append(row);(args.out/f"cycle-{old['id']}.json").write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(dict(id=old['id'],period=period,unstable=int(unstable),status='passed')),flush=True)
    sources={}
    for name in [Path(__file__).name,'check_soft_symbolic_ad.py']:
        data=Path(__file__).with_name(name).read_bytes();(args.out/name).write_bytes(data);sources[name]=hashlib.sha256(data).hexdigest()
    result=dict(status='passed',interval_digits=200,mpmath=mp.__version__,primary_report_sha256=sha(args.source/'report.json'),trap_sha256=sha(args.trap),sources=sources,cases=rows,seconds=time.time()-begin,scope='Independent direct-distance Gaussian formula, automatic derivative, larger unique-root boxes, determinant multiplier brackets and exact rational grid membership. Confirms the same three roots and their unstable indices. No hetero-chaos or transitivity claim.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--trap',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mpmath-path',type=Path);main(p.parse_args())
