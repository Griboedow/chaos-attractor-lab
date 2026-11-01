"""Classify all images of an affine window as regular or Whitney-fold capable.

Every accepted cell either has a nonzero determinant throughout its final
image enclosure, or has a nonzero 2x2 minor and nonzero determinant-gradient
contraction with the corresponding row-cross-product kernel candidate.
The latter certifies the fold criterion at EVERY critical point in the box;
it does not by itself prove that a critical point exists.
"""
import argparse
import gzip
import json
import time
from pathlib import Path
import mpmath as mp
import numpy as np
from numba import njit
from certify_raar_critical_doubleton import image
from certify_raar_coverings import Model
from certify_raar_cycles import exactfloat
from certify_raar_fold import jet, cross, sha
from search_raar_circles import step


@njit(cache=True)
def weights(x,p,S,steps):
    M=S.copy()
    for _ in range(steps):
        x,J,axis,flag=step(x,p)
        if flag:return np.ones(3)
        M=J@M
    return np.sqrt((M*M).sum(axis=0))


def strictly_nonzero(v):
    return bool(v.a>0 or v.b<0)


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    sources=['certify_raar_critical_fold_region.py','certify_raar_critical_doubleton.py','certify_raar_critical_crossing.py','certify_raar_fold.py','certify_raar_coverings.py','certify_raar_cycles.py','certify_raar_multiplier_events.py','search_raar_circles.py']
    for name in sources:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (args.out/'input.json').write_bytes(args.input.read_bytes())
    mp.mp.dps=args.digits;mp.iv.dps=args.digits;iv=mp.iv
    raw=json.loads(args.input.read_text());model=Model(raw['params']);p=[exactfloat(v,iv) for v in raw['params']]
    c,B,r=raw['centers'][args.window],raw['bases'][args.window],raw['radii'][args.window]
    center=[exactfloat(v,iv) for v in c]
    S=[[exactfloat(B[i][j],iv)*exactfloat(r[j],iv) for j in range(3)] for i in range(3)]
    numerical_S=np.array(B)*np.array(r);parameters=np.array(raw['params'])
    stack=[([0.,0.,0.],[1.,1.,1.],'')];examined=accepted=fold_boxes=depth=0
    begin=last=time.time();minimum=iv.mpf('inf')
    with gzip.open(args.out/'volume.jsonl.gz','wt',encoding='utf-8') as out:
        while stack and examined<args.max_cells:
            mid,rad,path=stack.pop();examined+=1;depth=max(depth,len(path)//2);good=False;record={}
            try:
                X,D,axis=image(model,center,S,mid,rad,args.steps)
                if strictly_nonzero(D):
                    good=True;record={'classification':'regular','determinant':str(D)}
                else:
                    D2,gradient,J,H,axes=jet(X,p,iv)
                    assert all(v.a>0 for v in axes)
                    for a,b in [(0,1),(0,2),(1,2)]:
                        kernel=cross(J[a],J[b]);components=[k for k in range(3) if strictly_nonzero(kernel[k])]
                        contraction=sum(gradient[k]*kernel[k] for k in range(3))
                        if components and strictly_nonzero(contraction):
                            good=True;fold_boxes+=1
                            record={'classification':'fold_if_critical','determinant':str(D),'second_determinant':str(D2),'rows':[a,b],'minor_component':components[0],'minor':str(kernel[components[0]]),'kernel':[str(v) for v in kernel],'gradient':[str(v) for v in gradient],'transversality':str(contraction)}
                            break
            except (ArithmeticError,ValueError,ZeroDivisionError,OverflowError):pass
            if good:
                accepted+=1;minimum=min(minimum,axis)
                record.update(path=path,image=[str(v) for v in X],minimum_axis=str(axis))
                out.write(json.dumps(record,separators=(',',':'))+'\n')
            else:
                sensitivity=weights(np.array(c)+numerical_S@np.array(mid),parameters,numerical_S,args.steps)
                assert np.isfinite(sensitivity).all()
                split=int(np.argmax((sensitivity+1e-12)*np.array(rad)))
                assert rad[split]>2**-40,('subdivision_limit',path)
                cr=rad.copy();cr[split]*=.5
                for sign in [-1,1]:
                    cm=mid.copy();cm[split]+=sign*cr[split]
                    stack.append((cm,cr,path+str(split)+('0' if sign<0 else '1')))
            if time.time()-last>20:
                print(json.dumps({'examined':examined,'accepted':accepted,'fold_boxes':fold_boxes,'pending':len(stack),'depth':depth,'seconds':time.time()-begin}),flush=True);last=time.time()
    if stack:(args.out/'pending.json').write_text(json.dumps(stack,indent=2)+'\n')
    report={'status':'passed' if not stack else 'incomplete','window':args.window,'steps':args.steps,'interval_digits':args.digits,'mpmath':mp.__version__,'numpy':np.__version__,'examined':examined,'accepted':accepted,'fold_boxes':fold_boxes,'regular_boxes':accepted-fold_boxes,'pending':len(stack),'maximum_depth':depth,'minimum_axis':str(minimum),'certificate_sha256':sha(args.out/'volume.jsonl.gz'),'input_sha256':sha(args.input),'sources':{name:sha(args.out/name) for name in sources},'seconds':time.time()-begin,'scope':'All critical points in the stated forward window image satisfy the standard smooth-fold criterion, conditional on complete coverage. This classification does not itself establish existence or invariant-set membership.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',type=Path,default=Path('results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'))
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--window',type=int,default=0)
    p.add_argument('--steps',type=int,default=42)
    p.add_argument('--digits',type=int,default=70)
    p.add_argument('--max-cells',type=int,default=20000)
    main(p.parse_args())
