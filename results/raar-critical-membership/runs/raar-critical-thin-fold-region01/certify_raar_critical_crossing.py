"""Enclose det DF after a declared number of iterates on both exit faces.

The input window is an exact affine image of [-1,1]^3. This program checks
face signs and regularity only. A separate topological argument and existing
self-covering/trapping certificates are needed for invariant-set membership.
"""
import argparse
import gzip
import hashlib
import json
import time
from pathlib import Path
import mpmath as mp
from certify_raar_coverings import Model
from certify_raar_cycles import exactfloat, matmul, matvec, det


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def image(model, center, matrix, mid, rad, steps):
    iv=mp.iv
    q=[exactfloat(v,iv) for v in mid]
    delta=[iv.mpf([-exactfloat(v,iv).b,exactfloat(v,iv).b]) for v in rad]
    offset=matvec(matrix,q)
    C=[center[i]+offset[i] for i in range(3)]
    T=matrix
    minimum=iv.mpf('inf')
    for k in range(steps+1):
        offset=matvec(T,delta)
        X=[C[i]+offset[i] for i in range(3)]
        _,J,axis=model.step(X)
        minimum=min(minimum,axis)
        if k==steps:
            return X,det(J),minimum
        C,_,_=model.step(C,False)
        T=matmul(J,T)


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    for name in ['certify_raar_critical_crossing.py','certify_raar_coverings.py','certify_raar_cycles.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (args.out/'input.json').write_bytes(args.input.read_bytes())
    mp.mp.dps=args.digits;mp.iv.dps=args.digits
    iv=mp.iv;raw=json.loads(args.input.read_text());model=Model(raw['params'])
    c,B,r=raw['centers'][args.window],raw['bases'][args.window],raw['radii'][args.window]
    center=[exactfloat(v,iv) for v in c]
    matrix=[[exactfloat(B[i][j],iv)*exactfloat(r[j],iv) for j in range(3)] for i in range(3)]
    scales=[sum(abs(matrix[i][j]).b for i in range(3)) for j in range(3)]
    rows=[];begin=time.time()
    faces=[-1,1] if args.face==0 else [args.face]
    for face in faces:
        sign=face*args.orientation
        stack=[([float(face),0.,0.],[0.,1.,1.],'')]
        accepted=examined=depth=0;minimum=iv.mpf('inf');dmin=iv.mpf('inf');last=time.time()
        with gzip.open(args.out/f'face-{face}.jsonl.gz','wt',encoding='utf-8') as out:
            while stack and examined<args.max_cells:
                mid,rad,path=stack.pop();examined+=1;depth=max(depth,len(path)//2)
                try:
                    X,D,axis=image(model,center,matrix,mid,rad,args.steps)
                    good=(sign*D).a>0
                except (ArithmeticError,ValueError,ZeroDivisionError,OverflowError):
                    good=False
                if good:
                    accepted+=1;minimum=min(minimum,axis);dmin=min(dmin,(sign*D).a)
                    out.write(json.dumps({'path':path,'image':[str(v) for v in X],'determinant':str(D),'minimum_axis':str(axis)},separators=(',',':'))+'\n')
                else:
                    split=max([1,2],key=lambda j:float(scales[j])*rad[j])
                    assert rad[split]>2**-40,('subdivision_limit',path)
                    cr=rad.copy();cr[split]*=.5
                    for s in [-1,1]:
                        cm=mid.copy();cm[split]+=s*cr[split]
                        stack.append((cm,cr,path+str(split)+('0' if s<0 else '1')))
                if time.time()-last>20:
                    print(json.dumps({'face':face,'examined':examined,'accepted':accepted,'pending':len(stack),'depth':depth,'seconds':time.time()-begin}),flush=True);last=time.time()
        row={'face':face,'sign':sign,'status':'passed' if not stack else 'incomplete','examined':examined,'accepted':accepted,'pending':len(stack),'maximum_depth':depth,'minimum_axis':str(minimum),'absolute_determinant_lower':str(dmin),'certificate_sha256':sha(args.out/f'face-{face}.jsonl.gz')}
        if stack:(args.out/f'face-{face}-pending.json').write_text(json.dumps(stack,indent=2)+'\n')
        rows.append(row);print(json.dumps(row),flush=True)
        (args.out/'progress.json').write_text(json.dumps(rows,indent=2)+'\n')
    report={'status':'passed' if len(rows)==2 and all(v['status']=='passed' for v in rows) else 'incomplete','window':args.window,'steps':args.steps,'orientation':args.orientation,'interval_digits':args.digits,'mpmath':mp.__version__,'input_sha256':sha(args.input),'source_sha256':sha(Path(__file__)),'checks':rows,'seconds':time.time()-begin,'scope':'Strict opposite determinant signs on images of the stated exit faces; invariant-set conclusion requires a separately checked crossing-continuum argument and trapping/self-covering hypotheses.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',type=Path,default=Path('results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'))
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--window',type=int,default=0)
    p.add_argument('--steps',type=int,default=42)
    p.add_argument('--orientation',type=int,choices=[-1,1],default=1)
    p.add_argument('--face',type=int,choices=[-1,0,1],default=0)
    p.add_argument('--digits',type=int,default=70)
    p.add_argument('--max-cells',type=int,default=50000)
    main(p.parse_args())
