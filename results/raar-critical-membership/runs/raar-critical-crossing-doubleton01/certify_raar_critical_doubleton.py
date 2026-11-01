"""Reconditioned mean-value enclosure for the critical-crossing face audit.

Maintain X subset C + A*delta + B*r. A and B are arbitrary selected point
matrices; their numerical selection need not be accurate. At each step,

 F(X) subset F(C) + Anew*delta
          + (J(X)*A-Anew)*delta + J(X)*B*r.

The last two terms are enclosed after applying an interval inverse of a
selected new B. This is a mean-value inclusion, not a sampled error estimate.
The original subdivision driver is shared with the unreconditioned audit.
"""
import argparse
import json
from pathlib import Path
import mpmath as mp
import numpy as np
import certify_raar_critical_crossing as driver
from certify_raar_coverings import inverse
from certify_raar_cycles import exactfloat, matmul, matvec, det, eye


def point_matrix(A):
    assert np.isfinite(A).all()
    return [[exactfloat(float(v),mp.iv) for v in row] for row in A]


def midpoint(A):
    return np.array([[float((v.a+v.b)/2) for v in row] for row in A])


def image(model, center, matrix, mid, rad, steps):
    iv=mp.iv
    q=[exactfloat(v,iv) for v in mid]
    delta=[iv.mpf([-exactfloat(v,iv).b,exactfloat(v,iv).b]) for v in rad]
    off=matvec(matrix,q);C=[center[i]+off[i] for i in range(3)]
    A=matrix;B=eye(iv);r=[iv.mpf(0)]*3;minimum=iv.mpf('inf')
    for k in range(steps+1):
        ad,br=matvec(A,delta),matvec(B,r)
        X=[C[i]+ad[i]+br[i] for i in range(3)]
        _,J,axis=model.step(X);minimum=min(minimum,axis)
        if k==steps:return X,det(J),minimum
        newC,Jc,_=model.step(C)
        An=point_matrix(midpoint(Jc)@midpoint(A))
        Q,_=np.linalg.qr(midpoint(Jc)@midpoint(B))
        Bn=point_matrix(Q);Bi=inverse(Bn)
        JA,JB=matmul(J,A),matmul(J,B)
        E=[[JA[i][j]-An[i][j] for j in range(3)] for i in range(3)]
        left,right=matvec(matmul(Bi,E),delta),matvec(matmul(Bi,JB),r)
        r=[left[i]+right[i] for i in range(3)]
        C,A,B=newC,An,Bn


def main(args):
    driver.image=image
    driver.main(args)
    path=args.out/'report.json';report=json.loads(path.read_text())
    report.update(enclosure='Reconditioned mean-value doubleton; arbitrary finite point matrices, interval inverse verified at every step.',enclosure_source_sha256=driver.sha(Path(__file__)),numpy=np.__version__)
    path.write_text(json.dumps(report,indent=2)+'\n')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps({'enclosure':report['enclosure'],'status':report['status'],'seconds':report['seconds']}),flush=True)


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
