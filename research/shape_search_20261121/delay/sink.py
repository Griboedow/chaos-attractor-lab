"""Certify an attracting period-39 orbit in the positive-entropy parameter box."""
import argparse
import json
from pathlib import Path

from flint import arb,arb_mat,ctx
import mpmath as mp
import numpy as np

from search import save

PAR=["4.52",".035",".025",".94",".00055"]


def step(x,p,exp):
    r,a,b,c,s=p
    u,v,w=x
    e=exp(r*(1-a*u-b*v-c*w))
    g=(1-s)*w*e
    nxt=[s*u+g,u,v]
    jac=[[s-r*a*g,-r*b*g,(1-s)*e*(1-r*c*w)],[1,0,0],[0,1,0]]
    return nxt,jac


def propose(point):
    mp.mp.dps=110
    p=list(map(mp.mpf,PAR))
    x=mp.matrix([mp.mpf(str(v)) for v in point])
    def advance(q):
        jac=mp.eye(3)
        for _ in range(39):
            q,j=step(q,p,mp.exp)
            jac=mp.matrix(j)*jac
        return mp.matrix(q),jac
    for _ in range(8):
        y,j=advance(x)
        x-=mp.lu_solve(j-mp.eye(3),y-x)
    y,j=advance(x)
    assert max(abs(v) for v in y-x)<mp.mpf("1e-90")
    best=None
    for phase in range(39):
        _,j=advance(x)
        e,v=np.linalg.eig(np.array(j.tolist(),dtype=float))
        if max(abs(e.imag))<1e-12 and np.linalg.cond(v.real)<1e9:
            score=np.linalg.cond(v.real)
            if best is None or score<best[0]:
                best=(score,{"parameters":PAR,"period":39,"center":[mp.nstr(vv,85) for vv in x],
                             "basis":[[format(float(vv),".17g") for vv in rr] for rr in v.real],
                             "radius":"1e-45"})
        x=mp.matrix(step(x,p,mp.exp)[0])
    assert best is not None
    return best[1]


def verify(data):
    ctx.prec=380
    period=data["period"]
    assert data["parameters"]==PAR and period==39
    p=list(map(arb,PAR))
    # The exact decimal center and basis define an affine cube. All matrix
    # inversions and bounds in this function use rigorous real-ball arithmetic.
    x0=arb_mat([[arb(v)] for v in data["center"]])
    v=arb_mat([[arb(c) for c in row] for row in data["basis"]])
    inv=v.inv()
    eps=arb(data["radius"])
    interval=arb(0,data["radius"])
    box=x0+v*arb_mat([[interval],[interval],[interval]])
    x=[box[i,0] for i in range(3)]
    z=[x0[i,0] for i in range(3)]
    jac=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    separations={}
    for k in range(1,period+1):
        x,j=step(x,p,lambda q:q.exp())
        z,_=step(z,p,lambda q:q.exp())
        jac=arb_mat(j)*jac
        if k in [1,3,13]:
            assert any(not x[i].overlaps(box[i,0]) for i in range(3)),k
            separations[k]=True
    deriv=inv*jac*v
    row_norms=[sum(abs(deriv[i,k]) for k in range(3)) for i in range(3)]
    residual=inv*(arb_mat([[t] for t in z])-x0)
    tests=[abs(residual[i,0])/eps+row_norms[i] for i in range(3)]
    assert all(t<1 for t in tests)
    assert all(x0[i,0]>arb("1e-22") and x0[i,0]<10 for i in range(3))
    assert all(box[i,0]>arb("1e-22") and box[i,0]<10 for i in range(3))
    return {"precision_bits":ctx.prec,"least_period":period,"proper_divisors_separated":separations,
            "contraction_row_bounds":[str(vv) for vv in row_norms],
            "strict_inclusion_row_bounds":[str(vv) for vv in tests]}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cloud",type=Path)
    ap.add_argument("--certificate",required=True,type=Path)
    args=ap.parse_args()
    if args.cloud:
        with np.load(args.cloud) as f:
            seed=f["points"][-1]
        data=propose(seed)
        data["verified"]=verify(data)
        save(args.certificate,data)
    else:
        data=json.loads(args.certificate.read_text(encoding="utf-8"))
        verify(data)
    print("Strict contraction and least period 39 verified")
