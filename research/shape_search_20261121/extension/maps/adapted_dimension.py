"""Bound return secants in phase-adapted bases and a weighted maximum norm."""
import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/"research/shape_search_20261121/maps"))
import numpy as np
from flint import arb,arb_mat,ctx
from scipy.spatial.transform import Rotation
from slice_geometry import affine_frame,cubes,inverse_jac,values
from verify_boxes import matrix
from slice_bound import ball

OLD=ROOT/"results/shape_search_20261121/maps/certificates"

def domains():
    d=json.loads((OLD/"2503155-dimension-phase-directed.json").read_text())["native_phase_hulls"]
    a=json.loads((OLD/"2503155-dimension-phase-arb.json").read_text())["phase_enclosures_binary64"]
    return [[[max(x[0],y[0]),min(x[1],y[1])] for x,y in zip(A,B)] for A,B in zip(a,d)]

def local(domain,p,n):
    mid,r=affine_frame(domain);I=matrix(np.eye(3).tolist())
    for x,h in cubes(n):
        yield inverse_jac([a+b*arb(u,v) for a,b,u,v in zip(mid,r,x,h)],p,I)

def approximate(J,F):
    result=np.eye(3)
    for i in (0,4,3,2,1):
        lo,delta=J[i];v=F[(i-1)%5].T;f=F[i]
        a=v@lo@f;b=abs(v)@delta@abs(f)
        low=(a-b).min(axis=0);high=(a+b).max(axis=0)
        result=np.maximum(abs(low),abs(high))@result
    return max(abs(np.linalg.eigvals(result)))

def basis(bounds,p):
    Js=[];F=[np.eye(3) for _ in bounds]
    for b in bounds:
        block=list(local(b,p,2))
        c=np.array([[[float(A[i,j].mid()) for j in range(3)] for i in range(3)] for A in block])
        r=np.array([[[float(A[i,j].rad()) for j in range(3)] for i in range(3)] for A in block])
        Js.append((c,r))
    for _ in range(15):
        for i in (0,4,3,2,1):
            F[(i-1)%5]=np.linalg.qr(Js[i][0].mean(axis=0)@F[i])[0]
    F=np.array(F);best=approximate(Js,F)
    rng=np.random.default_rng(6247508)
    for k in range(3100):
        trial=F.copy();i=rng.integers(5)
        trial[i]=trial[i]@Rotation.from_rotvec(rng.normal(size=3)*(.16*(1-k/3100)+.002)).as_matrix()
        size=approximate(Js,trial)
        if size<best:best,F=size,trial
    return F.tolist()

def verify(frames,grid):
    ctx.prec=150
    p=[arb(float(v),2.**-24) for v in json.loads((OLD/"2503155-box.json").read_text())["params"]]
    bounds=domains();F=[matrix(f) for f in frames]
    enclosures=[]
    for i,b in enumerate(bounds):
        H=[[None]*3 for _ in range(3)]
        V=F[(i-1)%5].inv()
        for A in local(b,p,grid):
            C=V*A*F[i]
            for j in range(3):
                for k in range(3):
                    H[j][k]=C[j,k] if H[j][k] is None else H[j][k].union(C[j,k])
        enclosures.append(arb_mat(H))
    R=matrix(np.eye(3).tolist())
    for i in (0,4,3,2,1):R=enclosures[i]*R
    B=arb_mat([[R[i,j].abs_upper() for j in range(3)] for i in range(3)])
    # Positive weights are only a proposal. The final Collatz bound is Arb.
    vals,vec=np.linalg.eig(np.array([[float(B[i,j]) for j in range(3)] for i in range(3)]))
    w=np.abs(vec[:,np.argmax(vals.real)].real);w=(w/max(w)).tolist()
    if min(w)<=0:raise ValueError("Weights must be strictly positive")
    L=max((sum(B[i,j]*arb(w[j]) for j in range(3))/arb(w[i])).abs_upper() for i in range(3))
    s=1+arb(2).log()/L.log();q=Fraction(int(float(s.lower())*1e6),1000000)
    if not(s>ball(q) and L>1):raise ArithmeticError("Rounding or expansion bound failed")
    return dict(case=2503155,block_radius_bits=24,dyadic_cells_per_phase=2**(3*grid),precision_bits=ctx.prec,
        frames=frames,return_frame_weights=w,phase_hulls=bounds,dimension_strictly_greater_than=str(q),
        decimal_lower=float(q),inverse_return_norm_upper_rounded_up=int(float(L.upper()))+1,
        source_certificate_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest()
        for f in OLD.glob("2503155-dimension-phase-*.json")})

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("output",type=Path)
    ap.add_argument("--grid",type=int,default=3);ap.add_argument("--frames",type=Path)
    a=ap.parse_args();ctx.prec=150
    p=[arb(float(v),2.**-24) for v in json.loads((OLD/"2503155-box.json").read_text())["params"]]
    F=json.loads(a.frames.read_text())["frames"] if a.frames else basis(domains(),p)
    result=verify(F,a.grid);a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result["decimal_lower"],result["inverse_return_norm_upper_rounded_up"],flush=True)
