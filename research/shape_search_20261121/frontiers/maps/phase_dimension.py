"""Certify a phase-restricted inverse norm for the Hausdorff slice lower bound."""
import argparse
import base64
import json
import zlib
from fractions import Fraction
from pathlib import Path
from affine_arb import test
from arb_entropy import arb_parameters
from entropy import INPUT,regions
from flint import arb,arb_mat,ctx
import numpy as np
from adapted_dimension import local
from directed_slices import derivative_hull
from directed_return import trig_constants
from entropy import parameters
from slice_bound import ball
from slice_geometry import endpoints
from verify_boxes import frames,matrix

def phase_hulls(box,record,bits):
    ctx.prec=144
    if record["bits"]!=bits or record["region"]!=regions()[-1]:raise ValueError("Expected entire return-box tree")
    t=regions()[-1];o,M,V=frames(box);o,M,V=o[0],M[0],V[0];m=box["case_order"]
    # The normalized source domain of this final region equals the full box.
    p=arb_parameters(box,bits);pc=arb_parameters(box,0)
    H=[[None]*3 for _ in range(m)]
    tree=zlib.decompress(base64.b64decode(record["partition_zlib_base64"]))
    stack=[([arb(0)]*3,[arb(1)]*3)];idx=0
    while stack:
        mid,rad=stack.pop();code=tree[idx];idx+=1
        if code==3:
            ok,_,_=test(mid,rad,M,o,o,V,p,pc,t,m,H)
            if not ok:raise RuntimeError("An archived trapping leaf failed")
            continue
        if code not in (0,1,2):raise ValueError("Invalid split")
        rr=rad.copy();rr[code]/=2
        for side in (-1,1):
            x=mid.copy();x[code]+=side*rr[code];stack.append((x,rr))
    if idx!=len(tree) or idx!=record["nodes"]:raise ValueError("Unconsumed tree")
    return [[list(endpoints(v)) for v in triple] for triple in H]

def from_matrix(J):
    B=arb_mat([[J[i,j].abs_upper() for j in range(3)] for i in range(3)])
    vals,v=np.linalg.eig(np.array([[float(B[i,j]) for j in range(3)] for i in range(3)]))
    w=abs(v[:,np.argmax(vals.real)].real);w=(w/max(w)).tolist()
    if min(w)<=0:raise ArithmeticError("The weighted norm needs positive weights")
    one=max((sum(B[i,j]*arb(w[j]) for j in range(3))/arb(w[i])).abs_upper() for i in range(3))
    L=one**2
    if not L>1:raise ArithmeticError("Invalid inverse bound")
    s=1+arb(2).log()/L.log()
    q=Fraction(int(float(s.lower())*10**6),10**6)
    if not s>ball(q):raise ArithmeticError("Incorrect rounding")
    return dict(phase_zero_weights=w,two_return_inverse_Lipschitz_upper_approximation=float(L.upper()),
        strict_dimension_lower_rational=str(q),dimension_lower=float(q))

def bound(box,record,bits,grid):
    bounds=phase_hulls(box,record,bits);m=len(bounds)
    p=arb_parameters(box,bits);pi=parameters(box,bits);trig=trig_constants()
    A=[];D=[]
    for b in bounds:
        H=[[None]*3 for _ in range(3)]
        for J in local(b,p,grid):
            for i in range(3):
                for j in range(3):H[i][j]=J[i,j] if H[i][j] is None else H[i][j].union(J[i,j])
        A.append(arb_mat(H))
        raw=derivative_hull(np.array(b),pi,grid,*trig)
        D.append(arb_mat([[arb(float(raw[i,j,0])).union(arb(float(raw[i,j,1]))) for j in range(3)] for i in range(3)]))
    out=[]
    for H in (A,D):
        J=matrix(np.eye(3).tolist())
        for i in (0,*reversed(range(1,m))):J=H[i]*J
        out.append(from_matrix(J))
    return dict(case=box["case"],block_radius_bits=bits,grid_cells_per_phase=2**(3*grid),
        starting_phase_count=m,full_inverse_iterations=2*m,common_Arb_verified_native_phase_hulls=bounds,
        Arb_inverse_derivative=out[0],directed_inverse_derivative=out[1],
        conservative_dimension_lower=min(x["dimension_lower"] for x in out))

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("tree",type=Path);ap.add_argument("output",type=Path)
    ap.add_argument("--source",type=Path,default=INPUT);ap.add_argument("--bits",type=int,required=True)
    ap.add_argument("--grid",type=int,default=3)
    a=ap.parse_args();result=bound(json.loads(a.source.read_text())["box"],json.loads(a.tree.read_text()),a.bits,a.grid)
    a.output.write_text(json.dumps(result,indent=2)+"\n");print(result["conservative_dimension_lower"],flush=True)
