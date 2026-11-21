"""Exclude collapse of a connected trapped component by isolating a second-return orbit."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/"research/shape_search_20261121/maps"))
from flint import arb,arb_mat,ctx
import numpy as np
from scipy.optimize import root
from search import START,sample
from periodic_points import iterate
from verify_boxes import frames,rotate,column,matrix,values

def repeat(x,p,n,I):
    J=I
    for _ in range(n):x,J=rotate(values(x),p,J)
    return x,J

def verify(x,p,m,o,V):
    ctx.prec=384
    I=matrix(np.eye(3).tolist());c=column(x)
    for _ in range(5):
        y,J=repeat(c,p,2*m,I)
        c=c-(J-I).inv()*(y-c);c=column([v.mid() for v in values(c)])
    y,J=repeat(c,p,2*m,I)
    P=(J-I).inv()
    for bits in (66,92,126,162):
        rad=arb(2)**(-bits);wide=c+column([arb(0,rad)]*3)
        if not all(abs(v)<1 for v in values(V*(wide-o))):continue
        halfway,D=repeat(wide,p,m,I)
        if not any(a>b or a<b for a,b in zip(values(halfway),values(wide))):continue
        z,K=repeat(halfway,p,m,D)
        N=I-P*(K-I);d=P*(y-c)
        if all(sum(N[i,j].abs_upper() for j in range(3))*rad+d[i,0].abs_upper()<rad for i in range(3)):
            return dict(distinct_second_return_orbit=True,native_center_dyadic=[[int(k) for k in v.man_exp()] for v in values(c)],
                radius_bits=bits,precision_bits=ctx.prec,trapped_phase_box_interior=True,first_return_disjoint=True,
                row_contraction_bounds_approx=[float(sum(N[i,j].abs_upper() for j in range(3))) for i in range(3)])

def prove(box):
    ctx.prec=1792
    oo,_,vv=frames(box);o,V=oo[0],vv[0]
    p=np.array(box["params"]);pars=[arb(float(v)) for v in p]
    m=box["case_order"]
    x=column(START);entry=None
    for n in range(650):
        if all(abs(v)<1 for v in values(V*(x-o))):entry=n;break
        x,_=rotate(values(x),pars)
        if any(not v.is_finite() or float(v.rad())>1e-3 for v in values(x)):break
    _,zz=sample(p[None,:],START[None,:],24000,8192,4,8192)
    O=np.array([float(v) for v in values(o)])
    W=np.array([[float(V[i,j]) for j in range(3)] for i in range(3)])
    z=zz[0];inside=np.max(abs((z-O)@W.T),axis=1)<1
    select=np.flatnonzero(inside[:-2*m])
    # Rank close recurrences, also seed every transverse part of the observed strand.
    score=np.linalg.norm(z[select]-z[select+2*m],axis=1)/(np.linalg.norm(z[select]-z[select+m],axis=1)+1e-16)
    select=np.r_[select[np.argsort(score)[:12]],select[::max(1,len(select)//40)]]
    for index in select:
        def residual(x):
            y,J=iterate(x,p,2*m)
            return y-x,J-np.eye(3)
        c=root(residual,z[index],jac=True,tol=1e-11).x
        if max(abs(W@(c-O)))>=1 or max(abs(residual(c)[0]))>1e-8:continue
        if max(abs(iterate(c,p,m)[0]-c))<1e-7:continue
        try:
            proof=verify(c,pars,m,o,V)
        except (ZeroDivisionError,ValueError):continue
        if proof:return dict(case=box["case"],period=2*m,canonical_start_entry_step=entry,**proof)
    return dict(case=box["case"],canonical_start_entry_step=entry,distinct_second_return_orbit=False)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("traps",type=Path);ap.add_argument("out",type=Path)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    for f in sorted(a.traps.glob("*.json")):
        box=json.loads(f.read_text())["box"];r=prove(box)
        (a.out/f.name).write_text(json.dumps(r,indent=2)+"\n")
        print(f.stem,r["distinct_second_return_orbit"],r["canonical_start_entry_step"],flush=True)
