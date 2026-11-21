"""Certify stable cycles selected from the one-parameter diagonal of the entropy box."""
import argparse
import csv
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, ctx
from scipy.linalg import solve_discrete_lyapunov

from confirm import native_spectrum
from continuation import parameters as scan_parameters
from search import save
from sink import step


def parameters(t):
    one = t*0+1
    a, b = (34+t)/1000,(24+t)/1000
    return [(448+4*t)/100,a,b,one-a-b,(45+10*t)/100000]


def advance(x,p,n):
    j=mp.eye(3)
    for _ in range(n):
        x,d=step(x,p,mp.exp)
        j=mp.matrix(d)*j
    return mp.matrix(x),j


def bases(j):
    result=[np.eye(3)]
    vals,vec=np.linalg.eig(j)
    cols=[]
    for e,v in zip(vals,vec.T):
        if abs(e.imag)<1e-16:
            cols.append(v.real)
        elif e.imag>0:
            cols.extend([v.real,v.imag])
    if len(cols)==3:
        result.append(np.column_stack(cols))
    try:
        energy=solve_discrete_lyapunov(j.T,np.eye(3))
        result.append(np.linalg.inv(np.linalg.cholesky(energy).T))
    except np.linalg.LinAlgError:
        pass
    for v in result:
        if np.linalg.cond(v)<1e10:
            d=np.linalg.solve(v,j@v)
            yield np.linalg.norm(d,ord=2),v


def propose(index,den,period,seed=None):
    mp.mp.dps=period+115
    par=parameters(mp.mpf(index)/den)
    p=list(map(float,par))
    if seed is None:
        scan_par=scan_parameters(np.linspace(0,1,den+1)[index])
        _,_,tail=native_spectrum(scan_par,np.array([.59,1.12,.47]),40000,6,9000,3)
        seed=tail[-1]
    x=seed
    x=mp.matrix([mp.mpf(float(k)) for k in x])
    for _ in range(10):
        z,j=advance(x,par,period)
        x-=mp.lu_solve(j-mp.eye(3),z-x)
    z,j=advance(x,par,period)
    assert max(abs(k) for k in (z-x))<mp.power(10,-period-60)
    xbest=None
    for phase in range(period):
        # The derivative is only a proposal. The ball verifier recomputes it from scratch.
        if phase % max(1,period//65):
            x=mp.matrix(step(x,par,mp.exp)[0])
            continue
        q=list(map(float,x))
        deriv=np.eye(3)
        for _ in range(period):
            q,d=step(q,p,math.exp)
            deriv=np.asarray(d)@deriv
        for rate,v in bases(deriv):
            if xbest is None or rate<xbest[0]:
                xbest=(rate,list(x),v)
        x=mp.matrix(step(x,par,mp.exp)[0])
    rate,x0,v=xbest
    print("proposal",period,rate,flush=True)
    assert rate<1
    return {"t_numerator":index,"t_denominator":den,"period":period,
            "center":[mp.nstr(q,period+85) for q in x0],
            "basis":[[format(q,".17g") for q in row] for row in v],
            "radius":"1e-"+str(period+32)}


def verify(row):
    period=row["period"]
    ctx.prec=4*(period+118)
    p=parameters(arb(row["t_numerator"])/row["t_denominator"])
    v=arb_mat([[arb(q) for q in line] for line in row["basis"]])
    inv=v.inv()
    center=arb_mat([[arb(q)] for q in row["center"]])
    eps=arb(row["radius"])
    ball=arb(0,row["radius"])
    cube=center+v*arb_mat([[ball],[ball],[ball]])
    ident=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    ds=ident
    absolute=[cube[i,0] for i in range(3)]
    at_center=[center[i,0] for i in range(3)]
    assert all(arb("1e-22")<q and q<10 for q in absolute)
    beyond=any(q>arb("2.85") for q in absolute)
    divisors=[]
    for k in range(1,period+1):
        absolute,j=step(absolute,p,lambda q:q.exp())
        at_center,_=step(at_center,p,lambda q:q.exp())
        ds=arb_mat(j)*ds
        if k<period and period%k==0:
            assert any(not q.overlaps(cube[i,0]) for i,q in enumerate(absolute)),k
            divisors.append(k)
        beyond=beyond or any(q>arb("2.85") for q in absolute)
    df=inv*ds*v
    # Gershgorin bounds the largest eigenvalue of D^T D. A Euclidean ball in
    # normalized coordinates is contained in the cube used for the derivatives.
    # The fixed-point ball and its center are both inside the global trapping cube.
    square=df.transpose()*df
    upper=max(square[i,i]+sum(abs(square[i,j]) for j in range(3) if j!=i) for i in range(3))
    residual=inv*(arb_mat([[q] for q in at_center])-center)
    closure=sum(residual[i,0]**2 for i in range(3)).sqrt()/eps
    if not upper.sqrt()+closure<1:
        # Gershgorin can lose sharpness off the eigenbasis. Certify an explicit
        # spectral norm through positive definiteness of h^2 I-D^T D instead.
        for bound in (".75",".9",".95",".99",".999",".9999",".99999",".999999",".9999999"):
            h=arb(bound)
            m=h*h*ident-square
            if m[0,0]>0 and m[0,0]*m[1,1]-m[0,1]*m[1,0]>0 and m.det()>0 and h+closure<1:
                upper=h*h
                break
    assert upper.sqrt()+closure<1,(period,str(upper),str(closure))
    return {"precision_bits":ctx.prec,"uniform_derivative_norm_bound":str(upper.sqrt()),
            "center_residual_over_radius":str(closure),"separated_divisors":divisors,
            "one_phase_above_symbolic_box":bool(beyond)}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--verify",action="store_true")
    ap.add_argument("--zoom",type=Path)
    args=ap.parse_args()
    rows=json.loads(args.out.read_text()) if args.verify else []
    if args.verify:
        for row in rows:
            verify(row)
    else:
        proposals=[(ix,1200,p,None) for ix,p in [(219,7),(192,21),(1187,78),(1182,156)]]
        if args.zoom:
            proposals=[]
            with args.zoom.open(newline="",encoding="utf-8") as f:
                grid=list(csv.DictReader(f))
            for p in sorted({int(v["period_multiple_of_39"]) for v in grid if int(v["period_multiple_of_39"])>156}):
                candidates=[]
                for v in grid:
                    if v["direction"]!="forward" or int(v["period_multiple_of_39"])!=p:
                        continue
                    ix,den=int(v["numerator"]),int(v["denominator"])
                    x=[float(v[k]) for k in "xyz"]
                    par=scan_parameters(ix/den)
                    end=x
                    d=np.eye(3)
                    for _ in range(p):
                        end,j=step(end,par,math.exp)
                        d=np.asarray(j)@d
                    candidates.append((max(abs(np.linalg.eigvals(d))),ix,den,x))
                radius,ix,den,x=min(candidates)
                assert radius<1,(p,radius)
                proposals.append((ix,den,p,x))
        for ix,den,p,seed in proposals:
            row=propose(ix,den,p,seed)
            row["verified"]=verify(row)
            rows.append(row)
            save(args.out,rows)
            print("certified",p,flush=True)
