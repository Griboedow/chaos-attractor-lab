"""Certify slit crossings and critical-value clearance of Fourier-driven cubic families."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import arb,ctx
from search import classify,basis

ctx.prec=112


def terms(row):
    h=row["h"]
    return [(j,row["coeffs"][j],row["coeffs"][h+j]) for j in range(1,h+1)
            if row["coeffs"][j] or row["coeffs"][h+j]]


def real_arb(theta,row):
    total=arb(row["coeffs"][0])
    for j,a,b in terms(row):
        if a: total+=a*(j*theta).cos()
        if b: total+=b*(j*theta).sin()
    return total/row["divisor"]


def drive(theta,row):
    theta=np.asarray(theta)
    real=np.full(theta.shape,float(row["coeffs"][0]))
    for n,c,s in terms(row):
        real+=c*np.cos(n*theta)+s*np.sin(n*theta)
    return real/row["divisor"]+1j*row["imag_amplitude"]*np.sin(row["k"]*theta)


def certify(row):
    k=row["k"];B=row["imag_amplitude"];parts=[]
    for j in range(2*k):
        x=real_arb(arb.pi()*j/k,row)
        if x>2: parts.append(3)
        elif x< -2: parts.append(-3)
        elif x> -2 and x<2: parts.append(0)
        else: raise ValueError(("Undecided slit",row["case"],j,str(x)))
    status,mu,crossings=classify(np.array(parts,dtype=float))
    assert status==0 and mu and mu==row["mu"],row["case"]
    theta=np.arange(32768)*2*np.pi/32768
    b=drive(theta,row)
    approximate=min(np.min(abs(b-2)),np.min(abs(b+2)))
    units=max(1,int(approximate*.85*128))
    lower=arb(units)/128
    todo=[(0,0)];depth=leaves=0
    while todo:
        n,j=todo.pop()
        t=arb.pi()*arb(2*j+1,1)/(2**n)
        x=real_arb(t,row);y=B*(k*t).sin()
        if all((x-v)**2+y**2>lower**2 for v in (-2,2)):
            depth=max(depth,n);leaves+=1
        else:
            assert n<35,(row["case"],j,units)
            todo.extend([(n+1,2*j),(n+1,2*j+1)])
    # R bounds |z| on the largest residual disk. The embedding radius is R+1/4.
    A=(sum(abs(c) for c in row["coeffs"])/row["divisor"])+B
    R=next(r for r in range(3,25) if r**3-3*r>A+units/128)
    delta=units/128
    # Dyadic choices allow exact checks of the gradient trapping and contraction bounds.
    flow_radius=2.**math.floor(math.log2(delta/(24*R*(R+2))))
    derivative_bound=sum(n*(abs(c)+abs(s)) for n,c,s in terms(row))/row["divisor"]+B*k
    minimum_gradient=9*(delta-flow_radius)/(R+2)
    gain=2**math.ceil(math.log2(2*derivative_bound/(flow_radius*minimum_gradient)))
    radius_q=Fraction(flow_radius)
    g_q=9*(Fraction(units,128)-radius_q)/(R+2)
    assert gain*radius_q*g_q>Fraction(derivative_bound) and g_q>6*R*radius_q
    assert R**3-3*R>Fraction(sum(abs(c) for c in row["coeffs"]),row["divisor"])+B+Fraction(units,128)
    row.update(cut_directions=parts,arb_bits=ctx.prec,clearance_units_128=units,
               accepted_angle_boxes=leaves,maximum_depth=depth,root_modulus_bound=R,
               embedding_radius=R+.25,flow_residual_radius=flow_radius,flow_gain=gain,
               drive_speed_bound=derivative_bound,abs_mu=abs(mu),minimal_drive_period_2pi=
               math.gcd(k,*[n for n,_,_ in terms(row)])==1,return_residual_radius=delta/2)
    return row


def prepare(broad,sparse):
    original=sorted(broad["best_by_degree_and_mu"],key=lambda r:-r["score"])[:6]
    special=sorted(sparse["best"],key=lambda r:r["case"])
    parent=next(r for r in special if r["case"]==47053641)
    examples=[]
    for turns in (1,2,3,5,7):
        h=3*turns
        c=[0]*(2*h+1);c[h+1]=1;c[h+3*turns]=-16
        row=dict(case=48000000+turns,h=h,k=6*turns,coeffs=c,divisor=4,imag_amplitude=3,
                 origin="frequency continuation with a primitive harmonic",parent=parent["case"])
        values=drive(np.arange(2*row["k"])*np.pi/row["k"],row)
        st,mu,_=classify(values.real.copy())
        assert st==0 and abs(mu)==turns
        row["mu"]=int(mu);examples.append(row)
    rows=[dict(r,origin="broad screen") for r in original]+[dict(r,origin="sparse grid") for r in special]+examples
    return [certify(r) for r in rows]


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("broad",type=Path);ap.add_argument("sparse",type=Path);ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args();b=json.loads(a.broad.read_text());s=json.loads(a.sparse.read_text())
    rows=prepare(b,s)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(dict(broad_draws=b["draws"],broad_totals=b["totals"],sparse_draws=s["draws"],
                                    sparse_hits=s["nonzero_mu_screen_hits"],rows=rows),indent=2)+"\n")
    print([(r["case"],r["abs_mu"],r["clearance_units_128"]) for r in rows])
