"""Certify long attracting cycles using short derivative blocks and balanced products."""
import argparse
import csv
import json
import math
from pathlib import Path

from flint import arb,arb_mat,ctx
import mpmath as mp
import numpy as np

from period_windows import advance,bases,parameters
from search import save
from sink import step


def choose(rows,n):
    candidates=[]
    for row in rows:
        if row["direction"]!="forward" or int(row["period_multiple_of_39"])!=n:
            continue
        candidates.append(row)
    if len(candidates)>35:
        candidates=[candidates[i] for i in np.linspace(0,len(candidates)-1,35).astype(int)]
    ranked=[]
    for row in candidates:
        ix,den=int(row["numerator"]),int(row["denominator"])
        p=parameters(ix/den)
        x=[float(row[k]) for k in "xyz"]
        z=x
        ds=np.eye(3)
        for _ in range(n):
            z,j=step(z,p,math.exp)
            ds=np.asarray(j)@ds
        rate,v=min(bases(ds),key=lambda r:r[0])
        ranked.append((rate,ix,den,x,v))
    best=min(ranked,key=lambda r:r[0])
    assert best[0]<.999
    return best


def propose(rows,n):
    rate,ix,den,x,v=choose(rows,n)
    mp.mp.dps=94
    p=parameters(mp.mpf(ix)/den)
    x=mp.matrix([mp.mpf(q) for q in x])
    for _ in range(6):
        end,j=advance(x,p,n)
        delta=mp.lu_solve(j-mp.eye(3),end-x)
        x-=delta
        if max(abs(q) for q in delta)<mp.mpf("1e-78"):
            break
    end,_=advance(x,p,n)
    assert max(abs(q) for q in end-x)<mp.mpf("1e-72")
    return {"period":n,"t_numerator":ix,"t_denominator":den,
            "center":[mp.nstr(q,85) for q in x],"basis":[[format(q,".17g") for q in t] for t in v],
            "radius":"1e-35","block_start_bound":"1e-21","block_length":8}


def identity(n):
    return arb_mat([[int(i==j) for j in range(n)] for i in range(n)])


def prefix_product(stack,m):
    # Each entry contains 2^level consecutive blocks. Balanced multiplication
    # limits dependency inflation when the short Jacobians are nonnormal.
    level=0
    while stack and stack[-1][0]==level:
        m=m*stack.pop()[1]
        level+=1
    stack.append((level,m))
    product=identity(4)
    for _,m in stack:
        product=m*product
    return product


def verify(row):
    mp.mp.dps=94
    ctx.prec=350
    n,k=row["period"],row["block_length"]
    ix,den=row["t_numerator"],row["t_denominator"]
    p=parameters(arb(ix)/den)
    proposal_params=parameters(mp.mpf(ix)/den)
    v=arb_mat([[arb(t) for t in r] for r in row["basis"]])
    eps=arb(row["radius"]); u=arb(0,row["radius"])
    tube=arb(row["block_start_bound"])
    err=arb(0,row["block_start_bound"])
    worst=arb(0)
    initial=v*arb_mat([[u],[u],[u]])
    assert all(abs(initial[j,0])<tube for j in range(3))
    x=list(map(mp.mpf,row["center"]))
    anchor=list(map(arb,row["center"]))
    assert all(a+err>arb("1e-22") and a+err<10 for a in anchor)
    first=list(anchor)
    stack=[]
    divisors=[]
    above=False
    for index in range(0,n,k):
        length=min(k,n-index)
        # The next anchor is untrusted and is checked by its block residual.
        for _ in range(length):
            x,_=step(x,proposal_params,mp.exp)
        nxt=first if index+length==n else [arb(mp.nstr(t,84)) for t in x]
        z=list(anchor)
        whole=[q+err for q in anchor]
        d=identity(3)
        for s in range(1,length+1):
            z,_=step(z,p,lambda t:t.exp())
            whole,j=step(whole,p,lambda t:t.exp())
            d=arb_mat(j)*d
            phase=index+s
            if phase<n and n%phase==0:
                assert any(not val.overlaps(first[h]+err) for h,val in enumerate(whole))
                divisors.append(phase)
            above=above or any(q>arb("2.85") for q in whole)
        b=[z[i]-nxt[i] for i in range(3)]
        mat=arb_mat([[d[i,j] for j in range(3)]+[b[i]] for i in range(3)]+[[0,0,0,1]])
        total=prefix_product(stack,mat)
        a=arb_mat([[total[i,j] for j in range(3)] for i in range(3)])
        offset=arb_mat([[total[i,3]] for i in range(3)])
        end=a*initial+offset
        bound=max(abs(end[i,0]) for i in range(3))
        assert bound<tube
        if bound.upper()>worst:
            worst=bound.upper()
        anchor=nxt
    # The linear part encloses the derivative on the whole initial ball.
    # The affine column encloses the transported endpoint error for its center.
    inv=v.inv()
    df=inv*a*v
    sq=df.transpose()*df
    upper=max(sq[i,i]+sum(abs(sq[i,j]) for j in range(3) if i!=j) for i in range(3))
    normalized=inv*offset
    deviation=sum(normalized[i,0]**2 for i in range(3)).sqrt()/eps
    assert upper.sqrt()+deviation<1,(str(upper),str(deviation))
    return {"precision_bits":ctx.prec,"blocks":math.ceil(n/k),"strict_contraction_bound":str(upper.sqrt()),
            "normalized_return_residual_bound":str(deviation),"separated_divisors":divisors,
            "one_phase_above_symbolic_box":bool(above),"highest_validated_block_radius":str(worst)}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zoom",type=Path)
    ap.add_argument("--periods",nargs="+",type=int,default=[2496,4992])
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--verify",action="store_true")
    args=ap.parse_args()
    if args.verify:
        for r in json.loads(args.out.read_text()):
            verify(r)
    else:
        with args.zoom.open(newline="",encoding="utf-8") as f:
            rows=list(csv.DictReader(f))
        certs=[]
        for n in args.periods:
            row=propose(rows,n)
            row["verified"]=verify(row)
            certs.append(row)
            save(args.out,certs)
            print("Certified block return of least period",n,flush=True)
