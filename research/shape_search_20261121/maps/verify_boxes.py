"""Verify strict cyclic trapping using outward-rounded Arb real-ball arithmetic.

All JSON coefficients are the exact binary64 values produced by parsing them.
For each box, transform the recurrence into its successor's normalized frame.
The mean-value enclosure uses the interval Jacobian over the entire subbox.
Partitions bisect the exact dyadic input cube and exhaust it before success.
"""
import argparse
import json
import math
from pathlib import Path
from flint import arb,arb_mat,ctx


def matrix(a):
    return arb_mat([[arb(x) for x in row] for row in a])


def column(x):
    return arb_mat([[arb(v)] for v in x])


def frames(packet):
    origins=[];bases=[]
    for c,b,lo,hi in zip(packet["center"],packet["basis"],packet["lo"],packet["hi"]):
        basis=matrix(b)
        mid=[(arb(h)+arb(l))/2 for l,h in zip(lo,hi)]
        rad=[(arb(h)-arb(l))/2 for l,h in zip(lo,hi)]
        if not all(r>0 for r in rad):
            raise ValueError("Nonpositive box extent")
        origins.append(column(c)+basis*column(mid))
        bases.append(basis*arb_mat([[rad[j] if i==j else 0 for j in range(3)] for i in range(3)]))
    return origins,bases,[m.inv() for m in bases]


def rotate(x,p,jac=None):
    x=[p[0]*v for v in x]
    x[0]+=1-p[0]
    if jac is not None:
        jac=jac*p[0]
    for i in range(3):
        j,k=(i+1)%3,(i+2)%3
        w=x[k];a=p[1+i]+p[4+i]*w+p[7+i]*w**3
        s,c=a.sin_cos()
        u,y=c*x[i]-s*x[j],s*x[i]+c*x[j]
        if jac is not None:
            slope=p[4+i]+3*p[7+i]*w*w
            rot=arb_mat([[int(t==q) for q in range(3)] for t in range(3)])
            rot[i,i],rot[i,j],rot[i,k]=c,-s,-y*slope
            rot[j,i],rot[j,j],rot[j,k]=s,c,u*slope
            jac=rot*jac
        x[i],x[j]=u,y
    return column(x),jac


def values(x):
    return [x[i,0] for i in range(3)]


def disjoint(origins,bases,inverses):
    # A single separating linear functional suffices for each pair. Each
    # candidate normal is made exact by taking its floating-point midpoint.
    candidates=[[float(i==j) for j in range(3)] for i in range(3)]
    clearances=[]
    for i in range(len(origins)):
        for j in range(i):
            norms=candidates+[[float(o[a,b].mid()) for b in range(3)]
                              for o in (inverses[i],inverses[j]) for a in range(3)]
            norms.append([float((origins[i][k,0]-origins[j][k,0]).mid()) for k in range(3)])
            best=0.
            for n in norms:
                row=matrix([n])
                distance=abs((row*(origins[i]-origins[j]))[0,0])
                radii=sum(abs((row*bases[h])[0,k]) for h in (i,j) for k in range(3))
                if distance>radii:
                    gap=(distance-radii)/sum(arb(v)**2 for v in n).sqrt()
                    best=max(best,float(gap.lower()))
            if best<=0:
                raise RuntimeError(f"Missing separating plane for phases {i}, {j}")
            clearances.append(best)
    return min(clearances)


def inclusion(origin,basis,next_origin,next_inverse,p,budget):
    stack=[([0.,0.,0.],[1.,1.,1.],0)]
    leaves=nodes=depth=0;worst=0.
    while stack:
        mid,rad,d=stack.pop()
        depth=max(depth,d);nodes+=1
        if nodes>budget or d>80:
            raise RuntimeError(f"Subdivision limit at {nodes} visited nodes, depth {d}")
        at=origin+basis*column(mid)
        f,_=rotate(values(at),p)
        centered=next_inverse*(f-next_origin)
        if any(abs(q)>1 for q in values(centered)):
            raise RuntimeError("A sampled center lies outside the target box")
        xx=origin+basis*column([arb(m,r) for m,r in zip(mid,rad)])
        natural,jac=rotate(values(xx),p,basis)
        direct=values(next_inverse*(natural-next_origin))
        jt=next_inverse*jac
        bounds=[]
        for k in range(3):
            remainder=sum(abs(jt[k,j])*rad[j] for j in range(3))
            derivative_bound=(abs(centered[k,0])+remainder).abs_upper()
            natural_bound=direct[k].abs_upper()
            # Either enclosure individually contains this output coordinate.
            bounds.append(natural_bound if natural_bound<derivative_bound else derivative_bound)
        if all(q<1 for q in bounds):
            leaves+=1;worst=max(worst,max(float(q) for q in bounds))
            continue
        # The split heuristic affects speed only. Both children are checked.
        weights=[max(float(jt[k,j].abs_upper())*rad[j] for k in range(3)) for j in range(3)]
        axis=max(range(3),key=lambda j:weights[j])
        child=rad.copy();child[axis]/=2
        if mid[axis]+child[axis]==mid[axis] or mid[axis]-child[axis]==mid[axis]:
            raise RuntimeError("Dyadic center exhausted binary64 precision")
        for sign in (-1,1):
            c=mid.copy();c[axis]+=sign*child[axis]
            # A binary64 center and radius must preserve the endpoints of the
            # dyadic parent exactly. Check in Arb, including the shared face.
            if (arb(c[axis])-sign*arb(child[axis]))!=arb(mid[axis]):
                raise RuntimeError("Inexact partition")
            stack.append((c,child,d+1))
    return dict(leaves=leaves,nodes=nodes,max_depth=depth,largest_accepted_bound_approximation=worst)


def entry(packet,origins,inverses,p,steps):
    x=[arb(.1),arb(.2),arb(.3)]
    for n in range(steps+1):
        q=column(x)
        for i,(o,inv) in enumerate(zip(origins,inverses)):
            u=values(inv*(q-o))
            if all(abs(v)<1 for v in u):
                return dict(iterate=n,phase=i,maximum_normalized_coordinate_upper=max(float(v.abs_upper()) for v in u))
        q,_=rotate(x,p);x=values(q)
        if max(float(v.rad()) for v in x)>.5:
            break
    return None


def verify(packet,prec,budget):
    ctx.prec=prec
    p=[arb(float(v)) for v in packet["params"]]
    origins,bases,inverses=frames(packet)
    gap=disjoint(origins,bases,inverses)
    checked=[]
    for i in range(len(bases)):
        j=(i+1)%len(bases)
        checked.append(inclusion(origins[i],bases[i],origins[j],inverses[j],p,budget))
    # Re-create exact-real enclosures at higher precision for finite entry.
    entered=None
    for bits in (256,1024,4096):
        ctx.prec=bits;origins,_,inverses=frames(packet)
        entered=entry(packet,origins,inverses,p,4096)
        if entered:
            entered["precision_bits"]=bits
            break
    ctx.prec=prec
    return dict(case=packet["case"],phase_regions=len(bases),precision_bits=prec,
                minimum_separating_gap_lower_approximation=gap,strict_inclusions=checked,
                primary_start_entry=entered)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("packet",type=Path)
    ap.add_argument("--precision",type=int,default=128)
    ap.add_argument("--budget",type=int,default=500_000)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    data=verify(json.loads(a.packet.read_text()),a.precision,a.budget)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(data,indent=2)+"\n")
    print(a.packet.stem,data["phase_regions"],sum(d["nodes"] for d in data["strict_inclusions"]),data["primary_start_entry"],flush=True)
