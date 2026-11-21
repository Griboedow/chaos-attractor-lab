"""Certify a return box and exclusion of every earlier iterate, using Arb."""
import argparse
import base64
import json
import zlib
from pathlib import Path
from flint import arb,ctx
from verify_boxes import column,entry,frames,rotate,values


def tighten(natural,center,jac,radii):
    result=[]
    for i in range(3):
        rem=sum(jac[i,j].abs_upper()*radii[j] for j in range(3))
        ball=center[i,0]+arb(0,rem.abs_upper())
        result.append(natural[i,0].intersection(ball))
    return column(result)


def check(packet,budget,prec,depth_limit,parameter_bits=None):
    ctx.prec=prec
    p=[arb(v) if parameter_bits is None else arb(v,2.**-parameter_bits) for v in packet["params"]]
    if not (p[0]>0 and p[0]<1):
        raise ValueError("The full contraction interval must lie in (0,1)")
    origins,bases,inverses=frames(packet)
    base,frame,inv=origins[0],bases[0],inverses[0]
    order=packet["case_order"]
    stack=[([0.,0.,0.],[1.,1.,1.],0)]
    partition=bytearray()
    count=leaves=depth=0;max_return_bound=0.;min_avoidance=float("inf")
    while stack:
        mid,rad,d=stack.pop()
        count+=1;depth=max(d,depth)
        if count>budget or d>depth_limit:
            raise RuntimeError(f"Subdivision limit: {count} nodes, depth {d}")
        xmid=base+frame*column(mid)
        xbox=base+frame*column([arb(m,r) for m,r in zip(mid,rad)])
        jac=frame
        scores=[0.,0.,0.]
        accepted=True;leaf_gap=float("inf")
        for k in range(1,order+1):
            xmid,_=rotate(values(xmid),p)
            xbox,jac=rotate(values(xbox),p,jac)
            xbox=tighten(xbox,xmid,jac,rad)
            cmid=inv*(xmid-base)
            jlocal=inv*jac
            u=tighten(inv*(xbox-base),cmid,jlocal,rad)
            bounds=values(u)
            if k<order:
                safe=any(v>1 or v< -1 for v in bounds)
                if all(abs(q)<1 for q in values(cmid)):
                    raise RuntimeError(f"Earlier image of a subbox center enters the reference box at iterate {k}")
                if safe:
                    leaf_gap=min(leaf_gap,max(float(v.abs_lower()) for v in bounds)-1)
            else:
                safe=all(abs(v)<1 for v in bounds)
                if any(abs(v)>1 for v in values(cmid)):
                    raise RuntimeError("Subbox center lies outside its return box")
                bound=max(float(v.abs_upper()) for v in bounds)
            if not safe:
                accepted=False
                for j in range(3):
                    scores[j]=max(scores[j],max(float(jlocal[i,j].abs_upper())*rad[j] for i in range(3)))
        if accepted:
            partition.append(3)
            leaves+=1;max_return_bound=max(max_return_bound,bound)
            min_avoidance=min(min_avoidance,leaf_gap)
            continue
        # Extremely wide derivative enclosures can overflow the heuristic at
        # coarse scales. Keep every axis eligible for eventual refinement.
        eligible=[j for j in range(3) if rad[j]>=max(rad)/64]
        axis=max(eligible,key=lambda j:scores[j])
        partition.append(axis)
        child=rad.copy();child[axis]/=2
        for direction in (-1,1):
            c=mid.copy();c[axis]+=direction*child[axis]
            if not (arb(c[axis])-direction*arb(child[axis])==arb(mid[axis])):
                raise RuntimeError("Inexact dyadic partition")
            stack.append((c,child,d+1))
    entered=None
    for bits in (256,1024,4096):
        ctx.prec=bits;o,_,invlist=frames(packet)
        entered=entry(packet,o,invlist,p,4096)
        if entered:
            entered["precision_bits"]=bits
            break
    ctx.prec=prec
    return dict(case=packet["case"],phase_regions=order,precision_bits=prec,method="strict return and earlier-image exclusion",
                independent_parameter_radius_bits=parameter_bits,
                verified_leaves=leaves,visited_nodes=count,max_depth=depth,
                partition_zlib_base64=base64.b64encode(zlib.compress(partition,9)).decode("ascii"),
                largest_normalized_return_upper_approximation=max_return_bound,
                smallest_normalized_earlier_exclusion_gap_approximation=min_avoidance,primary_start_entry=entered)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("packet",type=Path)
    ap.add_argument("--precision",type=int,default=128)
    ap.add_argument("--budget",type=int,default=300_000)
    ap.add_argument("--depth",type=int,default=70)
    ap.add_argument("--parameter-radius-bits",type=int,help="Enclose every coefficient within +/-2^(-bits)")
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    res=check(json.loads(a.packet.read_text()),a.budget,a.precision,a.depth,a.parameter_radius_bits)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(res,indent=2)+"\n")
    print(res,flush=True)
