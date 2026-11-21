"""Certify slab crossings and transverse containment for a two-symbol horseshoe."""
import argparse
import base64
import json
import zlib
from pathlib import Path
from flint import arb,arb_mat,ctx
from verify_boxes import column,frames,rotate,values
from verify_return import tighten


def monotone_bounds(bounds,jt,mid,rad,source_origin,source_linear,target_origin,target_inverse,p,steps,coordinate):
    signs=[]
    remainder=arb(0)
    for j in range(3):
        d=jt[coordinate,j]
        sign=1 if d>0 else (-1 if d<0 else 0)
        signs.append(sign)
        if sign==0:
            remainder+=d.abs_upper()*rad[j]
    if not any(signs):
        return bounds
    extremes=[]
    for direction in (-1,1):
        point=source_origin+source_linear*column([arb(v)+direction*s*arb(h) for v,s,h in zip(mid,signs,rad)])
        for _ in range(steps):
            point,_=rotate(values(point),p)
        extremes.append((target_inverse*(point-target_origin))[coordinate,0])
    enclosure=(extremes[0]-remainder).union(extremes[1]+remainder)
    return bounds.intersection(enclosure)


def check_region(packet,cover,lo,hi,kind,side,budget):
    origins,matrices,inverses=frames(packet)
    origin,matrix,inv=origins[0],matrices[0],inverses[0]
    p=[arb(v) for v in packet["params"]]
    axis=cover["axis"];lower,upper=arb(cover["faces"][0]),arb(cover["faces"][-1])
    mid0=column([(arb(l)+arb(h))/2 for l,h in zip(lo,hi)])
    radius=[(arb(h)-arb(l))/2 for l,h in zip(lo,hi)]
    o=origin+matrix*mid0
    dmat=matrix*arb_mat([[radius[i] if i==j else 0 for j in range(3)] for i in range(3)])
    stack=[([0.,0.,0.],[1.,1.,1.],0)]
    leaves=nodes=maxdepth=0;part=bytearray()
    while stack:
        m,d,depth=stack.pop();nodes+=1;maxdepth=max(depth,maxdepth)
        if nodes>budget or depth>80:
            raise RuntimeError(f"Subdivision limit for {kind} at {nodes} nodes, depth {depth}")
        point=o+dmat*column(m)
        wide=o+dmat*column([arb(a,b) for a,b in zip(m,d)])
        jac=dmat
        for _ in range(cover["map_steps"]):
            point,_=rotate(values(point),p)
            wide,jac=rotate(values(wide),p,jac)
            wide=tighten(wide,point,jac,d)
        at=inv*(point-origin);jt=inv*jac
        box=values(tighten(inv*(wide-origin),at,jt,d));w=values(at)
        active=[axis] if kind=="face" else [axis]+[i for i in range(3) if i!=axis]
        if kind=="face":
            failed=w[axis]<upper if side else w[axis]>lower
        else:
            failed=lower<w[axis] and w[axis]<upper and any(abs(w[i])>1 for i in active[1:])
        if failed:
            raise RuntimeError(f"Center witnesses a failed {kind} condition")
        for attempt in range(2):
            if kind=="face":
                passed=box[axis]>upper if side else box[axis]<lower
            else:
                passed=box[axis]>upper or box[axis]<lower or all(abs(box[i])<1 for i in active[1:])
            if passed or attempt==1:
                break
            for q in active:
                if q==axis or not (abs(box[q])<1):
                    box[q]=monotone_bounds(box[q],jt,m,d,o,dmat,origin,inv,p,cover["map_steps"],q)
        if passed:
            part.append(3);leaves+=1;continue
        allowed=[j for j in range(3) if lo[j]!=hi[j]]
        largest=max(d[j] for j in allowed)
        allowed=[j for j in allowed if d[j]>=largest/64]
        split=max(allowed,key=lambda j:max(float(jt[k,j].abs_upper())*d[j] for k in active))
        part.append(split);dr=d.copy();dr[split]/=2
        for sign in (-1,1):
            c=m.copy();c[split]+=sign*dr[split]
            if not (arb(c[split])-sign*arb(dr[split])==arb(m[split])):
                raise RuntimeError("Nonexact dyadic split")
            stack.append((c,dr,depth+1))
    return dict(kind=kind,side=side,lo=lo,hi=hi,leaves=leaves,nodes=nodes,depth=maxdepth,
                partition_zlib_base64=base64.b64encode(zlib.compress(part,9)).decode("ascii"))


def certify(packet,cover,out,budget):
    ctx.prec=128
    a=cover["axis"];f=cover["faces"]
    if not (-1<=f[0]<f[1]<f[2]<f[3]<=1 and cover["map_steps"]>=1):
        raise ValueError("Invalid slab coordinates")
    if not all(cover["side"][i]!=cover["side"][i+1] for i in (0,2)):
        raise ValueError("Each slab needs two opposite face images")
    tests=[]
    for j in range(4):
        lo=[-1.,-1.,-1.];hi=[1.,1.,1.]
        lo[a]=hi[a]=f[j]
        tests.append((lo,hi,"face",cover["side"][j]))
    for j in (0,2):
        lo=[-1.,-1.,-1.];hi=[1.,1.,1.]
        lo[a],hi[a]=f[j:j+2]
        tests.append((lo,hi,"slab",None))
    checks=[]
    for i,(lo,hi,kind,side) in enumerate(tests):
        partial=out.with_name(out.stem+f"-part-{i}.json")
        if partial.exists():
            record=json.loads(partial.read_text())
            if record.get("cover")!=cover or record.get("box")!=packet:
                raise RuntimeError("Incompatible partial certificate")
            check=record["check"]
        else:
            check=check_region(packet,cover,lo,hi,kind,side,budget)
            partial.write_text(json.dumps(dict(cover=cover,box=packet,check=check))+"\n")
        checks.append(check)
        print(packet["case"],kind,i,check["nodes"],flush=True)
    return dict(case=packet["case"],precision_bits=128,cover=cover,regions=checks,
                entropy_lower_bound=f"log(2)/{cover['map_steps']} per original iterate")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("proposals",type=Path)
    ap.add_argument("--index",type=int,default=0);ap.add_argument("--budget",type=int,default=350_000)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    packet=json.loads(a.box.read_text());cover=json.loads(a.proposals.read_text())[a.index]
    result=certify(packet,cover,a.output,a.budget)
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result["entropy_lower_bound"],flush=True)
