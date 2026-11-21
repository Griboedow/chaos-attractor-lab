"""Replay directed-interval partitions in Arb for cyclic or symbolic enclosures."""
import argparse
import base64
import hashlib
import json
import zlib
from pathlib import Path
from flint import arb,arb_mat,ctx
from verify_boxes import column,entry,frames,rotate,values
from verify_return import tighten
from verify_cover import monotone_bounds


def parameters(box,bits):
    return [arb(float(v)) if bits is None else arb(float(v),2.**-bits) for v in box["params"]]


def one(box,record,cover,budget,bits):
    origins,matrices,ins=frames(box)
    origin,matrix,inv=origins[0],matrices[0],ins[0]
    p=parameters(box,bits)
    kind=record["kind"];axis=record.get("axis",cover.get("axis",0));steps=record.get("map_steps",cover.get("map_steps",0))
    lo,hi=record["lo"],record["hi"]
    middle=[(arb(v)+arb(w))/2 for v,w in zip(lo,hi)]
    radii=[(arb(w)-arb(v))/2 for v,w in zip(lo,hi)]
    o=origin+matrix*column(middle)
    base=matrix*arb_mat([[radii[i] if i==j else 0 for j in range(3)] for i in range(3)])
    active=[j for j in range(3) if lo[j]!=hi[j]]
    threshold=cover.get("faces",[-1.,1.])
    l,h=arb(threshold[0]),arb(threshold[-1])
    tree=zlib.decompress(base64.b64decode(record["partition_zlib_base64"]))
    pos=extra=leaves=0
    stack=[([0.,0.,0.],[1.,1.,1.],True,0)]
    for_range=list(range(3))
    other=[i for i in for_range if i!=axis]
    def good(b):
        if kind=="face":
            return b[axis]>h if record["side"] else b[axis]<l
        if kind=="slab":
            return b[axis]>h or b[axis]<l or all(abs(b[i])<1 for i in other)
        if kind=="coordinate":
            return abs(b[axis])<1
        if kind=="avoid":
            return any(v>1 or v< -1 for v in b)
        raise ValueError(kind)
    while stack:
        mid,rad,frozen,depth=stack.pop()
        split=3
        if frozen:
            split=tree[pos];pos+=1
        if split==3:
            center=o+base*column(mid)
            wide=o+base*column([arb(v,w) for v,w in zip(mid,rad)])
            jac=base
            for _ in range(steps):
                center,_=rotate(values(center),p)
                wide,jac=rotate(values(wide),p,jac)
                wide=tighten(wide,center,jac,rad)
            at=inv*(center-origin);jt=inv*jac
            enclosures=values(tighten(inv*(wide-origin),at,jt,rad))
            outputs=[axis] if kind in ("face","coordinate") else for_range
            if not good(enclosures):
                for i in outputs:
                    enclosures[i]=monotone_bounds(enclosures[i],jt,mid,rad,o,base,origin,inv,p,steps,i)
            if good(enclosures):
                leaves+=1;continue
            extra+=1
            if extra>budget or depth>74:
                raise RuntimeError(f"Arb extra-subdivision limit for {kind} at step {steps}: {extra}")
            largest=max(rad[j] for j in active)
            eligible=[j for j in active if rad[j]>=largest/64]
            split=max(eligible,key=lambda j:max(float(jt[k,j].abs_upper())*rad[j] for k in outputs))
            frozen=False
        elif split not in for_range:
            raise ValueError("Invalid stored split")
        dr=rad.copy();dr[split]/=2
        for sign in (-1,1):
            m=mid.copy();m[split]+=sign*dr[split]
            if not (arb(m[split])-sign*arb(dr[split])==arb(mid[split])):
                raise RuntimeError("Nonexact dyadic split")
            stack.append((m,dr,frozen,depth+1))
    if pos!=len(tree):
        raise RuntimeError("Trailing partition data")
    return dict(kind=kind,map_steps=steps,axis=axis,original_nodes=pos,extra_splits=extra,accepted_leaves=leaves)


def replay(box,cert,budget):
    ctx.prec=128
    if cert["case"]!=box["case"]:
        raise ValueError("Case mismatch")
    if cert.get("box_digest") and cert["box_digest"]!=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest():
        raise ValueError("Box digest mismatch")
    records=cert.get("checks",cert.get("regions"))
    bits=cert.get("parameter_radius_bits")
    if "regions" in cert:
        cover=cert["cover"]
        faces=cover["faces"]
        axis=cover["axis"]
        if len(faces)!=4 or not (-1<=faces[0]<faces[1]<faces[2]<faces[3]<=1):
            raise ValueError("Inconsistent two-slab coordinates")
        if cover["map_steps"]!=box["case_order"] or cover["case"]!=box["case"]:
            raise ValueError("Symbolic return-map mismatch")
        sides=cover["side"]
        if any(sorted(sides[i:i+2])!=[0,1] for i in (0,2)):
            raise ValueError("Each slab must cross both outer faces")
        spec=[]
        for i in range(4):
            lo=[-1.]*3;hi=[1.]*3;lo[axis]=hi[axis]=faces[i]
            spec.append(("face",lo,hi,sides[i]))
        for i in (0,2):
            lo=[-1.]*3;hi=[1.]*3;lo[axis],hi[axis]=faces[i:i+2]
            spec.append(("slab",lo,hi,None))
        if [(r["kind"],r["lo"],r["hi"],r["side"]) for r in records]!=spec:
            raise ValueError("Incomplete symbolic certificate")
    out=[]
    for i,r in enumerate(records):
        out.append(one(box,r,cert.get("cover",{}),budget,bits))
        print(box["case"],i,r["kind"],out[-1]["accepted_leaves"],flush=True)
    report=dict(case=box["case"],precision_bits=128,independent_regions=out,parameter_radius_bits=bits)
    if "checks" in cert:
        known=sorted((r["map_steps"],r["kind"],r.get("axis",0)) for r in records)
        order=box["case_order"]
        expected=sorted([(j,"avoid",0) for j in range(1,order)]+[(order,"coordinate",i) for i in range(3)])
        if known!=expected:
            raise ValueError("Incomplete return certificate")
        got=None
        for prec in (256,1024,4096):
            ctx.prec=prec;o,_,inv=frames(box)
            got=entry(box,o,inv,parameters(box,bits),4096)
            if got:
                got["precision_bits"]=prec;break
        report.update(phase_regions=order,primary_start_entry=got)
    else:
        report.update(cover=cert["cover"],all_faces_and_full_slabs_verified=True)
    return report


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("certificate",type=Path)
    ap.add_argument("--budget",type=int,default=500_000)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=replay(json.loads(a.box.read_text()),json.loads(a.certificate.read_text()),a.budget)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+"\n")
