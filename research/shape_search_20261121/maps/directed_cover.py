"""Independent directed-interval search or replay of two symbolic branches."""
import argparse
import base64
import json
import zlib
from fractions import Fraction
from pathlib import Path
from directed_return import (enclosure,frame,intersect_mean,maximum_abs,minus,plus,rotate,times,
                             transform,trig_constants)
import numpy as np
from numba import njit


@njit(cache=True)
def evaluate(mid,rad,steps,p,base,origin,target,inv,pi2,cs,eps,axis,kind,side,lower,upper):
    centered=np.zeros((3,1,2));wide=np.zeros((3,1,2))
    for i in range(3):
        centered[i,0]=(mid[i],mid[i]);wide[i,0]=plus(centered[i,0],(-rad[i],rad[i]))
    centered=transform(base,centered);wide=transform(base,wide)
    jet=np.empty((3,4,2))
    for i in range(3):
        centered[i,0]=plus(centered[i,0],origin[i])
        jet[i,0]=plus(wide[i,0],origin[i])
        for j in range(3):
            jet[i,j+1]=base[i,j]
    for _ in range(steps):
        centered=rotate(centered,p,pi2,cs,eps,True)
        jet=rotate(jet,p,pi2,cs,eps,True)
        intersect_mean(jet,centered,rad)
    for i in range(3):
        jet[i,0]=minus(jet[i,0],target[i]);centered[i,0]=minus(centered[i,0],target[i])
    jet=transform(inv,jet);centered=transform(inv,centered)
    intersect_mean(jet,centered,rad)
    q=jet[axis,0];c=centered[axis,0]
    if kind==0:
        passed=q[0]>upper if side else q[1]<lower
        fail=c[1]<upper if side else c[0]>lower
    elif kind==1:
        passed=q[0]>upper or q[1]<lower
        tight=True;outside=False
        for i in range(3):
            if i!=axis:
                tight=tight and jet[i,0,0]>-1 and jet[i,0,1]<1
                outside=outside or centered[i,0,0]>1 or centered[i,0,1]<-1
        passed=passed or tight
        fail=c[0]>lower and c[1]<upper and outside
    elif kind==2:
        passed=False;fail=True
        for i in range(3):
            passed=passed or jet[i,0,0]>1 or jet[i,0,1]<-1
            fail=fail and centered[i,0,0]>-1 and centered[i,0,1]<1
    else:
        passed=q[0]>-1 and q[1]<1
        fail=c[0]>1 or c[1]<-1
    weight=np.zeros(3)
    for j in range(3):
        for i in range(3):
            if (kind==0 or kind==3) and i!=axis:
                continue
            weight[j]=max(weight[j],rad[j]*maximum_abs(jet[i,j+1]))
    return (1 if passed else (-1 if fail else 0)),weight


def region_frame(box,lo,hi):
    M,o,inv=frame(box)
    q=Fraction.from_float
    middle=[enclosure((q(a)+q(b))/2) for a,b in zip(lo,hi)]
    extent=[enclosure((q(b)-q(a))/2) for a,b in zip(lo,hi)]
    base=np.empty_like(M);start=o.copy()
    for i in range(3):
        for j in range(3):
            start[i]=plus(start[i],times(M[i,j],middle[j]))
            base[i,j]=times(M[i,j],extent[j])
    return base,start,o,inv


def region(box,cover,lo,hi,kind,side,frozen,budget,parameter_bits=None):
    base,o,target,inv=region_frame(box,lo,hi)
    p=np.array([[v,v] for v in box["params"]])
    if parameter_bits is not None:
        radius=Fraction(1,2**parameter_bits)
        p=np.array([[enclosure(Fraction.from_float(v)-radius)[0],
                     enclosure(Fraction.from_float(v)+radius)[1]] for v in box["params"]])
    pi2,cs,eps=trig_constants()
    tree=zlib.decompress(base64.b64decode(frozen)) if frozen else b""
    pos=accepted=extra=nodes=0
    whole=(0,0,0);stack=[(whole,whole,bool(frozen))];encoding=bytearray()
    axis=cover["axis"];faces=cover["faces"]
    active=[j for j in range(3) if lo[j]!=hi[j]]
    while stack:
        idx,depths,stored=stack.pop();nodes+=1
        if nodes>budget or max(depths)>48:
            raise RuntimeError(f"Budget for {kind}: {nodes} nodes, depths {depths}")
        split=3
        if stored:
            split=tree[pos];pos+=1
            if split>3:
                raise RuntimeError("Malformed tree")
        if split==3:
            rad=np.array([2.**-d for d in depths]);mid=np.array([(2*i+1)*h-1 for i,h in zip(idx,rad)])
            modes={"face":0,"slab":1,"avoid":2,"coordinate":3}
            state,weights=evaluate(mid,rad,cover["map_steps"],p,base,o,target,inv,pi2,cs,eps,axis,
                                   modes[kind],side or 0,faces[0],faces[-1])
            if state==1:
                encoding.append(3);accepted+=1;continue
            if state==-1:
                raise RuntimeError(f"Point witnesses a failed {kind} condition")
            largest=max(rad[j] for j in active)
            eligible=[j for j in active if rad[j]>=largest/64]
            split=max(eligible,key=lambda j:weights[j]);extra+=1;stored=False
        encoding.append(split)
        for b in (0,1):
            ii=list(idx);dd=list(depths)
            ii[split]=2*ii[split]+b;dd[split]+=1
            stack.append((tuple(ii),tuple(dd),stored))
    if pos!=len(tree):
        raise RuntimeError("Trailing tree nodes")
    return dict(kind=kind,side=side,lo=lo,hi=hi,nodes=nodes,leaves=accepted,new_splits=extra,
                partition_zlib_base64=base64.b64encode(zlib.compress(encoding,9)).decode("ascii"))


def run(box,cover,out,budget,primary=None,parameter_bits=None):
    axis=cover["axis"];faces=cover["faces"];spec=[]
    for j in range(4):
        lo=[-1.]*3;hi=[1.]*3;lo[axis]=hi[axis]=faces[j]
        spec.append((lo,hi,"face",cover["side"][j]))
    for j in (0,2):
        lo=[-1.]*3;hi=[1.]*3;lo[axis],hi[axis]=faces[j:j+2]
        spec.append((lo,hi,"slab",None))
    if primary and primary["cover"]!=cover:
        raise RuntimeError("Proposed faces do not match source certificate")
    all_regions=[]
    for k,(lo,hi,kind,side) in enumerate(spec):
        tree=primary["regions"][k]["partition_zlib_base64"] if primary else None
        hit=region(box,cover,lo,hi,kind,side,tree,budget,parameter_bits)
        all_regions.append(hit)
        print(box["case"],k,kind,hit["nodes"],flush=True)
    out.write_text(json.dumps(dict(case=box["case"],cover=cover,regions=all_regions,parameter_radius_bits=parameter_bits,
                                    method="outward binary64 intervals, rational trigonometric bounds"),indent=2)+"\n")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("proposals",type=Path)
    ap.add_argument("--index",type=int,default=0);ap.add_argument("--primary",type=Path)
    ap.add_argument("--budget",type=int,default=3_000_000);ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--parameter-radius-bits",type=int)
    a=ap.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    run(json.loads(a.box.read_text()),json.loads(a.proposals.read_text())[a.index],a.output,a.budget,
        json.loads(a.primary.read_text()) if a.primary else None,a.parameter_radius_bits)
