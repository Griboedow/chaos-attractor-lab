"""Replay a Hausdorff bound with directed intervals and direct derivative minors."""
import argparse
import base64
import hashlib
import json
import zlib
from fractions import Fraction
from pathlib import Path
from directed_return import enclosure,frame,intersect_mean,maximum_abs,minus,overlap,plus,rotate,times,transform,trig_constants
import numpy as np
from numba import njit


@njit(cache=True)
def minors(derivative,offset):
    result=np.empty((3,3,2))
    for i in range(3):
        a,b=(i+1)%3,(i+2)%3
        for j in range(3):
            x,y=(j+1)%3+offset,(j+2)%3+offset
            result[i,j]=minus(times(derivative[a,x],derivative[b,y]),times(derivative[a,y],derivative[b,x]))
    return result


@njit(cache=True)
def sufficient(mid,rad,p,M,o,pi2,cs,eps,steps,num,den):
    c=np.zeros((3,1,2));w=c.copy()
    for i in range(3):
        c[i,0]=(mid[i],mid[i])
        w[i,0]=plus(c[i,0],(-rad[i],rad[i]))
    c=transform(M,c);w=transform(M,w)
    jet=np.zeros((3,7,2))
    for i in range(3):
        c[i,0]=plus(c[i,0],o[i])
        jet[i,0]=plus(w[i,0],o[i])
        for j in range(3):
            jet[i,j+1]=M[i,j]
        jet[i,4+i]=(1.,1.)
    area_matrix=np.zeros((3,3,2))
    for i in range(3):
        area_matrix[i,i]=(1.,1.)
    for _ in range(steps):
        local=np.zeros((3,4,2))
        for i in range(3):
            local[i,0]=jet[i,0]
            local[i,1+i]=(1.,1.)
        local=rotate(local,p,pi2,cs,eps,True)
        area_matrix=transform(minors(local,1),area_matrix)
        c=rotate(c,p,pi2,cs,eps,True)
        jet=rotate(jet,p,pi2,cs,eps,True)
        intersect_mean(jet,c,rad)
    tangent=0.;area=0.
    direct=minors(jet,4)
    for i in range(3):
        for j in range(3):
            entry=maximum_abs(jet[i,j+4])
            minor=overlap(direct[i,j],area_matrix[i,j])
            wedge=maximum_abs(minor)
            tangent=plus((tangent,tangent),times((entry,entry),(entry,entry)))[1]
            area=plus((area,area),times((wedge,wedge),(wedge,wedge)))[1]
    bound=(1.,1.)
    for _ in range(den-num):
        bound=times(bound,(tangent,tangent))
    for _ in range(num):
        bound=times(bound,(area,area))
    return bound[1]


def check(box,source,budget):
    digest=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest()
    if source["box_digest"]!=digest or source["map_steps"]!=box["case_order"]:
        raise ValueError("Certificate and spatial frame differ")
    tree=zlib.decompress(base64.b64decode(source["partition_zlib_base64"]))
    dimension=Fraction(source["dimension_upper_rational"])-1
    num,den=dimension.numerator,dimension.denominator
    bits=source["parameter_radius_bits"]
    radius=Fraction(0) if bits is None else Fraction(1,2**bits)
    p=np.array([[enclosure(Fraction.from_float(v)-radius)[0],
                 enclosure(Fraction.from_float(v)+radius)[1]] for v in box["params"]])
    M,o,_=frame(box)
    scale=np.abs(M).max(axis=2).sum(axis=0)
    pi2,cs,eps=trig_constants()
    stack=[((0,0,0),(0,0,0),True)]
    pos=leaves=extra=0
    largest=0.
    while stack:
        idx,depth,stored=stack.pop()
        axis=3
        if stored:
            axis=tree[pos];pos+=1
        if axis==3:
            rad=np.array([2.**-d for d in depth])
            mid=np.array([(2*i+1)*h-1 for i,h in zip(idx,rad)])
            bound=sufficient(mid,rad,p,M,o,pi2,cs,eps,box["case_order"],num,den)
            if bound<1.:
                leaves+=1;largest=max(largest,bound)
                continue
            extra+=1
            if extra>budget or max(depth)>48:
                raise RuntimeError("Directed dimension refinement limit")
            axis=int(np.argmax(rad*scale))
            stored=False
        if axis not in (0,1,2):
            raise ValueError("Invalid tree")
        for side in (0,1):
            cc=list(idx);dd=list(depth)
            cc[axis]=2*cc[axis]+side;dd[axis]+=1
            stack.append((tuple(cc),tuple(dd),stored))
    if pos!=len(tree):
        raise ValueError("Trailing primary nodes")
    return dict(case=box["case"],parameter_radius_bits=bits,dimension_upper_rational=source["dimension_upper_rational"],
        independent_method="directed binary64 tangent matrix, direct two-by-two minors and products of single-step cofactor matrices",
        recomputed_nodes=len(tree),additional_splits=extra,accepted_leaves=leaves,
        maximum_accepted_power_upper=largest,box_digest=digest)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("box",type=Path);ap.add_argument("primary",type=Path)
    ap.add_argument("--budget",type=int,default=400_000);ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=check(json.loads(a.box.read_text()),json.loads(a.primary.read_text()),a.budget)
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
