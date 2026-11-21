"""Directed-interval replay of the derivative and exterior-power upper bound."""
import argparse
import base64
import json
import zlib
from fractions import Fraction
from pathlib import Path
import numpy as np
from numba import njit
from affine_directed import add,central,sub,transport
from entropy import INPUT,parameters
from directed_return import frame,maximum_abs,minus,overlap,plus,rotate,times,transform,trig_constants
from telescoped_return import product

@njit(cache=True)
def cof(A):
    C=np.empty((3,3,2))
    for i in range(3):
        for j in range(3):
            C[i,j]=minus(times(A[(i+1)%3,(j+1)%3],A[(i+2)%3,(j+2)%3]),times(A[(i+1)%3,(j+2)%3],A[(i+2)%3,(j+1)%3]))
    return C

@njit(cache=True)
def cell(mid,rad,M,src,p,pc,steps,pi2,cs,eps,num,den):
    I=np.zeros((3,3,2));u=np.empty((3,1,2));r=u.copy()
    for i in range(3):I[i,i]=(1.,1.);u[i,0]=(mid[i],mid[i]);r[i,0]=(-rad[i],rad[i])
    init=add(src,transform(M,u));c=central(init);seed=sub(init,c);err=add(seed,transform(M,r))
    Js=np.empty((steps,3,3,2));Es=Js.copy();Cs=Js.copy();Bs=Js.copy();Q=np.empty((steps,3,1,2))
    D=I.copy();A=I.copy()
    for k in range(steps):
        local=np.empty((3,4,2));local[:,1:,:]=I;at=local.copy();at[:,0,:]=c[:,0,:]
        for i in range(3):local[i,0]=plus(c[i,0],(min(0.,err[i,0,0]),max(0.,err[i,0,1])))
        local=rotate(local,p,pi2,cs,eps,True);at=rotate(at,pc,pi2,cs,eps,True);yc=rotate(c,p,pi2,cs,eps,True)
        nn=central(at[:,:1,:]);Js[k]=central(at[:,1:,:])
        J=local[:,1:,:];dev=sub(J,Js[k])
        Q[k]=add(transform(dev,err),sub(yc,nn))
        err,_=transport(I,Js,Q,k,M,seed,r)
        direct=sub(local[:,:1,:],nn)
        for i in range(3):err[i,0]=overlap(err[i,0],direct[i,0])
        c=nn
        Es[k]=transform(dev,D);D=product(I,Js,Es,k,I)
        Cs[k]=central(cof(Js[k]));Bs[k]=transform(sub(cof(J),Cs[k]),A);A=product(I,Cs,Bs,k,I)
    alt=cof(D);norm=0.;area=0.
    for i in range(3):
        for j in range(3):
            a=maximum_abs(D[i,j]);b=maximum_abs(overlap(A[i,j],alt[i,j]))
            norm=plus((norm,norm),times((a,a),(a,a)))[1];area=plus((area,area),times((b,b),(b,b)))[1]
    score=(1.,1.)
    for i in range(den-num):score=times(score,(norm,norm))
    for i in range(num):score=times(score,(area,area))
    return score[1]

def replay(box,record,budget=150000):
    bits=record["block_radius_bits"];m=box["case_order"]
    if record["map_steps"]!=m or record["case"]!=box["case"]:raise ValueError("Wrong model")
    a=Fraction(record["dimension_upper"])-1;num,den=a.numerator,a.denominator
    if not 0<num<den:raise ValueError("Invalid singular-value dimension")
    M,o,_=frame(box);p=parameters(box,bits);pc=parameters(box,0);pi2,cs,eps=trig_constants()
    dim=np.abs(M).max(axis=2).sum(axis=0)
    tree=zlib.decompress(base64.b64decode(record["partition_zlib_base64"]))
    stack=[((0,0,0),(0,0,0),True)];nodes=extra=leaves=0;largest=0.
    while stack:
        ids,ds,stored=stack.pop()
        if max(ds)>46 or extra>budget:raise RuntimeError("Additional directed splits exhausted")
        rad=np.array([2.**-x for x in ds]);mid=np.array([(2*i+1)*r-1 for i,r in zip(ids,rad)])
        split=tree[nodes] if stored else 3
        if stored:nodes+=1
        if split==3:
            s=cell(mid,rad,M,o[:,None,:],p,pc,m,pi2,cs,eps,num,den)
            if s<1:
                largest=max(largest,s);leaves+=1;continue
            extra+=1;split=int(np.argmax(dim*rad));stored=False
        if split not in (0,1,2):raise ValueError("Illegal subdivision")
        for side in (0,1):
            d=list(ds);i=list(ids);d[split]+=1;i[split]=2*i[split]+side;stack.append((tuple(i),tuple(d),stored))
    if nodes!=len(tree) or nodes!=record["nodes"]:raise ValueError("Incomplete subdivision")
    return dict(case=box["case"],block_radius_bits=bits,dimension_upper=record["dimension_upper"],
        stored_nodes=nodes,additional_splits=extra,accepted_leaves=leaves,maximum_power_upper=largest)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("primary",type=Path);ap.add_argument("output",type=Path)
    ap.add_argument("--source",type=Path,default=INPUT);a=ap.parse_args()
    r=replay(json.loads(a.source.read_text())["box"],json.loads(a.primary.read_text()))
    a.output.write_text(json.dumps(r,indent=2)+"\n");print(r["dimension_upper"],r["accepted_leaves"],flush=True)
