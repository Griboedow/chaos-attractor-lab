"""Verify all phases with a centered expansion of the tangent-matrix product."""
import argparse
import base64
import hashlib
import json
import sys
import zlib
from pathlib import Path
PARENT=Path(__file__).resolve().parents[2]/"maps"
sys.path.insert(0,str(PARENT))
from directed_return import frame,plus,minus,times,overlap,maximum_abs,transform,rotate,trig_constants
import numpy as np
from numba import njit


@njit(cache=True)
def product(V,C,E,last,M):
    P=V.copy()
    result=np.zeros((3,3,2))
    for i in range(last,-1,-1):
        add=transform(P,E[i])
        for a in range(3):
            for b in range(3):
                result[a,b]=plus(result[a,b],add[a,b])
        P=transform(P,C[i])
    base=transform(P,M)
    for a in range(3):
        for b in range(3):
            result[a,b]=plus(result[a,b],base[a,b])
    return result


@njit(cache=True)
def mean(center,J,rad):
    result=np.empty((3,2))
    for i in range(3):
        width=0.
        for j in range(3):
            a=maximum_abs(J[i,j])
            width=plus((width,width),times((a,a),(rad[j],rad[j])))[1]
        result[i]=plus(center[i,0],(-width,width))
    return result


@njit(cache=True)
def test(mid,rad,M,o,V,p,period,pi2,cs,eps):
    I=np.zeros((3,3,2))
    for i in range(3):
        I[i,i]=(1.,1.)
    center=np.empty((3,1,2));wide=center.copy()
    for i in range(3):
        center[i,0]=(mid[i],mid[i]);wide[i,0]=plus(center[i,0],(-rad[i],rad[i]))
    center=transform(M,center);wide=transform(M,wide)
    for i in range(3):
        center[i,0]=plus(center[i,0],o[i]);wide[i,0]=plus(wide[i,0],o[i])
    C=np.empty((period,3,3,2));E=C.copy()
    D=M.copy()
    for n in range(period):
        local=np.empty((3,4,2));middle=local.copy()
        local[:,0]=wide[:,0];middle[:,0]=center[:,0]
        local[:,1:]=I;middle[:,1:]=I
        local=rotate(local,p,pi2,cs,eps,True)
        middle=rotate(middle,p,pi2,cs,eps,True)
        center=middle[:,:1,:].copy()
        dev=np.empty((3,3,2))
        for i in range(3):
            for j in range(3):
                v=middle[i,j+1,0]/2+middle[i,j+1,1]/2
                C[n,i,j]=(v,v);dev[i,j]=minus(local[i,j+1],C[n,i,j])
        E[n]=transform(dev,D)
        D=product(I,C,E,n,M)
        state=mean(center,D,rad)
        for j in range(3):
            wide[j,0]=overlap(local[j,0],state[j])
        delta=center.copy()
        for j in range(3):
            delta[j,0]=minus(delta[j,0],o[j])
        v=transform(V,delta)
        weights=product(V,C,E,n,M)
        bounds=mean(v,weights,rad)
        if n+1<period:
            passed=False
            for axis in range(3):
                if bounds[axis,0]>1 or bounds[axis,1]<-1:
                    passed=True
            if passed:
                continue
        else:
            passed=True
            for axis in range(3):
                if bounds[axis,0]<=-1 or bounds[axis,1]>=1:
                    passed=False
                if v[axis,0,0]>1 or v[axis,0,1]<-1:
                    return -1,np.zeros(3)
            if passed:
                return 1,np.zeros(3)
        wv=np.zeros(3)
        for j in range(3):
            for i in range(3):
                wv[j]=max(wv[j],rad[j]*maximum_abs(weights[i,j]))
        return 0,wv
    return 1,np.zeros(3)


def run(box,budget):
    M,o,V=frame(box)
    p=np.array([[v,v] for v in box["params"]])
    pi2,cs,eps=trig_constants()
    m=box["case_order"];stack=[((0,0,0),(0,0,0))]
    nodes=leaves=0;encoding=bytearray()
    while stack:
        idx,depth=stack.pop();nodes+=1
        if nodes>budget or max(depth)>45:
            raise RuntimeError(f"subdivision budget or depth reached at {nodes}: {depth}")
        rad=np.array([2.**-d for d in depth])
        mid=np.array([(2*i+1)*r-1 for i,r in zip(idx,rad)])
        status,weight=test(mid,rad,M,o,V,p,m,pi2,cs,eps)
        if status==1:
            leaves+=1;encoding.append(3);continue
        if status<0:
            raise RuntimeError("A center witnesses failed strict return")
        eligible=np.where(rad>=rad.max()/64)[0]
        split=int(eligible[np.argmax(weight[eligible])]);encoding.append(split)
        for side in (0,1):
            t=list(idx);d=list(depth);t[split]=2*t[split]+side;d[split]+=1
            stack.append((tuple(t),tuple(d)))
    return dict(case=box["case"],map_order=m,nodes=nodes,leaves=leaves,
        box_digest=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest(),
        partition_zlib_base64=base64.b64encode(zlib.compress(encoding,9)).decode("ascii"),
        method="directed IEEE interval arithmetic, centered telescoping tangent products",
        complete_avoidance_for_every_earlier_iterate=True,strict_full_return=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("boxes",nargs="+",type=Path)
    ap.add_argument("--output",type=Path,required=True);ap.add_argument("--budget",type=int,default=12000)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    for f in a.boxes:
        try:
            b=json.loads(f.read_text());r=run(b,a.budget)
            (a.output/(f.stem+".json")).write_text(json.dumps(r,indent=2)+"\n")
            print(f.stem,b["case_order"],"proved",r["nodes"],flush=True)
        except RuntimeError as exc:
            print(f.stem,str(exc),flush=True)
