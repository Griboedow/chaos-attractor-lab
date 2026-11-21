"""Certify interval transfers between closed boxes for a two-branch return."""
import argparse
import base64
import hashlib
import json
import sys
import zlib
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
BASE=ROOT/"research/shape_search_20261121/maps"
EXT=ROOT/"research/shape_search_20261121/extension/maps"
sys.path[:0]=[str(BASE),str(EXT)]
import numpy as np
from numba import njit
from directed_cover import region_frame
from directed_return import enclosure,maximum_abs,minus,overlap,plus,rotate,transform,trig_constants
from telescoped_return import product,mean

INPUT=ROOT/"results/shape_search_20261121/extension/maps/traps/2300491-0.json"
FACES=[-.5625,-.25,-.1875,.28125]
SIDES=[-1,1,1,-1]
# Strict intermediate transfer boxes in the original normalized trapping frame.
BRIDGES=[
    [[-.20,-.355,-.417],[-.147,-.283,-.308]],
    [[-.018,-.096,-.940],[.082,.094,-.787]],
    [[-.010,-.107,-.961],[.101,.097,-.797]],
    [[-.122,-.723,-.461],[-.069,-.590,-.296]],
]

def regions(faces=FACES,bridges=BRIDGES,axis=2):
    tests=[]
    for i,(face,(bridge_lo,bridge_hi),side) in enumerate(zip(faces,bridges,SIDES)):
        lo=[-1.]*3;hi=[1.]*3;lo[axis]=hi[axis]=face
        tests.append(dict(name=f"face-{i}-first-return",lo=lo,hi=hi,target_lo=bridge_lo,target_hi=bridge_hi))
        lower=[None]*3;upper=[None]*3
        if side<0:upper[axis]=faces[0]
        else:lower[axis]=faces[-1]
        tests.append(dict(name=f"face-{i}-second-return",lo=bridge_lo,hi=bridge_hi,target_lo=lower,target_hi=upper))
    tests.append(dict(name="trapping-and-prior-phase-exclusions",lo=[-1.]*3,hi=[1.]*3,target_lo=[-1.]*3,target_hi=[1.]*3))
    return tests

def parameters(box,bits):
    radius=Fraction(1,2**bits) if bits else Fraction(0)
    return np.array([[enclosure(Fraction.from_float(x)-radius)[0],enclosure(Fraction.from_float(x)+radius)[1]] for x in box["params"]])

@njit(cache=True)
def evaluate(mid,rad,M,o,t,V,p,period,pi2,cs,eps,lower,upper,exclusions):
    I=np.zeros((3,3,2))
    for i in range(3):I[i,i]=(1.,1.)
    center=np.empty((3,1,2));wide=center.copy()
    for i in range(3):
        center[i,0]=(mid[i],mid[i]);wide[i,0]=plus(center[i,0],(-rad[i],rad[i]))
    center=transform(M,center);wide=transform(M,wide)
    for i in range(3):
        center[i,0]=plus(center[i,0],o[i]);wide[i,0]=plus(wide[i,0],o[i])
    C=np.empty((period,3,3,2));E=C.copy();D=M.copy()
    for n in range(period):
        local=np.empty((3,4,2));middle=local.copy()
        local[:,0]=wide[:,0];middle[:,0]=center[:,0];local[:,1:]=I;middle[:,1:]=I
        local=rotate(local,p,pi2,cs,eps,True);middle=rotate(middle,p,pi2,cs,eps,True)
        center=middle[:,:1,:].copy()
        dev=np.empty((3,3,2))
        for i in range(3):
            for j in range(3):
                v=middle[i,j+1,0]/2+middle[i,j+1,1]/2
                C[n,i,j]=(v,v);dev[i,j]=minus(local[i,j+1],C[n,i,j])
        E[n]=transform(dev,D);D=product(I,C,E,n,M)
        state=mean(center,D,rad)
        for j in range(3):wide[j,0]=overlap(local[j,0],state[j])
        if n+1<period and not exclusions:continue
        delta=center.copy()
        for j in range(3):delta[j,0]=minus(delta[j,0],t[j])
        v=transform(V,delta);W=product(V,C,E,n,M);b=mean(v,W,rad)
        if n+1<period:
            if any_away(b):continue
        else:
            good=True
            for j in range(3):
                if not (b[j,0]>lower[j] and b[j,1]<upper[j]):good=False
            if good:return 1,np.zeros(3)
        w=np.zeros(3)
        for i in range(3):
            for j in range(3):
                if n+1<period or lower[j]>-np.inf or upper[j]<np.inf:
                    w[i]=max(w[i],rad[i]*maximum_abs(W[j,i]))
        return 0,w
    return 0,np.zeros(3)

@njit(cache=True)
def any_away(b):
    for j in range(3):
        if b[j,0]>1 or b[j,1]<-1:return True
    return False

def directed(box,t,bits,budget):
    M,o,target,V=region_frame(box,t["lo"],t["hi"])
    p=parameters(box,bits);pi2,cs,eps=trig_constants()
    lower=np.array([-np.inf if x is None else x for x in t["target_lo"]])
    upper=np.array([np.inf if x is None else x for x in t["target_hi"]])
    active=np.array(t["lo"])!=np.array(t["hi"])
    stack=[((0,0,0),(0,0,0))];part=bytearray();leaves=0
    while stack:
        ids,ds=stack.pop()
        if len(part)>=budget or max(ds)>46:raise RuntimeError(f"{t['name']} budget {len(part)}, depths {ds}")
        rad=np.array([2.**-x for x in ds]);mid=np.array([(2*i+1)*r-1 for i,r in zip(ids,rad)])
        good,w=evaluate(mid,rad,M,o,target,V,p,box["case_order"],pi2,cs,eps,lower,upper,t["name"].startswith("trapping"))
        if good:
            part.append(3);leaves+=1;continue
        eligible=np.flatnonzero(active&(rad>=rad[active].max()/64))
        s=int(eligible[np.argmax(w[eligible])]);part.append(s)
        for side in (0,1):
            d=list(ds);i=list(ids);d[s]+=1;i[s]=2*i[s]+side
            stack.append((tuple(i),tuple(d)))
    return dict(nodes=len(part),leaves=leaves,partition_zlib_base64=base64.b64encode(zlib.compress(part,9)).decode())

def generate(source,out,bits,budget):
    box=json.loads(source.read_text())["box"];out.mkdir(parents=True,exist_ok=True)
    for t in regions():
        path=out/(t["name"]+".json")
        if path.exists():continue
        record=directed(box,t,bits,budget)
        path.write_text(json.dumps(dict(region=t,box_sha256=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest(),bits=bits,**record),indent=2)+"\n")
        print(t["name"],record["nodes"],flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("output",type=Path)
    ap.add_argument("--source",type=Path,default=INPUT);ap.add_argument("--bits",type=int,default=0)
    ap.add_argument("--budget",type=int,default=350000)
    a=ap.parse_args();generate(a.source,a.output,a.bits,a.budget)
