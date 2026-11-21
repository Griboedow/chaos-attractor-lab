"""Replay affine-remainder trees with directed binary64 intervals and rational trigonometry."""
import argparse
import base64
import hashlib
import json
import zlib
from pathlib import Path
import numpy as np
from numba import njit
from entropy import INPUT,regions,parameters
from directed_cover import region_frame
from directed_return import maximum_abs,minus,overlap,plus,rotate,transform,trig_constants

@njit(cache=True)
def add(A,B):
    out=A.copy()
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            out[i,j]=plus(A[i,j],B[i,j])
    return out

@njit(cache=True)
def sub(A,B):
    out=A.copy()
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            out[i,j]=minus(A[i,j],B[i,j])
    return out

@njit(cache=True)
def central(A):
    out=A.copy()
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            x=A[i,j,0]/2+A[i,j,1]/2;out[i,j]=(x,x)
    return out

@njit(cache=True)
def transport(V,C,Q,n,M,seed,r):
    P=V.copy();R=np.zeros((3,1,2))
    for i in range(n,-1,-1):
        R=add(R,transform(P,Q[i]));P=transform(P,C[i])
    W=transform(P,M)
    return add(add(R,transform(P,seed)),transform(W,r)),W

@njit(cache=True)
def evaluate(mid,rad,M,src,dst,V,p,pc,m,pi2,co,eps,lower,upper,exclude):
    I=np.zeros((3,3,2));u=np.empty((3,1,2));r=u.copy()
    for i in range(3):
        I[i,i]=(1.,1.);u[i,0]=(mid[i],mid[i]);r[i,0]=(-rad[i],rad[i])
    init=add(src,transform(M,u));c=central(init);seed=sub(init,c)
    err=add(seed,transform(M,r))
    C=np.empty((m,3,3,2));Q=np.empty((m,3,1,2));minimum=np.inf
    for k in range(m):
        local=np.empty((3,4,2));local[:,1:,:]=I;at=local.copy()
        at[:,0,:]=c[:,0,:]
        for i in range(3):
            local[i,0]=plus(c[i,0],(min(0.,err[i,0,0]),max(0.,err[i,0,1])))
        local=rotate(local,p,pi2,co,eps,True)
        at=rotate(at,pc,pi2,co,eps,True)
        yc=rotate(c,p,pi2,co,eps,True)
        nn=central(at[:,:1,:]);C[k]=central(at[:,1:,:])
        Q[k]=add(transform(sub(local[:,1:,:],C[k]),err),sub(yc,nn))
        err,_=transport(I,C,Q,k,M,seed,r)
        direct=sub(local[:,:1,:],nn)
        for j in range(3):err[j,0]=overlap(err[j,0],direct[j,0])
        c=nn
        if k+1<m and not exclude:continue
        v=transform(V,sub(c,dst));b,W=transport(V,C,Q,k,M,seed,r);b=add(b,v)
        if k+1<m:
            slack=-np.inf
            for j in range(3):slack=max(slack,minus(b[j,0],(1.,1.))[0],minus((-1.,-1.),b[j,0])[0])
        else:
            slack=np.inf
            for j in range(3):
                if np.isfinite(lower[j]):slack=min(slack,minus(b[j,0],(lower[j],lower[j]))[0])
                if np.isfinite(upper[j]):slack=min(slack,minus((upper[j],upper[j]),b[j,0])[0])
        if not slack>0:
            w=np.zeros(3)
            for i in range(3):
                for j in range(3):
                    if k+1<m or np.isfinite(lower[j]) or np.isfinite(upper[j]):
                        w[i]=max(w[i],rad[i]*maximum_abs(W[j,i]))
            return False,0.,w
        minimum=min(minimum,slack)
    return True,minimum,np.zeros(3)

def replay(box,t,record,bits):
    if record["region"]!=t or record["bits"]!=bits or record["box_sha256"]!=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest():raise ValueError("Certificate mismatch")
    tree=zlib.decompress(base64.b64decode(record["partition_zlib_base64"]))
    M,src,dst,V=region_frame(box,t["lo"],t["hi"])
    pc=parameters(box,0);p=parameters(box,bits);pi2,co,eps=trig_constants()
    l=np.array([-np.inf if x is None else x for x in t["target_lo"]])
    h=np.array([np.inf if x is None else x for x in t["target_hi"]])
    active=np.array(t["lo"])!=np.array(t["hi"])
    stack=[((0,0,0),(0,0,0),True)];nodes=extra=leaves=0;slack=np.inf
    while stack:
        ids,ds,stored=stack.pop()
        if max(ds)>46 or extra>150000:raise RuntimeError("Directed replay split budget")
        rad=np.array([2.**-x for x in ds]);mid=np.array([(2*i+1)*r-1 for i,r in zip(ids,rad)])
        code=tree[nodes] if stored else 3
        if stored:nodes+=1
        if code==3:
            ok,margin,w=evaluate(mid,rad,M,src[:,None,:],dst[:,None,:],V,p,pc,box["case_order"],pi2,co,eps,l,h,t["name"].startswith("trapping"))
            if ok:
                leaves+=1;slack=min(slack,margin);continue
            extra+=1
            available=np.flatnonzero(active&(rad>=rad[active].max()/64));code=int(available[np.argmax(w[available])]);stored=False
        if code not in (0,1,2) or not active[code]:raise ValueError("Illegal split")
        for side in (0,1):
            d=list(ds);i=list(ids);d[code]+=1;i[code]=2*i[code]+side;stack.append((tuple(i),tuple(d),stored))
    if nodes!=len(tree) or nodes!=record["nodes"]:raise ValueError("Inconsistent tree size")
    return dict(name=t["name"],original_nodes=nodes,extra_splits=extra,leaves=leaves,
        normalized_margin_lower_directed=slack)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("primary",type=Path);ap.add_argument("out",type=Path)
    ap.add_argument("--bits",type=int,required=True);ap.add_argument("--source",type=Path,default=INPUT)
    a=ap.parse_args();box=json.loads(a.source.read_text())["box"];a.out.mkdir(parents=True,exist_ok=True)
    for t in regions():
        f=a.out/(t["name"]+".json")
        if f.exists():continue
        record=json.loads((a.primary/f.name).read_text());r=replay(box,t,record,a.bits)
        f.write_text(json.dumps(r,indent=2)+"\n");print(t["name"],r["leaves"],flush=True)
