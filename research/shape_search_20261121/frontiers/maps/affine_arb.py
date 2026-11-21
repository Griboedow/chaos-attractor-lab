"""Verify larger parameter blocks by transporting interval state remainders."""
import argparse
import base64
import hashlib
import json
import zlib
from pathlib import Path
from arb_entropy import arb_parameters
from entropy import INPUT,regions
from flint import arb,arb_mat,ctx
from verify_boxes import column,frames,matrix,rotate,values

def centered(v):
    return column([x.mid() for x in values(v)])

def transport(V,C,Q,M,seed,rad):
    P=V;R=arb_mat(3,1)
    for A,b in zip(reversed(C),reversed(Q)):
        R+=P*b;P=P*A
    W=P*M
    R+=P*seed+W*column([arb(0,r) for r in rad])
    return R,W

def test(mid,rad,M,src,dst,V,p,pc,t,nsteps,phase_hulls=None,differentials=None):
    I=matrix([[int(i==j) for j in range(3)] for i in range(3)])
    init=src+M*column(mid);c=centered(init);seed=init-c
    err=seed+M*column([arb(0,r) for r in rad])
    C=[];Q=[];margin=None
    for k in range(nsteps):
        if phase_hulls is not None:
            for j,x in enumerate(values(c+err)):
                before=phase_hulls[k][j]
                phase_hulls[k][j]=x if before is None else before.union(x)
        e=column([v.union(arb(0)) for v in values(err)])
        native,J=rotate(values(c+e),p,I)
        best,Jc=rotate(values(c),pc,I)
        y,_=rotate(values(c),p)
        nn=centered(best);Jc=arb_mat([[Jc[i,j].mid() for j in range(3)] for i in range(3)])
        if differentials is not None:differentials.append((J,Jc,values(c+e),values(c)))
        C.append(Jc);Q.append((J-Jc)*err+(y-nn))
        err,_=transport(I,C,Q,M,seed,rad)
        err=column([a.intersection(b) for a,b in zip(values(err),values(native-nn))]);c=nn
        if differentials is not None:continue
        if k+1<nsteps and not t["name"].startswith("trapping"):continue
        bound,W=transport(V,C,Q,M,seed,rad);projected=values(V*(c-dst)+bound)
        if k+1<nsteps:
            good=max(max((value-1).lower(),(-1-value).lower()) for value in projected);active=range(3)
        else:
            active=[j for j in range(3) if t["target_lo"][j] is not None or t["target_hi"][j] is not None]
            a=[]
            for j in active:
                if t["target_lo"][j] is not None:a.append((projected[j]-t["target_lo"][j]).lower())
                if t["target_hi"][j] is not None:a.append((t["target_hi"][j]-projected[j]).lower())
            good=min(a)
        if not good>0:
            weights=[float(rad[i])*max(float(W[j,i].abs_upper()) for j in active) for i in range(3)]
            return False,None,weights
        margin=good if margin is None else min(margin,good)
    return True,margin,[]

def region(box,t,bits,budget):
    ctx.prec=144
    oo,mm,vv=frames(box);o,M,V=oo[0],mm[0],vv[0]
    hw=[(arb(h)-arb(l))/2 for l,h in zip(t["lo"],t["hi"])]
    center=[(arb(h)+arb(l))/2 for l,h in zip(t["lo"],t["hi"])]
    src=o+M*column(center);A=M*arb_mat([[hw[i]*int(i==j) for j in range(3)] for i in range(3)])
    p=arb_parameters(box,bits);pc=arb_parameters(box,0)
    active=[k for k in range(3) if t["lo"][k]!=t["hi"][k]]
    stack=[([arb(0)]*3,[arb(1)]*3,0)];tree=bytearray();leaves=0;small=None
    while stack:
        mid,rad,depth=stack.pop()
        if depth>66 or len(tree)>budget:raise RuntimeError(f"{t['name']} unresolved at {len(tree)}, depth {depth}")
        ok,margin,w=test(mid,rad,A,src,o,V,p,pc,t,box["case_order"])
        if ok:
            tree.append(3);leaves+=1;small=margin if small is None else min(small,margin);continue
        large=max(float(rad[i]) for i in active)
        split=max((i for i in active if float(rad[i])>=large/64),key=lambda i:w[i])
        tree.append(split);rr=rad.copy();rr[split]/=2
        for sign in (-1,1):
            x=mid.copy();x[split]+=sign*rr[split];stack.append((x,rr,depth+1))
    margin=2.**(int(float(small.log()/arb(2).log())//1)-1)
    if not small>arb(margin):raise ArithmeticError("Incorrect margin rounding")
    return dict(region=t,bits=bits,box_sha256=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest(),
        nodes=len(tree),leaves=leaves,precision_bits=ctx.prec,normalized_strict_margin_above=margin,
        partition_zlib_base64=base64.b64encode(zlib.compress(tree,9)).decode())

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("out",type=Path)
    ap.add_argument("--source",type=Path,default=INPUT);ap.add_argument("--bits",type=int,required=True)
    ap.add_argument("--budget",type=int,default=30000)
    a=ap.parse_args();box=json.loads(a.source.read_text())["box"];a.out.mkdir(parents=True,exist_ok=True)
    for t in regions():
        f=a.out/(t["name"]+".json")
        if f.exists():continue
        r=region(box,t,a.bits,a.budget);f.write_text(json.dumps(r,indent=2)+"\n");print(t["name"],r["nodes"],flush=True)
