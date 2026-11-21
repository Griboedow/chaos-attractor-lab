"""Replay complete region transfers with independently rounded Arb arithmetic."""
import argparse
import base64
import hashlib
import json
import zlib
from pathlib import Path
from entropy import INPUT,regions
from flint import arb,arb_mat,ctx
from replay_traps import product,enclosure
from verify_boxes import column,frames,matrix,rotate,values

def arb_parameters(box,bits):
    return [arb(float(v),2.**-bits) if bits else arb(float(v)) for v in box["params"]]

def capture(mid,rad,src,M,dst,V,p,t,order):
    identity=matrix([[int(i==j) for j in range(3)] for i in range(3)])
    c=src+M*column(mid);w=src+M*column([arb(y,r) for y,r in zip(mid,rad)])
    D=M;Cs=[];Es=[];slack=None;weights=[]
    for n in range(order):
        c,C=rotate(values(c),p,identity);w,J=rotate(values(w),p,identity)
        C=arb_mat([[C[i,j].mid() for j in range(3)] for i in range(3)])
        Es.append((J-C)*D);Cs.append(C)
        D=product(identity,Cs,Es,M);new=enclosure(c,D,rad)
        w=column([a.intersection(b) for a,b in zip(values(w),values(new))])
        if n+1<order and not t["name"].startswith("trapping"):continue
        W=product(V,Cs,Es,M)
        v=values(enclosure(V*(c-dst),W,rad))
        if n+1<order:
            margin=max(max((value-1).lower(),(-1-value).lower()) for value in v);active=range(3)
        else:
            margins=[]
            active=[j for j in range(3) if t["target_lo"][j] is not None or t["target_hi"][j] is not None]
            for j in active:
                if t["target_lo"][j] is not None:margins.append((v[j]-arb(float(t["target_lo"][j]))).lower())
                if t["target_hi"][j] is not None:margins.append((arb(float(t["target_hi"][j]))-v[j]).lower())
            margin=min(margins)
        slack=margin if slack is None else min(slack,margin)
        if not margin>0:
            weights=[max(float(W[j,i].abs_upper())*float(rad[i]) for j in active) for i in range(3)]
            return False,None,weights
    return True,slack,weights

def replay(box,t,record,bits,budget=60000):
    ctx.prec=160
    if record["region"]!=t or record["bits"]!=bits:raise ValueError("Different region or parameter block")
    if record["box_sha256"]!=hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest():raise ValueError("Different coefficient vector or chart")
    tree=zlib.decompress(base64.b64decode(record["partition_zlib_base64"]))
    o,M,V=frames(box);o,M,V=o[0],M[0],V[0]
    small=[(arb(h)-arb(l))/2 for l,h in zip(t["lo"],t["hi"])]
    mid=[(arb(h)+arb(l))/2 for l,h in zip(t["lo"],t["hi"])]
    start=o+M*column(mid);base=M*arb_mat([[small[i]*int(i==j) for j in range(3)] for i in range(3)])
    p=arb_parameters(box,bits);active=[k for k in range(3) if t["lo"][k]!=t["hi"][k]]
    stack=[([arb(0)]*3,[arb(1)]*3,True,0)];nodes=extra=leaves=0;minimum=None
    while stack:
        mid,rad,stored,depth=stack.pop()
        if depth>80:raise RuntimeError("Maximum dyadic depth")
        code=tree[nodes] if stored else 3
        if stored:nodes+=1
        if code==3:
            ok,margin,w=capture(mid,rad,start,base,o,V,p,t,box["case_order"])
            if ok:
                minimum=margin if minimum is None else min(minimum,margin);leaves+=1;continue
            extra+=1
            if extra>budget:raise RuntimeError("Arb supplementary split limit")
            maximum=max(float(rad[i]) for i in active)
            code=max((i for i in active if float(rad[i])>=maximum/64),key=lambda i:w[i]);stored=False
        if code not in active:raise ValueError("Invalid split")
        rr=rad.copy();rr[code]/=2
        for sign in (-1,1):
            cc=mid.copy();cc[code]+=sign*rr[code];stack.append((cc,rr,stored,depth+1))
    if nodes!=len(tree) or nodes!=record["nodes"]:raise ValueError("Partition indices differ")
    lower=2.**(int(float(minimum.log()/arb(2).log())//1)-1)
    if not minimum>arb(lower):raise ArithmeticError("Unverified margin rounding")
    return dict(name=t["name"],stored_nodes=nodes,extra_splits=extra,accepted_leaves=leaves,
        normalized_strict_margin_above=lower)

def run(source,folder,out,bits):
    ctx.prec=160
    box=json.loads(source.read_text())["box"]
    out.mkdir(parents=True,exist_ok=True)
    for t in regions():
        path=out/(t["name"]+".json")
        if path.exists():continue
        record=json.loads((folder/(t["name"]+".json")).read_text())
        proof=replay(box,t,record,bits)
        path.write_text(json.dumps(proof,indent=2)+"\n");print(t["name"],proof["accepted_leaves"],flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("directed",type=Path);ap.add_argument("out",type=Path)
    ap.add_argument("--bits",type=int,default=0);ap.add_argument("--source",type=Path,default=INPUT)
    a=ap.parse_args();run(a.source,a.directed,a.out,a.bits)
