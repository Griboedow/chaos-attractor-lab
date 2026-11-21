"""Independently replay a telescoping trapping proof with 150-bit Arb matrices."""
import argparse
import base64
import hashlib
import json
import sys
import zlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/"research/shape_search_20261121/maps"))
from flint import arb,arb_mat,ctx
from verify_boxes import column,frames,matrix,rotate,values

def product(V,C,E,M):
    P=V;R=arb_mat(3,3)
    for A,B in zip(reversed(C),reversed(E)):
        R+=P*B;P=P*A
    return R+P*M

def enclosure(c,D,rad):
    return column([x+arb(0,sum(D[i,j].abs_upper()*rad[j] for j in range(3))) for i,x in enumerate(values(c))])

def test(mid,rad,o,M,V,p,steps):
    ident=matrix([[int(i==j) for j in range(3)] for i in range(3)])
    c=o+M*column(mid);w=o+M*column([arb(y,r) for y,r in zip(mid,rad)])
    D=M;Cs=[];Es=[]
    for n in range(steps):
        c,C=rotate(values(c),p,ident)
        w,J=rotate(values(w),p,ident)
        C=arb_mat([[C[i,j].mid() for j in range(3)] for i in range(3)])
        Es.append((J-C)*D);Cs.append(C)
        D=product(ident,Cs,Es,M);extra=enclosure(c,D,rad)
        w=column([a.intersection(b) for a,b in zip(values(w),values(extra))])
        b=values(enclosure(V*(c-o),product(V,Cs,Es,M),rad))
        if n+1==steps:
            return all(abs(x)<1 for x in b)
        if not any(x>1 or x< -1 for x in b):return False
    raise ValueError("Empty return")

def replay(box,proof):
    ctx.prec=150
    if hashlib.sha256(json.dumps(box,sort_keys=True).encode()).hexdigest()!=proof["box_digest"]:raise ValueError("Box digest")
    o,M,V=frames(box);p=[arb(float(x)) for x in box["params"]]
    m=box["case_order"];tree=zlib.decompress(base64.b64decode(proof["partition_zlib_base64"]))
    stack=[([arb(0)]*3,[arb(1)]*3,True,0)];nodes=extra=leaves=0
    while stack:
        mid,rad,stored,depth=stack.pop()
        code=tree[nodes] if stored else 3
        if stored:nodes+=1
        if code==3:
            if test(mid,rad,o[0],M[0],V[0],p,m):
                leaves+=1;continue
            extra+=1
            if extra>30000 or depth>68:raise RuntimeError("Extra Arb subdivisions exhausted")
            code=max(range(3),key=lambda j:float(rad[j]));stored=False
        if code not in (0,1,2):raise ValueError("Invalid tree")
        r=rad.copy();r[code]/=2
        for sign in (-1,1):
            c=mid.copy();c[code]+=r[code]*sign
            stack.append((c,r,stored,depth+1))
    if nodes!=len(tree):raise ValueError("Unvisited nodes")
    return dict(case=box["case"],components=m,precision_bits=ctx.prec,replayed_nodes=nodes,
                accepted_leaves=leaves,extra_splits=extra,all_earlier_phases_disjoint=True,strict_last_return=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("cache",type=Path);ap.add_argument("out",type=Path)
    ap.add_argument("--proof-dir",default="proofs")
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    for f in sorted((a.cache/a.proof_dir).glob("*.json")):
        proof=json.loads(f.read_text());box=json.loads((a.cache/"return-refined"/f.name).read_text())
        record=replay(box,proof)
        out=a.out/f.name;out.write_text(json.dumps(dict(box=box,directed=proof,arb=record),indent=2)+"\n")
        print(f.stem,record["components"],record["accepted_leaves"],flush=True)
