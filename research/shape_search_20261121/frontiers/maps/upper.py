"""Certify Hausdorff upper bounds using recentered orbit and exterior derivatives."""
import argparse
import base64
import json
import zlib
from fractions import Fraction
from pathlib import Path
from affine_arb import arb_parameters,test
from entropy import INPUT
from flint import arb,arb_mat,ctx
from replay_traps import product
from verify_boxes import frames,matrix
from dimension import cofactors,exterior_step,frobenius_squared

def evaluate(mid,rad,M,o,V,p,pc,m,num,den):
    ident=matrix([[int(i==j) for j in range(3)] for i in range(3)])
    items=[]
    test(mid,rad,M,o,o,V,p,pc,{},m,differentials=items)
    Js=[];Es=[];Cs=[];Bs=[];D=ident;A=ident
    for J,Jc,x,xc in items:
        Es.append((J-Jc)*D);Js.append(Jc);D=product(ident,Js,Es,ident)
        Tc=cofactors(Jc)
        T=exterior_step(x,p,ident)
        byminors=cofactors(J);T=arb_mat([[T[i,j].intersection(byminors[i,j]) for j in range(3)] for i in range(3)])
        Bs.append((T-Tc)*A);Cs.append(Tc);A=product(ident,Cs,Bs,ident)
    alt=cofactors(D);A=arb_mat([[A[i,j].intersection(alt[i,j]) for j in range(3)] for i in range(3)])
    return frobenius_squared(D)**(den-num)*frobenius_squared(A)**num

def verify(box,bits,s,budget):
    ctx.prec=144
    frac=Fraction(s)-1;num,den=frac.numerator,frac.denominator
    if not 0<num<den:raise ValueError("Dimension between one and two required")
    p=arb_parameters(box,bits);pc=arb_parameters(box,0);oo,mm,vv=frames(box);o,M,V=oo[0],mm[0],vv[0]
    dims=[sum(float(M[k,j].abs_upper()) for k in range(3)) for j in range(3)]
    stack=[([arb(0)]*3,[arb(1)]*3,0)];tree=bytearray();leaves=0;largest=arb(0)
    while stack:
        mid,rad,depth=stack.pop()
        if len(tree)>budget or depth>100:raise RuntimeError(f"Exterior derivative unresolved at {len(tree)} nodes, depth {depth}")
        upper=evaluate(mid,rad,M,o,V,p,pc,box["case_order"],num,den)
        if upper<1:
            leaves+=1;tree.append(3);largest=max(largest,upper.upper());continue
        split=max(range(3),key=lambda j:float(rad[j])*dims[j]);tree.append(split)
        rr=rad.copy();rr[split]/=2
        for sign in (-1,1):
            cc=mid.copy();cc[split]+=sign*rr[split];stack.append((cc,rr,depth+1))
    return dict(case=box["case"],block_radius_bits=bits,map_steps=box["case_order"],dimension_upper=str(s),
        precision_bits=ctx.prec,nodes=len(tree),leaves=leaves,maximum_accepted_power_approx=float(largest),
        partition_zlib_base64=base64.b64encode(zlib.compress(tree,9)).decode())

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("out",type=Path)
    ap.add_argument("--source",type=Path,default=INPUT);ap.add_argument("--bits",type=int,required=True)
    ap.add_argument("--dimension",default="17/10");ap.add_argument("--budget",type=int,default=40000)
    a=ap.parse_args();r=verify(json.loads(a.source.read_text())["box"],a.bits,a.dimension,a.budget)
    a.out.write_text(json.dumps(r,indent=2)+"\n");print(r["dimension_upper"],r["nodes"],flush=True)
