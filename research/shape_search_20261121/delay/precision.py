"""Estimate three exponents in MPFR with independently evolved finite perturbations."""
import argparse
import json
from pathlib import Path

import gmpy2 as g

from search import save


def measure(par,burn,blocks,width):
    g.get_context().precision=224
    M=g.mpfr
    r,a,b,c=map(lambda x:M(str(x)),par[:4])
    s=M(str(par[4])) if len(par)==5 else None
    def f(q):
        x,y,z=q
        e=g.exp(r*(1-a*x-b*y-c*z))
        return ((s*x+(1-s)*z*e) if s is not None else x*e,x,y)
    q=tuple(map(M,[".59","1.12",".47"]))
    for _ in range(burn):
        q=f(q)
    eps=M("1e-25")
    basis=[[M(i==j) for j in range(3)] for i in range(3)]
    result=[]
    for block in range(blocks+1):
        accum=[M(0) for _ in range(3)]
        for _ in range(width):
            nextq=f(q)
            new=[]
            for k in range(3):
                w=f([q[i]+eps*basis[k][i] for i in range(3)])
                v=[(w[i]-nextq[i])/eps for i in range(3)]
                for orth in new:
                    dot=sum(v[i]*orth[i] for i in range(3))
                    v=[v[i]-dot*orth[i] for i in range(3)]
                norm=g.sqrt(sum(t*t for t in v))
                accum[k]+=g.log(norm)
                new.append([t/norm for t in v])
            q,basis=nextq,new
        if block:
            result.append([float(t/width) for t in accum])
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest",type=Path)
    ap.add_argument("--cases",type=int,nargs="+",required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--width",type=int,default=20000)
    args=ap.parse_args()
    cases=json.loads(args.manifest.read_text(encoding="utf-8"))
    rows=[]
    for ident in args.cases:
        case=next(r for r in cases if r["case"]==ident)
        windows=measure(case["parameters"],10000,5,args.width)
        rows.append({"case":ident,"windows":windows,"precision_bits":224,"burn":10000,
                     "tangent_burn":args.width,"window_length":args.width,
                     "perturbation":"1e-25","parameter_input":"Decimal strings of saved parameters."})
        save(args.out,rows)
        print(ident,[min(w[k] for w in windows) for k in range(3)],flush=True)
