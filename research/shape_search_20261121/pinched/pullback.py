"""Check finite pullbacks by independent high-precision summation at fixed base points."""
import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from search import save


def additive_log(row,q,depth):
    """Forward log-height solver, initialized at height one at the past endpoint."""
    w=(math.sqrt(2)-1,math.sqrt(3)-1)
    t=(float(q[0])-depth*w[0])%1
    u=(float(q[1])-depth*w[1])%1
    h=0.
    for _ in range(depth):
        a=math.log(row["kappa"])
        for k,l,phase in row["waves"]:
            v=k*t+l*u+phase/32
            a+=2*math.log(abs(math.sin(math.pi*(v-round(v)))))
        h+=a-math.log1p((row["kappa"]-1)*math.exp(h))
        t=(t+w[0])%1
        u=(u+w[1])%1
    return h


def exact_series(row,q,depths):
    """Form the finite reciprocal products using absolute phases and guard digits."""
    w=(mp.sqrt(2)-1,mp.sqrt(3)-1)
    kappa=row["kappa"]
    product=mp.mpf(1)
    total=mp.mpf(0)
    answers={}
    for j in range(1,max(depths)+1):
        t,u=q[0]-j*w[0],q[1]-j*w[1]
        a=mp.mpf(kappa)
        for k,l,ph in row["waves"]:
            v=mp.fmod(k*t+l*u+mp.mpf(ph)/32,1)
            a*=mp.sinpi(v)**2
        product/=a
        total+=(kappa-1)*product
        if j in depths:
            # The final summand represents the finite pullback starting at h=1.
            answers[j]=float(-mp.log(total+product))
    return answers


def check(manifest,out):
    mp.mp.prec=180
    rows=json.loads(manifest.read_text())["selected"]
    data=[]
    depths=[256,1024,2048]
    phases=[(mp.mpf(13)/37,mp.mpf(31)/43),(mp.mpf(19)/59,mp.mpf(23)/41)]
    for row in rows:
        checks=[]
        for q in phases:
            answers=exact_series(row,q,depths)
            checks.append({
                "log_heights":[answers[d] for d in depths],
                "absolute_difference_at_2048":abs(answers[2048]-additive_log(row,q,2048)),
                "depth_1024_to_2048_difference":abs(answers[2048]-answers[1024])})
        data.append({"case":row["case"],"checks":checks})
    mx=max(v["absolute_difference_at_2048"] for x in data for v in x["checks"])
    assert mx<1e-4,mx
    save(out,{"precision_bits":mp.mp.prec,"depths":depths,"initial_height":1,
              "phase_samples":[["13/37","31/43"],["19/59","23/41"]],
              "maximum_absolute_error_in_log_height":mx,"records":data})
    print("Largest absolute log-height discrepancy",mx,flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest",type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    check(args.manifest,args.out)
