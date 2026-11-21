"""Search one- and two-term real Fourier drives on an integer coefficient grid."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from search import basis,batch


def sparse(hmax):
    choices=(-4,-3,-2,-1,1,2,3,4)
    rows=[]
    for n in (1,2):
        for positions in itertools.combinations(range(1,1+2*hmax),n):
            for amplitudes in itertools.product(choices,repeat=n):
                for bias in (-1,0,1):
                    c=np.zeros(1+2*hmax,int);c[0]=bias
                    for i,a in zip(positions,amplitudes):c[i]=a
                    rows.append(c)
    coeff=np.stack(rows)
    best={}
    total=0;hit=0
    for k in range(3,21):
        real=coeff@basis(np.arange(2*k)*np.pi/k,hmax).T
        m=np.min(abs(abs(real)-2),axis=1)
        out=batch(real)
        idx=np.flatnonzero((out[:,0]==0)&(out[:,1]!=0)&(m>.03))
        total+=len(coeff);hit+=len(idx)
        for i in idx:
            c=coeff[i];mu=int(out[i,1])
            highest=max([0]+[j for j in range(1,hmax+1) if c[j] or c[hmax+j]])
            score=float(m[i]/(1+.10*k+.1*highest))
            key=(abs(mu),highest)
            if key not in best or score>best[key]["score"]:
                packed=np.r_[c[:1+highest],c[1+hmax:1+hmax+highest]]
                best[key]=dict(case=int(47000000+(k-3)*len(coeff)+i),h=highest,k=k,
                               coeffs=packed.tolist(),divisor=1,imag_amplitude=3,mu=mu,
                               crossings=int(out[i,2]),score=score,sampled_threshold_margin=float(m[i]))
    return dict(draws=total,nonzero_mu_screen_hits=hit,best=list(best.values()))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args();d=sparse(7);a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2)+"\n");print(d["draws"],d["nonzero_mu_screen_hits"],len(d["best"]),flush=True)
