"""Independent high-precision finite-separation growth for selected map cases."""
import argparse
import json
from pathlib import Path
import mpmath as mp


def step(x,p):
    r,a,b,c,k,l,m,h,j,n=p
    x,y,z=(r*x[0]+1-r,r*x[1],r*x[2])
    t=a+k*z+h*z**3
    co,si=mp.cos(t),mp.sin(t)
    x,y=co*x-si*y,si*x+co*y
    t=b+l*x+j*x**3
    co,si=mp.cos(t),mp.sin(t)
    y,z=co*y-si*z,si*y+co*z
    t=c+m*y+n*y**3
    co,si=mp.cos(t),mp.sin(t)
    z,x=co*z-si*x,si*z+co*x
    return x,y,z


def probe(row,digits,burn,observed):
    mp.mp.dps=digits
    p=[mp.mpf(float(v)) for v in row["params"]]
    x=tuple(mp.mpf(float(v)) for v in row["primary_start"])
    eps=mp.mpf("1e-25")
    direction=[x/mp.sqrt(6) for x in (1,2,-1)]
    values=[];acc=mp.mpf(0)
    for n in range((burn+observed)//4):
        y=tuple(a+eps*v for a,v in zip(x,direction))
        for _ in range(4):
            x,y=step(x,p),step(y,p)
        delta=[b-a for a,b in zip(x,y)]
        norm=mp.sqrt(sum(b*b for b in delta))
        direction=[b/norm for b in delta]
        done=4*(n+1)
        if done>burn:
            acc+=mp.log(norm/eps)
            if (done-burn)%(observed//4)==0:
                values.append(float(acc/(observed//4)))
                acc=mp.mpf(0)
    return dict(decimal_digits=digits,burn=burn,observed=observed,
                rescale_every=4,separation="1e-25",leading_windows=values,mean=sum(values)/len(values))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    default=Path(__file__).resolve().parents[3]/"results/shape_search_20261121/maps"
    p.add_argument("--results",type=Path,default=default)
    p.add_argument("--cases",type=int,nargs="+",required=True)
    p.add_argument("--digits",type=int,nargs="+",default=[70,95])
    p.add_argument("--burn",type=int,default=8000)
    p.add_argument("--observed",type=int,default=32_000)
    a=p.parse_args()
    if a.burn%4 or a.observed%16 or a.observed<=0:
        raise ValueError("Need four equal windows and four-step rescaling")
    (a.results/"precision").mkdir(parents=True,exist_ok=True)
    for case in a.cases:
        row=json.loads((a.results/"checks"/f"{case}.json").read_text())
        packet=dict(case=case,probes=[])
        for digits in a.digits:
            v=probe(row,digits,a.burn,a.observed);packet["probes"].append(v)
            (a.results/"precision"/f"{case}.json").write_text(json.dumps(packet,indent=2)+"\n")
            print(case,digits,v["mean"],flush=True)
