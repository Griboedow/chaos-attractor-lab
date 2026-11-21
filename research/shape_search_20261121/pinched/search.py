"""Search pinched fractional-linear maps forced by a translation of a two-torus."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numba import njit

OMEGA=(math.sqrt(2)-1,math.sqrt(3)-1)
STARTS=((math.sqrt(5)/7,math.sqrt(7)/9,.81),
        (math.sqrt(11)/7,math.sqrt(13)/9,.011),
        (math.sqrt(17)/7,math.sqrt(19)/9,.42))


def cases(count):
    rng=np.random.default_rng(2026112112)
    modes=np.array([(k,l) for k in range(0,6) for l in range(-5,6)
                    if (k>0 or l>0) and math.gcd(k,abs(l))==1])
    result=[]
    while len(result)<count:
        groups=2 if len(result)%4 else 3
        waves=modes[rng.choice(len(modes),size=groups,replace=False)]
        if all(waves[0,0]*l==waves[0,1]*k for k,l in waves[1:]):
            continue
        factors=np.column_stack((waves,rng.integers(0,32,size=groups)))
        kappa=int(rng.choice([6,8,12,16]))*4**(groups-1)
        # All isometries of the base are kept fixed across the survey.
        result.append({"case":len(result)+1,"waves":factors.tolist(),"kappa":kappa})
    return result


@njit(cache=True)
def forward(waves,kappa,init,burn,n):
    t,u,h=init
    logh=math.log(h)
    out=np.empty((n,3))
    heights=np.empty(n)
    base,physical=0.,0.
    for i in range(n+burn):
        loga=math.log(kappa)
        for k,l,phase in waves:
            angle=k*t+l*u+phase/32
            sine=math.sin(math.pi*(angle-round(angle)))
            if sine==0:
                raise ValueError("An evaluated phase hit a pinch exactly")
            loga+=2*math.log(abs(sine))
        logden=math.log1p((kappa-1)*math.exp(logh))
        if i>=burn:
            base+=loga/n
            physical+=(loga-2*logden)/n
        logh+=loga-logden
        h=math.exp(logh)
        t=(t+OMEGA[0])%1
        u=(u+OMEGA[1])%1
        if i>=burn:
            v=.20+.95*h
            radial=1.32+v*math.cos(2*math.pi*u)
            out[i-burn]=radial*math.cos(2*math.pi*t),radial*math.sin(2*math.pi*t),v*math.sin(2*math.pi*u)
            heights[i-burn]=h
    return out,heights,base,physical


def sample(row,init,burn,n):
    return forward(np.array(row["waves"],dtype=np.int64),row["kappa"],init,burn,n)


def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def screen(count,out,short):
    rows=cases(count)
    desc=[]
    for row in rows:
        xyz,h,base,fiber=sample(row,STARTS[0],2000,short)
        analytic=math.log(row["kappa"]/(4**len(row["waves"])))
        # Summarize shape in fixed physical coordinates. Success or failure of a
        # short orbit's exponent estimate has no bearing on the analytic theorem.
        x=(xyz-xyz.mean(0))
        d=np.linalg.norm(x-x[np.random.default_rng(47).permutation(len(x))],axis=1)
        v=np.r_[np.quantile(d,[.05,.15,.3,.5,.7,.85,.95]),
                np.quantile(h,[.05,.15,.3,.5,.65,.8,.9,.96,.99])]
        desc.append(v)
        row.update(analytic_zero_rate=analytic,observed_zero_rate=base,observed_upper_rate=fiber)
    desc=np.array(desc)
    desc=(desc-desc.mean(0))/np.maximum(.015,desc.std(0))
    nearest=np.full(len(rows),np.inf)
    selected=[]
    idx=int(np.argmax(desc[:,-2]))
    for _ in range(min(36,len(rows))):
        selected.append(idx)
        nearest=np.minimum(nearest,np.linalg.norm(desc-desc[idx],axis=1))
        nearest[selected]=-1
        idx=int(np.argmax(nearest))
    save(out/"screen.json",{"seed":2026112112,"count":count,"steps":short,
                           "selected":[rows[k] for k in selected]})
    print(count,"parameter sets and",len(selected),"shape representatives",flush=True)


def replay(selection,out,steps,ids):
    rows=json.loads(selection.read_text())["selected"]
    if ids:
        rows=[r for r in rows if r["case"] in ids]
        assert len(rows)==len(set(ids))
    out.mkdir(parents=True,exist_ok=True)
    data=[]
    for row in rows:
        runs=[]
        for index,initial in enumerate(STARTS):
            xyz,h,base,fiber=sample(row,initial,8000,steps)
            assert np.all((h>=0)&(h<=1+1e-14)) and np.isfinite([base,fiber]).all()
            if index==0:
                name=f"case-{row['case']}.npz"
                np.savez_compressed(out/name,points=xyz)
            runs.append({"initial":list(initial),"zero_graph_rate":base,"upper_graph_rate":fiber,
                         "height_underflows_in_display":int((h==0).sum()),
                         "height_quantiles":np.quantile(h,[.1,.5,.9]).tolist()})
        data.append({**row,"runs":runs,"steps":steps,"burn":8000,"cloud":name})
    save(out/"manifest.json",data)
    print(len(data),"sampled geometries",flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--count",type=int,default=768)
    ap.add_argument("--short",type=int,default=16000)
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--replay",type=Path)
    ap.add_argument("--steps",type=int,default=240000)
    ap.add_argument("--cases",nargs="+",type=int)
    args=ap.parse_args()
    if args.replay:
        replay(args.replay,args.out,args.steps,args.cases)
    else:
        screen(args.count,args.out,args.short)
