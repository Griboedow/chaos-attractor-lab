"""Look inside each cyclic cloud for a higher-order positive-expansion partition."""
import argparse
import json
from pathlib import Path
from search import CACHE, RESULTS, START, sample
from cyclic import phase_metrics
from numba import njit
import numpy as np

BASES=[(2202255,34),(2603621,29),(2601546,13),(3603818,11),(2303546,8),(2603268,3),
       (20002293,170),(20100964,145),(20202568,78)]


@njit(cache=True)
def residuals(z,m):
    centers=np.zeros((m,3));count=np.zeros(m);sums=np.zeros(m)
    for n in range(len(z)):
        j=n%m;count[j]+=1
        for i in range(3):
            centers[j,i]+=z[n,i]
    for j in range(m):
        for i in range(3):
            centers[j,i]/=count[j]
    for n in range(len(z)):
        j=n%m
        for i in range(3):
            sums[j]+=(z[n,i]-centers[j,i])**2
    return sums


def refinement(z,m):
    total=residuals(z,1)[0];initial=residuals(z,m)
    if total<1e-12 or 1-initial.sum()/total<.65:
        return None
    scores=[]
    for q in (2,3,4,5,6,7,8):
        resid=residuals(z,m*q).reshape(q,m).sum(axis=0)
        relative=1-resid/np.maximum(initial,1e-30)
        scores.append(float(np.median(relative)))
    highest=max(scores)
    if highest<.50:
        return None
    q=next(i+2 for i,s in enumerate(scores) if s>highest-.04)
    return dict(base_order=m,multiplier=q,order=m*q,median_extra_variance=scores[q-2],
                original_variance_explained=float(1-initial.sum()/total))


def local_design(p,group,count):
    rng=np.random.default_rng(np.random.SeedSequence([20261121,990,group]))
    delta=rng.uniform(-1,1,(count,10))
    # Half the lanes vary one active coefficient, half perturb several.
    active=np.flatnonzero(np.r_[np.ones(7,dtype=bool),p[7:10]!=0])
    for i in range(0,count,2):
        j=active[(i//2)%len(active)]
        v=delta[i,j];delta[i,:]=0;delta[i,j]=v
    delta[:,7:10][:,p[7:10]==0]=0
    exponent_range=(-7.5,-1.8) if group>=6 else (-4.8,-1.0)
    scale=10**rng.uniform(exponent_range[0],exponent_range[1],count)
    scales=np.array([.15,1,1,1,5,5,5,10,10,10])
    return p+delta*scale[:,None]*scales


def run(cache,results,group,count,name):
    base,m=BASES[group]
    row=json.loads((results/"checks"/f"{base}.json").read_text())
    params=local_design(np.array(row["params"]),group,count)
    rows=[];positive=0
    for begin in range(0,count,64):
        p=params[begin:begin+64]
        protocol=(48_000,65_536,8,20_480) if group>=6 else (12_000,16_384,8,8192)
        w,z=sample(p,np.tile(START,(len(p),1)),*protocol)
        cutoff=.00025 if group>=6 else .002
        for i in np.flatnonzero(w.min(axis=1)>cutoff):
            positive+=1
            refined=refinement(z[i],m)
            if refined:
                rows.append(dict(case=20_000_000+group*100_000+begin+int(i),params=p[i].tolist(),
                                 rates=w[i].tolist(),**refined))
    out=cache/name;out.mkdir(parents=True,exist_ok=True)
    (out/f"screen-{group}.json").write_text(json.dumps(dict(base=base,parameters=count,positive=positive,candidates=rows))+"\n")
    print(base,count,positive,len(rows),max([r["order"] for r in rows],default=0),flush=True)


def continue_search(cache,name,deep):
    out=cache/name
    rows=[]
    for f in sorted(out.glob("screen-*.json")):
        data=json.loads(f.read_text())
        for q in range(2,9):
            group=[r for r in data["candidates"] if r["multiplier"]==q]
            group.sort(key=lambda r:-(r["median_extra_variance"]*np.mean(r["rates"])))
            # Include both the clearest and intermediates on the diagnostic scale.
            rows.extend(group[:4])
            if len(group)>4:
                rows.extend(group[j] for j in np.linspace(4,len(group)-1,min(4,len(group)-4),dtype=int))
    retained=[]
    for begin in range(0,len(rows),12):
        batch=rows[begin:begin+12];params=np.array([r["params"] for r in batch])
        protocol=(200_000,1_200_000,16,200_000) if deep else (100_000,400_000,16,100_000)
        w,z=sample(params,np.tile(START,(len(params),1)),*protocol)
        for r,ww,zz in zip(batch,w,z):
            ref=refinement(zz,r["base_order"])
            if ww.min()>(.00012 if deep else .0008) and ref and ref["order"]==r["order"]:
                info=phase_metrics(zz,r["order"])
                retained.append(dict(case=r["case"],params=r["params"],rates=ww.tolist(),**ref,
                    minimum_center_95radius_ratio=info["minimum_center_95radius_ratio"],
                    median_nearest_ratio=info["median_nearest_ratio"]))
    (out/"preview.json").write_text(json.dumps(retained)+"\n")
    print("Long positive continuations",len(retained),"out of",len(rows),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE);p.add_argument("--results",type=Path,default=RESULTS)
    p.add_argument("--count",type=int,default=3072)
    p.add_argument("--groups",type=int,nargs="+",default=list(range(6)))
    p.add_argument("--name",default="nested-phase")
    p.add_argument("--continue-only",action="store_true")
    a=p.parse_args()
    if not a.continue_only:
        for g in a.groups:
            run(a.cache,a.results,g,a.count,a.name)
    continue_search(a.cache,a.name,min(a.groups)>=6)
