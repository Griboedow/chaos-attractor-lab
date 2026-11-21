"""Explore a positive delay map with weighted density feedback on three successive samples."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
from numba import njit

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"research"))
from compare_attractor_shapes import fingerprint


@njit(cache=True)
def orbit(par, initial, burn, windows, length, retain):
    r,a,b,c = par[:4]
    survival = par[4] if len(par)==5 else -1.
    x,y,z = np.log(initial)
    u,v,w = 0.234,0.342,-0.934
    rates = np.zeros(windows)
    cloud = np.empty((retain,3))
    minlog = 0.
    total = burn+windows*length
    for i in range(total):
        ex,ey,ez = math.exp(x),math.exp(y),math.exp(z)
        crowd = r*(1-a*ex-b*ey-c*ez)
        if survival >= 0:
            # In log coordinates the Jacobian row is a weighted average of
            # the two summands' logarithmic derivatives.
            live = survival*ex
            recruits = (1-survival)*math.exp(z+crowd)
            xn = math.log(live+recruits)
            fraction = recruits/(live+recruits)
            un = (1-fraction)*u+fraction*(-r*a*ex*u-r*b*ey*v+(1-r*c*ez)*w)
        else:
            xn = x+crowd
            un = (1-r*a*ex)*u-r*b*ey*v-r*c*ez*w
        norm = math.sqrt(un*un+u*u+v*v)
        if norm < 1e-280:
            return rates*math.nan,cloud,math.nan
        u,v,w = un/norm,u/norm,v/norm
        x,y,z = xn,x,y
        if i>=burn:
            rates[(i-burn)//length]+=math.log(norm)/length
            minlog=min(minlog,x)
        if i>=total-retain:
            cloud[i-total+retain]=math.exp(x),math.exp(y),math.exp(z)
    return rates,cloud,minlog


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def run(count,out,recruitment):
    out.mkdir(parents=True,exist_ok=True)
    seed = 2026112108 if recruitment else 2026112107
    rng = np.random.default_rng(seed)
    initial = np.array([.59,1.12,.47])
    passed = []
    with (out/"screen.csv").open("w",newline="",encoding="utf-8") as stream:
        csvout=csv.writer(stream)
        csvout.writerow(["case","r","a","b","c"]+(["survival"] if recruitment else [])+
                        ["min_log_abundance","min_window_rate"])
        for j in range(1,count+1):
            weights=rng.dirichlet([1.25,1.05,1.05])*.91 + np.array([.045,.02,.025])
            par=np.r_[rng.uniform(1.05,4.60),weights]
            if recruitment:
                par=np.r_[rng.uniform(2.2,5.1),weights,rng.uniform(.035,.68)]
            rates,z,low = orbit(par,initial,1100,4,450,400)
            csvout.writerow([j,*par,low,rates.min()])
            if rates.min()<.012 or not np.isfinite(rates).all() or low < -55:
                continue
            eig=np.linalg.eigvalsh(np.cov(z.T))
            if eig[0] < .003*eig[-1]:
                continue
            passed.append((j,par,float(rates.mean())))
    # A deterministic pilot subset represents the screen, then farthest-first
    # selection uses the archive fingerprint in native abundance coordinates.
    chosen = rng.choice(len(passed), min(1200,len(passed)),replace=False)
    records,descriptors=[],[]
    for idx in chosen:
        j,par,_=passed[idx]
        rates,z,low=orbit(par,initial,4800,4,1800,4000)
        if rates.min()<.01 or low< -55:
            continue
        if recruitment:
            scaled=(z-z.min(0))/(z.max()-z.min()+1e-20)
            cells=np.floor(scaled*22).astype(np.int32)
            _,counts=np.unique(cells[:,0]*529+cells[:,1]*23+cells[:,2],return_counts=True)
            if len(counts)<80 or counts.max()>.2*len(z):
                continue
        records.append({"case":j,"parameters":par.tolist(),"min_rate":float(rates.min()),
                        "mean_rate":float(rates.mean()),"minimum_log":low})
        descriptors.append(fingerprint(z))
    f=np.array(descriptors)
    if len(f):
        closest=np.full(len(f),np.inf)
        ids=[]
        k=int(np.argmax([r["min_rate"] for r in records]))
        for _ in range(min(60,len(f))):
            ids.append(k)
            ds=np.linalg.norm(f[:,:65]-f[k,:65],axis=1)+.5*np.linalg.norm(f[:,65:68]-f[k,65:68],axis=1)+.25*np.linalg.norm(f[:,68:]-f[k,68:],axis=1)
            closest=np.minimum(closest,ds)
            closest[ids]=-1
            k=int(np.argmax(closest))
    else:
        ids=[]
    selected=[records[k] for k in ids]
    for item in selected:
        rates,z,low=orbit(np.array(item["parameters"]),initial,12000,8,2500,16000)
        item["refined_min_rate"]=float(rates.min())
        if rates.min()>.009 and low> -55:
            np.savez_compressed(out/f"case-{item['case']}.npz",points=z)
    save(out/"screen.json",{"draws":count,"seed":seed,"recruitment_delay":recruitment,"early_passes":len(passed),
                            "replayed":len(chosen),"middle_passes":len(f),"diversity_sample":selected})
    print(count,len(passed),len(f),len(ids),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--count",type=int,default=16000)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--recruitment",action="store_true")
    args=ap.parse_args()
    run(args.count,args.out,args.recruitment)
