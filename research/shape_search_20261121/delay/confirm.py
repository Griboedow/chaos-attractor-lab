"""Replay selected delay maps in native abundance and logarithmic coordinates."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numba import njit

from search import orbit,save


@njit(cache=True)
def native_spectrum(par,initial,burn,windows,length,retain):
    r,a,b,c=par[:4]
    s=par[4] if len(par)==5 else -1.
    x,y,z=initial
    q=np.eye(3)
    powers=np.zeros((windows,3))
    determinants=np.zeros(windows)
    samples=np.zeros((retain,3))
    for time in range(burn+length*windows):
        e=math.exp(r*(1-a*x-b*y-c*z))
        if s>=0:
            recruit=(1-s)*z*e
            n=s*x+recruit
            j1,j2,j3=s-r*a*recruit,-r*b*recruit,(1-s)*e*(1-r*c*z)
        else:
            n=x*e
            j1,j2,j3=e*(1-r*a*x),-r*b*n,-r*c*n
        v=np.empty((3,3))
        for k in range(3):
            v[0,k]=j1*q[0,k]+j2*q[1,k]+j3*q[2,k]
            v[1,k]=q[0,k]
            v[2,k]=q[1,k]
        for k in range(3):
            for _ in range(2):
                for l in range(k):
                    dot=v[0,k]*v[0,l]+v[1,k]*v[1,l]+v[2,k]*v[2,l]
                    for h in range(3):
                        v[h,k]-=dot*v[h,l]
            norm=math.sqrt(v[0,k]**2+v[1,k]**2+v[2,k]**2)
            if norm < 1e-280 or abs(j3)<1e-280:
                raise ValueError("Native tangent underflow")
            for h in range(3):
                v[h,k]/=norm
            if time>=burn:
                powers[(time-burn)//length,k]+=math.log(norm)/length
        if time>=burn:
            determinants[(time-burn)//length]+=math.log(abs(j3))/length
        q=v
        x,y,z=n,x,y
        if time>=burn+length*windows-retain:
            samples[time-burn-length*windows+retain]=x,y,z
    return powers,determinants,samples


def confirm(row,out,recruitment,long):
    p=np.array(row["parameters"])
    initial=np.array([.59,1.12,.47])
    n=250000 if long else 12000
    windows=8 if long else 4
    burn=60000 if long else 6000
    retain=180000 if long else 1000
    rates,logs,cloud=native_spectrum(p,initial,burn,windows,n,retain)
    primary,_,minimum=orbit(p,initial,burn,windows,n,0)
    det_error=float(np.max(np.abs(rates.sum(1)-logs)))
    other=[]
    if long:
        rng=np.random.default_rng(2026112109)
        for z in rng.uniform(.12,2.4,size=(12,3)):
            rates1,_,low=orbit(p,z,12000,4,12000,0)
            other.append({"initial":z.tolist(),"window_rates":rates1.tolist()})
        np.savez_compressed(out/f"case-{row['case']}.npz",points=cloud)
    return {"case":row["case"],"parameters":p.tolist(),"initial":initial.tolist(),
            "family":"Delayed recruitment" if recruitment else "Delayed crowding",
            "window_length":n,"windows":windows,"discarded":burn,
            "native_spectra":rates.tolist(),"log_coordinate_rates":primary.tolist(),
            "native_mean_spectrum":rates.mean(0).tolist(),"maximum_log_determinant_error":det_error,
            "full_qr_determinant_check":det_error<1e-5,
            "minimum_log_abundance":minimum,"other_initials":other,
            "positive_extra_starts":sum(min(x["window_rates"])>0 for x in other),
            "cloud":f"case-{row['case']}.npz" if long else None}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("screen",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--cases",nargs="+",type=int)
    ap.add_argument("--long",action="store_true")
    args=ap.parse_args()
    data=json.loads(args.screen.read_text(encoding="utf-8"))
    rows=data if isinstance(data,list) else data["diversity_sample"]
    if args.cases:
        rows=[v for v in rows if v["case"] in args.cases]
        assert len(rows)==len(set(args.cases))
    args.out.mkdir(parents=True,exist_ok=True)
    reports=[]
    for row in rows:
        try:
            result=confirm(row,args.out,len(row["parameters"])==5,args.long)
        except ValueError as exc:
            result={"case":row["case"],"failure":str(exc)}
        reports.append(result)
        if args.long:
            save(args.out/"manifest.json",reports)
    save(args.out/"manifest.json",reports)
    print(f"{len(reports)} cases replayed",flush=True)
