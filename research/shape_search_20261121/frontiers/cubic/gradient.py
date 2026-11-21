"""Track actual periodic responses of the single unswitched forced cubic gradient equation."""
import argparse
import json
from pathlib import Path

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp

from model import all_roots,forcing,geometric_word,invariants


@njit(cache=True)
def field(z,b,gain):
    return -gain*(3*z*z-3).conjugate()*(z*z*z-3*z+b)


@njit(cache=True)
def integrate(starts,b,gain,tours):
    ns=(len(b)-1)//2;dt=2*np.pi/ns
    paths=np.empty((ns+1,len(starts)),np.complex128)
    ends=np.empty((tours+1,len(starts)),np.complex128)
    for j in range(len(starts)):
        z=starts[j];ends[0,j]=z
        for turn in range(tours):
            if turn==tours-1:paths[0,j]=z
            for i in range(ns):
                a=field(z,b[2*i],gain)
                v=field(z+.5*dt*a,b[2*i+1],gain)
                w=field(z+.5*dt*v,b[2*i+1],gain)
                y=field(z+dt*w,b[2*i+2],gain)
                z+=dt/6*(a+2*v+2*w+y)
                if turn==tours-1:paths[i+1,j]=z
            ends[turn+1,j]=z
    return paths,ends


def record(paths,ends,gain):
    periods=[]
    for j in range(3):
        periods.append(next((k for k in range(1,13) if np.max(abs(ends[-5:,j]-ends[-5-k:-k,j]))<3e-7),None))
    gap=min(float(np.min(abs(paths[:,i]-paths[:,j]))) for i,j in ((0,1),(0,2),(1,2)))
    row=dict(gain=gain,least_observed_periods=periods,minimum_same_phase_pair_distance=gap,
             closure=float(np.max(abs(paths[0]-paths[-1]))),
             starts_last_period=np.column_stack((paths[0].real,paths[0].imag)).tolist())
    if periods==[1,1,1] and gap>1.e-4:
        views=[]
        for angle in (.129,.593):
            word=geometric_word(paths,angle)
            lk,mu=invariants(word);views.append(dict(angle=angle,pair_linking=lk,mu=mu))
        row["views"]=views
    return row


def sweep(row,levels,tours,ns):
    b=forcing(np.linspace(0,2*np.pi,2*ns+1),row)
    roots=all_roots(np.array([b[0]]))[0]
    result=[]
    for gain in levels:
        paths,ends=integrate(roots,b,float(gain),tours)
        result.append(record(paths,ends,float(gain)))
    return result


def independent(row,gain,initial,full_matrix=False,samples=4097):
    h=row["h"];c=np.asarray(row["coeffs"])/row["divisor"]
    cs=[(n,c[n],c[h+n]) for n in range(1,h+1) if c[n] or c[h+n]]
    def rhs(t,q):
        z=q[0]+1j*q[1]
        b=c[0]+sum(a*np.cos(n*t)+s*np.sin(n*t) for n,a,s in cs)+1j*row["imag_amplitude"]*np.sin(row["k"]*t)
        pd=3*z*z-3;res=z*z*z-3*z+b
        v=-gain*np.conj(pd)*res
        lin=-gain*abs(pd)**2
        anti=-gain*6*np.conj(z)*res
        J=np.array([[lin+anti.real,anti.imag],[anti.imag,lin-anti.real]])
        return np.r_[v.real,v.imag,(J@q[2:].reshape(2,2)).ravel()]
    result=solve_ivp(rhs,[0,2*np.pi],np.r_[initial,np.eye(2).ravel()],rtol=4e-11,atol=2e-13,method="DOP853",
                     t_eval=np.linspace(0,2*np.pi,samples))
    assert result.success
    q=result.y
    matrix=q[2:,-1].reshape(2,2)
    return q[0]+1j*q[1],matrix if full_matrix else max(abs(np.linalg.eigvals(matrix)))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("certificates",type=Path)
    ap.add_argument("--out",type=Path,required=True);ap.add_argument("--adaptive",action="store_true")
    ap.add_argument("--cases",type=int,nargs="+",default=[47087142,48000001,47053641]);a=ap.parse_args()
    certs=json.loads(a.certificates.read_text())["rows"];result=[]
    for row in (r for r in certs if r["case"] in a.cases):
        levels=[.01,.03,.06,.1,.15,.23,.34,.5,.75,1.,1.5,2.,3.,4.]
        runs=sweep(row,levels,42,8192)
        if a.adaptive:
            examples=[]
            for condition in (lambda s:all(v==1 for v in s["least_observed_periods"]) and s.get("views") and s["views"][0]["mu"]==0,
                              lambda s:s.get("views") and abs(s["views"][0]["mu"] or 0)>0):
                chosen=next((s for s in runs if condition(s)),None)
                if chosen:
                    pairs=[independent(row,chosen["gain"],v) for v in chosen["starts_last_period"]]
                    paths=np.column_stack([z for z,m in pairs]);lk,mu=invariants(geometric_word(paths,.413))
                    examples.append(dict(gain=chosen["gain"],linking=lk,mu=mu,
                                        dominant_multiplier_moduli=[m for z,m in pairs],
                                        closure=float(abs(paths[0]-paths[-1]).max())))
        else:examples=[]
        result.append(dict(case=row["case"],runs=runs,adaptive=examples))
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+"\n")
        print(row["case"],[(q["gain"],q.get("views",[{}])[0].get("mu"),q["least_observed_periods"]) for q in runs],flush=True)
