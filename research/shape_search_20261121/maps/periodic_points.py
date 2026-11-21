"""Locate and rigorously isolate representative periodic orbits in the two slabs."""
import argparse
import itertools
import json
import sys
from pathlib import Path
from search import RESULTS,REPO
from cover_scout import charts
from validate import write_json
from verify_boxes import column,frames,rotate,values
from verify_return import tighten
from flint import arb,arb_mat,ctx
from numba import njit
import numpy as np
from scipy.optimize import brentq, root

sys.path.insert(0,str(REPO/"research"))
from twist_map import step_tangent


@njit(cache=True)
def iterate(x,p,n):
    v=np.eye(3)
    for _ in range(n):
        x,v=step_tangent(x,v,p)
    return x,v


def words(n):
    for w in itertools.product((0,1),repeat=n):
        if any(w==w[d:]+w[:d] for d in range(1,n)):
            continue
        if w==min(w[d:]+w[:d] for d in range(n)):
            yield w


def locate(w,observed,o,M,inv,p,cover):
    n=len(w);f=cover["faces"];axis=cover["axis"];k=cover["map_steps"]
    total=len(observed)-n-1
    score=np.linalg.norm(observed[:total]-observed[n:n+total],axis=1)
    for j,b in enumerate(w):
        u=observed[j:j+total,axis];lo,hi=f[2*b:2*b+2]
        score+=10*np.maximum(0,np.maximum(lo-u,u-hi))
    def both(flat):
        points=flat.reshape(n,3)
        res=np.empty_like(points);jac=np.zeros((3*n,3*n))
        for j,u in enumerate(points):
            nextj=(j+1)%n
            y,J=iterate(o+M@u,p,k)
            res[j]=inv@(y-o)-points[nextj]
            a=slice(3*j,3*j+3);b=slice(3*nextj,3*nextj+3)
            jac[a,a]=inv@J@M
            jac[a,b]-=np.eye(3)
        return res.ravel(),jac
    def guesses():
        for ix in np.argsort(score)[:24]:
            yield observed[ix:ix+n]
        # The certified face crossings bracket any target longitudinal value.
        # Inverting only this coordinate and propagating transverse coordinates
        # gives seeds even for rare words absent from the reference trajectory.
        trial=np.zeros((n,3))
        trial[:,axis]=[(f[2*b]+f[2*b+1])/2 for b in w]
        for sweep in range(36):
            for i,b in enumerate(w):
                j=(i+1)%n
                target=trial[j,axis]
                def residual(u):
                    point=trial[i].copy()
                    point[axis]=u
                    y,_=iterate(o+M@point,p,k)
                    return (inv@(y-o))[axis]-target
                lo,hi=f[2*b:2*b+2]
                trial[i,axis]=brentq(residual,lo,hi,xtol=1e-13)
                y,_=iterate(o+M@trial[i],p,k)
                previous=trial[j,axis]
                trial[j]=np.clip(inv@(y-o),-1,1)
                trial[j,axis]=previous
            if sweep in (1,5,12,24,35):
                yield trial.copy()
    for seed in guesses():
        sol=root(both,seed.ravel(),jac=True,tol=1e-11)
        candidate=sol.x.reshape(n,3)
        if np.max(np.abs(both(sol.x)[0]))>1e-9:
            continue
        if all(f[2*b]<u[axis]<f[2*b+1] and np.max(abs(u))<1 for u,b in zip(candidate,w)):
            return candidate
    return None


def certify_multiple(candidate,w,packet,cover):
    """Precondition the full cyclic system so a long single iterate cannot inflate the enclosure."""
    ctx.prec=192
    origins,bases,ivs=frames(packet)
    o,M,inv=origins[0],bases[0],ivs[0]
    p=[arb(float(v)) for v in packet["params"]]
    n=len(w)
    dim=3*n
    identity=arb_mat([[int(i==j) for j in range(dim)] for i in range(dim)])
    ident3=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    axis,k=cover["axis"],cover["map_steps"]
    faces=cover["faces"]
    for exponent in (34,38,42,30):
        rad=2.**-exponent
        full=arb_mat(dim,dim)
        residual=arb_mat(dim,1)
        monodromy=ident3
        valid=True
        for i,(u,b) in enumerate(zip(candidate,w)):
            local=[arb(v,rad) for v in u]
            if not (all(abs(v)<1 for v in local) and local[axis]>faces[2*b] and local[axis]<faces[2*b+1]):
                valid=False
                break
            center=o+M*column(u)
            wide=o+M*column(local)
            jac=M
            for _ in range(k):
                center,_=rotate(values(center),p)
                wide,jac=rotate(values(wide),p,jac)
                wide=tighten(wide,center,jac,[rad]*3)
            j=(i+1)%n
            r=inv*(center-o)-column(candidate[j])
            block=inv*jac
            monodromy=block*monodromy
            for h in range(3):
                residual[3*i+h,0]=r[h,0]
                for l in range(3):
                    full[3*i+h,3*i+l]+=block[h,l]
                full[3*i+h,3*j+h]-=1
        if not valid:
            continue
        midpoint=arb_mat([[full[i,j].mid() for j in range(dim)] for i in range(dim)])
        raw=midpoint.inv()
        pre=arb_mat([[raw[i,j].mid() for j in range(dim)] for i in range(dim)])
        if pre.det().contains(0):
            continue
        error=identity-pre*full
        slopes=[sum(error[i,j].abs_upper() for j in range(dim)) for i in range(dim)]
        offsets=values(pre*residual)
        if not all(s<1 and delta.abs_upper()+s*rad<rad for delta,s in zip(offsets,slopes)):
            continue
        trace=sum(monodromy[i,i] for i in range(3))
        c2=sum(monodromy[i,i]*monodromy[j,j]-monodromy[i,j]*monodromy[j,i] for i in range(3) for j in range(i))
        det=p[0]**(3*k*n)
        margin=trace.abs_lower()-(1+c2.abs_upper()+det.abs_upper())
        return dict(normalized_center=candidate[0].tolist(),normalized_itinerary_centers=candidate.tolist(),
            radius=rad,contraction_norm_upper=max(float(s) for s in slopes),method="multiple shooting",
            unique_root=True,word="".join(map(str,w)),exact_period=k*n,
            native_center_midpoint=[float(q.mid()) for q in values(o+M*column(candidate[0]))],
            trace_midpoint=float(trace.mid()),rouche_margin_lower_approximation=float(margin.lower()),
            certified_unstable_multipliers=(1 if margin>0 else None))
    return None


def certify(u,w,packet,cover):
    ctx.prec=192
    origins,bases,ivs=frames(packet);o,M,inv=origins[0],bases[0],ivs[0]
    p=[arb(float(v)) for v in packet["params"]]
    I=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    origin=column(u)
    axis,k=cover["axis"],cover["map_steps"];faces=cover["faces"]
    for exponent in (34,38,42,46,50):
        rad=2.**-exponent
        x=o+M*origin
        wide=o+M*column([arb(v,rad) for v in u])
        jac=M
        center_jac=M
        itinerary=True
        for t in range(k*len(w)):
            if t%k==0:
                normalized=values(tighten(inv*(wide-o),inv*(x-o),inv*jac,[rad]*3))
                b=w[t//k]
                if not (normalized[axis]>faces[2*b] and normalized[axis]<faces[2*b+1] and all(abs(v)<1 for v in normalized)):
                    itinerary=False;break
            x,center_jac=rotate(values(x),p,center_jac)
            wide,jac=rotate(values(wide),p,jac)
            wide=tighten(wide,x,jac,[rad]*3)
        if not itinerary:
            continue
        fc=inv*(x-o)-origin
        mid_jac=inv*center_jac-I
        precondition=arb_mat([[mid_jac.inv()[i,j].mid() for j in range(3)] for i in range(3)])
        if precondition.det().contains(0):
            continue
        monodromy=inv*jac
        error=I-precondition*(monodromy-I)
        slopes=[sum(error[i,j].abs_upper() for j in range(3)) for i in range(3)]
        offsets=values(precondition*fc)
        if not all(s<1 and d.abs_upper()+s*rad<rad for d,s in zip(offsets,slopes)):
            continue
        trace=sum(monodromy[i,i] for i in range(3))
        c2=sum(monodromy[i,i]*monodromy[j,j]-monodromy[i,j]*monodromy[j,i] for i in range(3) for j in range(i))
        det=p[0]**(3*k*len(w))
        margin=trace.abs_lower()-(1+c2.abs_upper()+det.abs_upper())
        return dict(normalized_center=u.tolist(),radius=rad,contraction_norm_upper=max(float(s) for s in slopes),
            unique_root=True,word="".join(map(str,w)),exact_period=k*len(w),native_center_midpoint=[float(q.mid()) for q in values(o+M*origin)],
            trace_midpoint=float(trace.mid()),rouche_margin_lower_approximation=float(margin.lower()),
            certified_unstable_multipliers=(1 if margin>0 else None))
    return None


def build(results,maxword):
    cert=results/"certificates"
    box=json.loads((cert/"2503155-symbolic-box.json").read_text())
    cover=json.loads((cert/"2503155-faces.json").read_text())[0]
    row=json.loads((results/"checks/2503155.json").read_text())
    o,M,inv=charts(box);p=np.array(row["params"])
    z=np.load(row["cloud"])["points"].astype(float)
    normalized=(z-o)@inv.T
    stride=cover["map_steps"]
    phase=max(range(stride),key=lambda j:(np.abs(normalized[j::stride]).max(axis=1)<1).mean())
    observed=normalized[phase::stride]
    successes=[];unresolved=[]
    for n in range(1,maxword+1):
        for w in words(n):
            found=locate(w,observed,o,M,inv,p,cover)
            checked=None if found is None else certify_multiple(found,w,box,cover)
            if not checked and found is not None:
                checked=certify(found[0],w,box,cover)
            if checked:
                successes.append(checked)
            else:
                unresolved.append("".join(map(str,w)))
    write_json(cert/"2503155-periodic-points.json",dict(case=2503155,arb_precision_bits=192,maximum_word_length=maxword,
        interval_certified_orbits=successes,numerically_unresolved_representatives=unresolved))
    print("Certified cycles",len(successes),"with one unstable multiplier",sum(x["certified_unstable_multipliers"]==1 for x in successes),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--results",type=Path,default=RESULTS)
    ap.add_argument("--maxword",type=int,default=5)
    a=ap.parse_args();build(a.results,a.maxword)
