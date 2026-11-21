"""Propose and verify cyclic ellipsoids using 192-bit Arb enclosures."""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/"research/shape_search_20261121/maps"))
import numpy as np
from flint import arb,arb_mat,ctx
from scipy.spatial.distance import cdist
from periodic_points import iterate
from search import sample
from verify_boxes import column,matrix,rotate,values


def isolate(p,start,requested=None):
    w,z=sample(p[None,:],np.array(start)[None,:],72000,4096,4,1024)
    t=z[0];ret=np.flatnonzero(np.linalg.norm(t[-1]-t[-513:-1][::-1],axis=1)<1e-8)+1
    if not len(ret):return
    m=int(ret[0])
    if requested is not None and m!=requested:return
    if np.max(np.abs(t[-2*m:-m]-t[-m:]))>1e-8:return
    xs=t[-m:].copy()
    js=np.array([iterate(x,p,1)[1] for x in xs])
    # Iterating the covariance equations is only a source of trial charts.
    # Arb tests all conclusions against the original map independently.
    P=np.eye(3)
    for cycle in range(300):
        prev=P.copy()
        for J in js:
            P=J@P@J.T+np.eye(3)
        if not np.isfinite(P).all() or abs(P).max()>1e18:return
        if np.linalg.norm(P-prev)<=1e-12*np.linalg.norm(P):break
    else:return
    frames=[]
    for J in js:
        frames.append(np.linalg.cholesky((P+P.T)/2))
        P=J@P@J.T+np.eye(3)
    f=np.array(frames)
    norms=[np.linalg.norm(np.linalg.solve(f[(i+1)%m],js[i]@f[i]),2) for i in range(m)]
    if max(norms)>=1:return
    factors=((1+np.array(norms))/2).tolist()
    for exp in (33,36,39,30,42,27,45):
        cert=dict(period=m,centers=xs.tolist(),frames=f.tolist(),radius_bits=exp,factors=factors,
            initial_state=list(start),numerical_leading_rate=float(w.mean()))
        result=verify(p,cert)
        if result:
            cert["verification"]=result
            return cert
    # Recenter in high precision when weak contraction or highly anisotropic
    # charts make the roundoff of binary64 orbit proposals too large.
    ctx.prec=192
    pars=[arb(float(v)) for v in p]
    x=[arb(float(v)) for v in xs[0]]
    for turn in range(260):
        x0=x[:];phases=[]
        for i in range(m):
            phases.append([[int(a) for a in v.man_exp()] for v in x])
            y,_=rotate(x,pars);x=[v.mid() for v in values(y)]
        if all(abs(v-w)<arb(2.)**-116 for v,w in zip(x,x0)):break
    else:return
    for exp in (50,60,72,84,96):
        cert=dict(period=m,centers_dyadic=phases,frames=f.tolist(),radius_bits=exp,factors=factors,
            initial_state=list(start),numerical_leading_rate=float(w.mean()))
        result=verify(p,cert)
        if result:
            cert["verification"]=result
            return cert


def verify(p,cert):
    ctx.prec=192
    p=[arb(float(v)) for v in p]
    m=cert["period"];rad=arb(2.)**(-cert["radius_bits"])
    if "centers_dyadic" in cert:
        X=[column([arb(tuple(v)) for v in x]) for x in cert["centers_dyadic"]]
    else:
        X=[column(x) for x in cert["centers"]]
    M=[matrix(f) for f in cert["frames"]]
    if not (len(X)==len(M)==len(cert["factors"])==m and m>0):raise ValueError("Missing phases")
    native=[];worst=0.;slack=1.
    for i,(x,S) in enumerate(zip(X,M)):
        j=(i+1)%m
        q=arb(float(cert["factors"][i]))
        if not (q>0 and q<1):return
        V=M[j].inv()
        enclosure=x+S*column([arb(0,rad) for _ in range(3)])
        y,_=rotate(values(x),p)
        _,K=rotate(values(enclosure),p,S)
        K=V*K
        symmetric=arb_mat([[q*q*int(a==b)-sum(K[k,a]*K[k,b] for k in range(3))
            for b in range(3)] for a in range(3)])
        minors=[symmetric[0,0],symmetric[0,0]*symmetric[1,1]-symmetric[0,1]*symmetric[1,0],symmetric.det()]
        shift=sum(v.abs_upper()**2 for v in values(V*(y-X[j]))).sqrt()
        if not (all(d>0 for d in minors) and shift+q*rad<rad):return
        worst=max(worst,float(q));slack=min(slack,float(((1-q)-shift/rad).lower()))
        native.append(values(enclosure))
    # One native coordinate separates each pair. Every actual ball lies in
    # the associated cubical over-enclosure used here.
    for i in range(m):
        for j in range(i):
            if not any(native[i][k]>native[j][k] or native[i][k]<native[j][k] for k in range(3)):return
    return dict(precision_bits=192,all_phase_pairs_disjoint=True,maximum_exact_binary64_factor=worst,
        minimum_relative_trapping_slack_approximation=slack)


def run(cache,out):
    hits=json.loads((cache/"scout.json").read_text())
    out.mkdir(parents=True,exist_ok=True)
    # Only one parameter per distinct pair of certified periods.
    seen=set();done=[]
    for m in sorted({h["period"] for h in hits},reverse=True):
        for h in sorted([h for h in hits if h["period"]==m],key=lambda h:h["rate"])[:6]:
            p=np.array(h["params"]);c=isolate(p,h["start"],m)
            if not c:continue
            certified=[c]
            for start in h["alternate_starts"]:
                d=isolate(p,start)
                if d and d["period"] not in [c["period"] for c in certified]:
                    certified.append(d)
            orders=tuple(sorted(x["period"] for x in certified))
            if orders in seen:continue
            seen.add(orders)
            packet=dict(id=h["id"],parent=h["parent"],index=h["index"],params=h["params"],cycles=certified)
            f=cache/("certificate-"+h["id"]+".json")
            f.write_text(json.dumps(packet,separators=(",",":"))+"\n")
            done.append(dict(id=h["id"],parent=h["parent"],periods=list(orders),certificate=f.name))
            (out/"cycles.json").write_text(json.dumps(done,indent=2)+"\n")
            print(h["id"],orders,flush=True)
            break


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("cache",type=Path)
    ap.add_argument("output",type=Path);a=ap.parse_args();run(a.cache,a.output)
