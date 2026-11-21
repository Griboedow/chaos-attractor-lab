"""Certify attracting periodic responses at finite gains by an a posteriori moving-disk bound."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from flint import arb,acb,ctx
from certify import real_arb,terms
from gradient import independent
from model import invariants

ctx.prec=104


def dyadic_word(curves,direction=(5,1)):
    ratios=[[(float(v.real).as_integer_ratio(),float(v.imag).as_integer_ratio()) for v in row] for row in curves]
    denominator=max(d for row in ratios for point in row for _,d in point)
    a,b=direction
    def frame(row):
        out=[]
        for (x,dx),(y,dy) in row:
            x*=denominator//dx;y*=denominator//dy
            out.append((a*x+b*y,a*y-b*x))
        return out
    before=frame(ratios[0]);ordering=sorted(range(3),key=lambda i:before[i][0]);word=[]
    for row in ratios[1:]:
        after=frame(row);order=sorted(range(3),key=lambda i:after[i][0])
        assert len({a for a,b in after})==3
        ix=[i for i in range(3) if order[i]!=ordering[i]]
        if ix:
            assert len(ix)==2 and ix[1]==ix[0]+1
            j=ix[0];left,right=ordering[j:j+2]
            x0=before[left][0]-before[right][0];x1=after[left][0]-after[right][0]
            y0=before[left][1]-before[right][1];y1=after[left][1]-after[right][1]
            height=y0*x1-y1*x0
            assert x1>x0 and height!=0
            word.append((j+1)*(1 if height>0 else -1))
        before,ordering=after,order
    return word


def bounds(row,gain,curves,cap,start):
    from fractions import Fraction
    kappa=arb(gain.numerator)/gain.denominator
    n=len(curves)-1;dt=2*arb.pi()/n
    B=arb(row["drive_speed_bound"])
    radii=[arb(start)]*3;worst=0.;total=[arb(0)]*3;mingap=1.e99
    # All floating input vertices and the radii are exact binary numbers in Arb.
    for j,(u,v) in enumerate(zip(curves[:-1],curves[1:])):
        theta=arb.pi()*(2*j+1)/n
        b=acb(real_arb(theta,row),row["imag_amplitude"]*(row["k"]*theta).sin())
        mid=[];half=[]
        for i in range(3):
            z0=acb(float(u[i].real),float(u[i].imag));z1=acb(float(v[i].real),float(v[i].imag))
            m=(z1+z0)/2;l=abs(z1-z0)/2
            mid.append(m);half.append(l)
            J=3*m*m-3;P=m**3-3*m+b;am=abs(m);aj=abs(J);ap=abs(P)
            def remainder(rad):
                return (aj*rad+3*am*rad**2+rad**3+B*dt/2,6*am*rad+3*rad**2)
            dp,dj=remainder(l)
            error=abs((z1-z0)/dt+kappa*J.conjugate()*P)+kappa*((aj+dj)*dp+dj*ap)
            dp,dj=remainder(l+cap)
            lower=max(arb(0),(aj-dj).lower())
            mu=(kappa*(-lower**2+6*(am+l+cap)*(ap+dp))).upper()
            total[i]+=mu*dt
            growth=(mu*dt).exp()
            end=(growth*radii[i]+max(arb(1),growth)*dt*error).upper()
            if not max(end,radii[i])<cap:
                return dict(success=False,failed_at=j,component=i,radius_upper=float(end),cap=cap)
            radii[i]=end;worst=max(worst,float(end))
        for i,k in ((0,1),(0,2),(1,2)):
            gap=(abs(mid[i]-mid[k])-half[i]-half[k]-2*cap).lower()
            assert gap>0
            mingap=min(mingap,float(gap))
    success=all(x<start for x in radii) and all(x<0 for x in total)
    return dict(success=success,start_radius=start,cap=cap,max_propagated_radius=worst,
                end_radius_upper=[float(x) for x in radii],log_norm_integral_upper=[float(x) for x in total],
                minimum_gap_after_disk_radii=mingap)


if __name__=="__main__":
    from fractions import Fraction
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates",type=Path);ap.add_argument("survey",type=Path)
    ap.add_argument("--out",type=Path,required=True);ap.add_argument("--cache",type=Path,required=True)
    ap.add_argument("--case",type=int,default=48000001);ap.add_argument("--cells",type=int,default=65536)
    args=ap.parse_args()
    row=next(r for r in json.loads(args.certificates.read_text())["rows"] if r["case"]==args.case)
    survey=next(s for s in json.loads(args.survey.read_text()) if s["case"]==args.case)
    records=[]
    for gain in (Fraction(3,100),Fraction(34,100)):
        s=next(s for s in survey["runs"] if s["gain"]==float(gain))
        curves=np.column_stack([independent(row,float(gain),u,samples=args.cells+1)[0] for u in s["starts_last_period"]])
        curves[-1]=curves[0]
        checks=[]
        choices=[(.015625,.0088)] if gain==Fraction(3,100) else [(.0078125,.0055)]
        for cap,start in choices:
            assert row["embedding_radius"]+curves.real.min()-cap>0
            check=bounds(row,gain,curves,cap,start)
            checks.append(check)
            if check["success"]:break
        word=dyadic_word(curves);lk,mu=invariants(word)
        import hashlib
        args.cache.mkdir(parents=True,exist_ok=True)
        path=args.cache/f"centers-{args.case}-{gain.numerator}-{gain.denominator}.npz"
        np.savez_compressed(path,centers=curves)
        records.append(dict(gain=[gain.numerator,gain.denominator],precision_bits=ctx.prec,cells=args.cells,
                            checks=checks,absolute_mu=abs(mu),pairwise_linking=lk,dyadic_crossing_word=word,
                            center_file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        args.out.write_text(json.dumps(dict(case=row["case"],records=records),indent=2)+"\n")
        print(gain,checks,flush=True)
        if not checks[-1]["success"]:
            raise RuntimeError(f"Moving disks failed to certify the three cycles at gain {gain}")
