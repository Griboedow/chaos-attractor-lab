"""Outward-rounded Arb/Picard-Taylor enclosure of individual points under the C7 chart."""
import argparse
import json
import math
import sys
from pathlib import Path
from time import monotonic
from flint import arb,arb_series,ctx
from native_section import P,RAW,RESULT,packed
from return_tube import graph
from reduced import reduced
import numpy as np


def binary(x):
    n,d=float(x).as_integer_ratio()
    return arb(n)/arb(d)


def field(a,p):
    x,y,z=a
    r=x*x+y*y
    b=2/(1+r)
    scale=b*b*b.sqrt()
    u=x*(x*x-3*y*y);v=y*(3*x*x-y*y)
    c=u*u-v*v;s=2*u*v
    alpha=p[1]-r+p[2]*z;w=p[3]+p[4]*z+p[5]*r
    return [x*alpha-y*w+p[6]*scale*c,
            y*alpha+x*w-p[6]*scale*s,
            -p[7]*z+p[8]*(r-1)/(1+r)+p[9]*scale*b*(p[11]*(x*c-y*s)-p[12]*(y*c+x*s))]


def jets(a,p,n):
    coeff=[[q] for q in a]
    for k in range(n):
        t=[arb_series(v,k+1) for v in coeff]
        f=field(t,p)
        for c,j in zip(coeff,f):c.append(j[k]/(k+1))
    return coeff


def advance(a,h,p,n):
    f=field(a,p);sw=arb(h/2,h/2)
    pads=[2*h*w.abs_upper()+arb('1e-38') for w in f]
    for _ in range(16):
        tube=[arb(v.mid(),v.rad()+s) for v,s in zip(a,pads)]
        ft=field(tube,p)
        b=[v+sw*w for v,w in zip(a,ft)]
        if all(q.contains_interior(v) for q,v in zip(tube,b)):break
        pads=[max(arb(5)/4*s,2*(v-w).abs_upper()) for s,v,w in zip(pads,b,a)]
    else:raise ArithmeticError('Failed Picard inclusion')
    jac=[]
    for k in range(3):
        row=field([arb_series([q,int(i==k)],2) for i,q in enumerate(tube)],p)
        jac.append([v[1].abs_upper() for v in row])
    lipschitz=h*max(sum(u) for u in zip(*jac))
    if not lipschitz<1:raise ArithmeticError('Failed contraction')
    cn=jets(a,p,n);rem=jets(tube,p,n+1)
    b=[]
    for c,r in zip(cn,rem):
        v=c[-1]
        for d in reversed(c[:-1]):v=v*h+d
        b.append(v+r[-1]*(h**(n+1)))
    return b,tube,ft,float(lipschitz.upper())


def record(v):
    m,e=v.mid().man_exp();s,f=v.rad().man_exp()
    return dict(display=v.str(27),midpoint=[str(m),int(e)],radius=[str(s),int(f)])


def run(val,delta,rad,steps,n,precision):
    ctx.prec=precision;ctx.cap=n+2
    p=[binary(q) for q in P];p += [p[10].cos(),p[10].sin()]
    z=arb(val);g=[binary(x) for x in graph().coef];v=z-binary(.99);x=g[-1]
    for q in reversed(g[:-1]):x=x*v+q
    x=x+arb(delta)
    trial=reduced(np.array([float(x)]),np.array([float(z)]),.00625,packed())[0]
    guess=round(trial[0],8)
    center=arb(str(guess));width=arb('0.00000003')
    left,right=center-width,center+width
    phi=8*arb.pi()/7
    x=arb(x,arb(rad));z=arb(z,arb(rad))
    a=[x*phi.cos(),x*phi.sin(),z];maxrad=0.;max_lip=0.
    gate_clear=True;h=left/steps
    for k in range(steps):
        try:a,tube,ft,lip=advance(a,h,p,n)
        except ArithmeticError as exc:raise ArithmeticError((k,str(exc),[v.str(7) for v in a])) from exc
        maxrad=max(maxrad,max(float(q.rad()) for q in a));max_lip=max(max_lip,lip)
        if (k+1)*h>=arb('5.7'):
            if not (tube[1]<0 or tube[0]<arb('.77') or tube[0]>arb('.85') or
                    tube[2]<arb('.945') or tube[2]>arb('1.04') or (ft[1]>0 and a[1]<0)):
                gate_clear=False
    start=a
    a,tube,ft,lip=advance(a,right-left,p,n)
    passed=start[1]<0 and a[1]>0 and ft[1]>0 and gate_clear
    # Intersect x,z from a Picard enclosure with the chart limits.
    zret=tube[2];xret=tube[0];v=zret-binary(.99);j=g[-1]
    for q in reversed(g[:-1]):j=j*v+q
    delta_ret=xret-j
    inside=zret>arb('.98071') and zret<arb('.99943') and delta_ret.abs_upper()<arb('1e-6')
    return dict(status='certified_return' if passed else 'inconclusive',reduced_tube_membership_certified=bool(inside and passed),
            initial_z=val,initial_x_offset=delta,independent_native_initial_radius=rad,
            precision_bits=precision,steps=steps,taylor_order=n,max_step_state_radius=maxrad,max_picard_contraction=max(max_lip,lip),
            event_slab=[record(left),record(right)],negative_before=record(start[1]),positive_after=record(a[1]),
            positive_normal=record(ft[1]),event_x=record(xret),event_z=record(zret),event_delta=record(delta_ret),
            preceding_gate_excluded=bool(gate_clear),graph_binary64=[float(q) for q in graph().coef],
            parameter_binary64=P.tolist())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--z',default='.99');ap.add_argument('--delta',default='0');ap.add_argument('--radius',default='1e-44')
    ap.add_argument('--steps',type=int,default=560);ap.add_argument('--order',type=int,default=27)
    ap.add_argument('--bits',type=int,default=224);ap.add_argument('--id',default='center')
    a=ap.parse_args();start=monotonic()
    try:r=run(a.z,a.delta,a.radius,a.steps,a.order,a.bits)
    except Exception as exc:r=dict(status='inconclusive',error=repr(exc))
    r['seconds']=monotonic()-start
    print(r.get('status'),r.get('reduced_tube_membership_certified'),r.get('error'),r['seconds'],flush=True)
    (RAW/f'arb-{a.id}.json').write_text(json.dumps(r,indent=2)+'\n')
