"""Prove a critically periodic parameter and fifteen basin witnesses near one pole with Arb balls."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from flint import acb,arb,ctx,__version__ as flint_version
from selected import D,K,CASE,APPROX
from search import OUTPUT,RAW,roots

def bb(v):
    n,e=v.mid().man_exp()
    return [str(n),int(e)]
def exact_xy(v):return [bb(v.real),bb(v.imag)]
def from_xy(z):return acb(arb((int(z[0][0]),int(z[0][1]))),arb((int(z[1][0]),int(z[1][1]))))
def box(v,rad):return acb(arb(v.real.mid(),rad),arb(v.imag.mid(),rad))
def cp(v):return complex(float(v.real.mid()),float(v.imag.mid()))
def nmap(z,a):
    if abs(z)>3:
        c=arb(D-1)/D
        return c*z+a*(c*z-1)/(D*z**(D-1)-a)
    return ((D-1)*z**D-a)/(D*z**(D-1)-a)
def dzmap(z,a):
    p=z**D-a*z+a
    return D*(D-1)*z**(D-2)*p/(D*z**(D-1)-a)**2
def parameter_eq(a):
    z=acb(0);der=acb(0);minpole=None
    for _ in range(K):
        B=D*z**(D-1)-a;A=(D-1)*z**D-a
        gap=abs(B).lower();assert gap>0
        minpole=gap if minpole is None else min(minpole,gap)
        der=dzmap(z,a)*der+(A-B)/B**2;z=A/B
    return z,der,minpole
def refine(seed,func):
    c=acb(seed).mid()
    for _ in range(15):
        value,der=func(c)[:2];c=(c-value/der).mid()
    return c
def kcheck(c,rho,fun):
    f,d=fun(c)[:2];Y=(1/d).mid();u=box(c,rho)
    fbox,db=fun(u)[:2]
    q=abs(1-Y*db).upper();shift=abs(Y*f).upper()
    assert q<1 and shift+q*rho<rho,(str(q),str(shift),str(rho))
    return dict(center_dyadic=exact_xy(c),radius=str(rho),q_upper=float(q),shift_upper=float(shift))
def configure():
    rho=arb('1e-142')
    c=refine(APPROX,parameter_eq);A=box(c,rho)
    par=kcheck(c,rho,parameter_eq)
    cyc=[acb(0)]
    for i in range(1,K):cyc.append(nmap(cyc[-1],A))
    separation=min(abs(u-v).lower() for j,u in enumerate(cyc) for v in cyc[:j]);assert separation>0
    # Every return denominator has a strictly positive lower bound on the full parameter box.
    _,_,pole_margin=parameter_eq(A)
    root_data=[];allroots=[]
    for seed in roots(D,cp(c)):
        r=refine(seed,lambda z:(z**D-c*z+c,D*z**(D-1)-c))
        rr=arb('1e-112')
        rc=kcheck(r,rr,lambda z:(z**D-A*z+A,D*z**(D-1)-A))
        rb=box(r,rr);allroots.append(rb)
        local=abs(dzmap(box(r,arb('2e-7')),A)).upper()
        assert local<arb('.5')
        rc['capture_radius']='1e-7';rc['derivative_upper']=float(local);root_data.append(rc)
    rs=min(abs(u-v).lower() for j,u in enumerate(allroots) for v in allroots[:j]);assert rs>0
    z=box(acb(0),arb('1e-7'));ds=arb(1)
    for _ in range(K):
        ds*=abs(dzmap(z,A)).upper();z=nmap(z,A)
    assert ds<arb('.5')
    return c,A,cyc,allroots,dict(parameter_certificate=par,critical_orbit_expanded=[[float(t.real.mid()),float(t.imag.mid())] for t in cyc],
         critical_orbit_minimum_separation=float(separation),critical_denominator_minimum=float(pole_margin),
         pairwise_root_minimum=float(rs),root_certificates=root_data,cycle_zero_capture_radius='1e-7',
         cycle_return_derivative_upper=float(ds.upper()))
def one_inverse(target,c,stage):
    w=target
    if abs(w)>3:
        seed=acb(arb(D)/arb(D-1)*w)
    else:
        ww=cp(w);aa=cp(c);pol=np.zeros(D+1,dtype=complex);pol[0]=D-1;pol[1]=-D*ww;pol[-1]=aa*(ww-1)
        possible=np.roots(pol);possible=possible[abs(possible-ww)>1e-5]
        # Deterministic angular tie-break at the origin.
        seed=acb(max(possible,key=lambda z:(round(abs(z),8),z.imag,z.real)))
    x=seed.mid()
    for j in range(7):
        z=x
        F=(D-1)*z**D-D*w*z**(D-1)+c*(w-1)
        grad=D*(D-1)*z**(D-2)*(z-w)
        x=(z-F/grad).mid()
    return x
def near_pole(w,c,pole):
    v=1/w;x=pole.mid()
    for _ in range(6):
        F=(D-1)*x**D*v-D*x**(D-1)+c*(1-v)
        der=D*(D-1)*x**(D-2)*(x*v-1)
        x=(x-F/der).mid()
    return x
def verify_witness(record,A,allroots,p0):
    original=from_xy(record['center_dyadic']);u=box(original,arb('1e-101'))
    dist=abs(u-p0).upper();assert dist<arb('3e-32')
    z=u
    for i in range(record['outer_chain_length']+24):
        for j,r in enumerate(allroots):
            if abs(z-r)<arb('1e-7'):
                label=j;break
        else:
            if abs(z)<arb('1e-7'):label=D+((-i)%K)
            else:
                den=D*z**(D-1)-A
                assert abs(den).lower()>0,(i,record['outer_chain_length'])
                z=nmap(z,A);continue
        return int(label),int(i),float(dist),float(abs(z).upper())
    raise RuntimeError('No validated capture')
def make_witnesses(c,A,cyc,allroots):
    pole=(A/D).root(D-1)
    candidates=[]
    for j,target in enumerate([*allroots,acb(0)]):
        w=target.mid();stage=0
        while abs(w)<arb('2e30'):
            w=one_inverse(w,c,stage);stage+=1
            if stage>700:raise RuntimeError('Unable to reach large inverse branch')
        for extra in range(K if j==D else 1):
            w=one_inverse(w,c,stage);stage+=1
            x=near_pole(w,c,pole)
            q=dict(center_dyadic=exact_xy(x),outer_chain_length=stage,
                   target='critical cycle' if j==D else 'root '+str(j))
            label,i,dist,_=verify_witness(q,A,allroots,pole)
            q.update(certified_outcome=label,capture_iteration=i,upper_native_distance_from_pole=dist)
            candidates.append(q);print('witness',j,extra,label,i,dist,flush=True)
    assert {q['certified_outcome'] for q in candidates}==set(range(D+K))
    return dict(selected_pole='principal (a/6)^(1/5)',approximate_pole=[float(pole.real.mid()),float(pole.imag.mid())],
                common_disk_radius='3e-32',witness_start_box_native_coordinate_halfwidth='1e-101',witnesses=candidates)
def main():
    ctx.prec=656
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--verify',type=Path)
    args=ap.parse_args()
    c,A,cyc,allroots,report=configure()
    if args.verify:
        frozen=json.loads(args.verify.read_text());p0=(A/D).root(D-1)
        assert exact_xy(c)==frozen['parameter_certificate']['center_dyadic']
        labels=[]
        for item in frozen['witnesses']:
            label,step,_,_=verify_witness(item,A,allroots,p0)
            assert (label,step)==(item['certified_outcome'],item['capture_iteration']);labels.append(label)
        assert set(labels)==set(range(D+K))
        print('Verified exact parameter, roots, cycle, capture sets and all',len(labels),'witness boxes');return
    report.update(make_witnesses(c,A,cyc,allroots));report.update(case=CASE,d=D,k=K,a_approximate=[float(c.real.mid()),float(c.imag.mid())],
         python_flint=flint_version,precision_bits=ctx.prec,scope='Statements refer to the unique complex root inside the parameter disk, not to its rounded decimal approximation.')
    (OUTPUT/'certificate.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
