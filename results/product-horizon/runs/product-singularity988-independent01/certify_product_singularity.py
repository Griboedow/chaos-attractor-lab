"""Isolate a zero of the internal reflected vector using interval Krawczyk.

This is a singularity of the normalization, not an equilibrium of F. The
inverse-function theorem gives a spherical cluster set of nearby images.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
from product_binary_interval import parameters


def main(args):
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    import verify_product_coverings as full
    full.mp=mp;mp.mp.dps=mp.iv.dps=args.digits;iv=mp.iv
    exact=lambda x,ctx:ctx.mpf(Fraction(x).numerator)/Fraction(x).denominator
    enc=lambda v:dict(decimal=str(v),binary=[list(t) for t in v._mpi_])
    args.out.mkdir(parents=True,exist_ok=False)
    raw=json.loads(args.input.read_text(encoding='utf-8'))
    beta,tau,c,a=parameters(raw['params'],np.array(raw['points']),mp.mp)
    roots=[mp.findroot(lambda x:2*(c[j]+a[j]*mp.tanh(2*a[j]*(x-c[j])/tau))-x,c[j],solver='newton',tol=mp.mpf(10)**(-args.digits+20)) for j in range(3)]
    center=[mp.nstr(v,args.digits+10) for v in roots];r=exact(args.radius,iv)
    xi=[exact(v,iv) for v in center];box=[v+iv.mpf([-r.b,r.b]) for v in xi]
    bi,ti,ci,ai=parameters(raw['params'],np.array(raw['points']),iv)
    def binary(x):
        values=[];jac=[[iv.mpf(0) for _ in range(3)] for _ in range(3)]
        for j in range(3):
            k=2*ai[j]/ti;w=iv.exp(2*k*(x[j]-ci[j]));response=(w-1)/(w+1)
            values.append(2*(ci[j]+ai[j]*response)-x[j]);jac[j][j]=8*ai[j]*k*w/(w+1)**2-1
        return values,jac
    def all_sites(x):
        xx=[full.Dual(v,[iv.mpf(int(i==j)) for j in range(3)]) for i,v in enumerate(x)]
        sites=[[exact(t,iv) for t in p] for p in raw['points']]
        weights=[full.exp(-sum(full.square(xx[j]-p[j]) for j in range(3))/ti) for p in sites]
        total=sum(weights);m=[sum(w*p[j] for w,p in zip(weights,sites))/total for j in range(3)]
        v=[2*m[j]-xx[j] for j in range(3)]
        return [q.v for q in v],[q.d for q in v]
    field=all_sites if args.engine=='full' else binary
    f0,_=field(xi);_,J=field(box)
    inverse=[mp.nstr(1/(4*a[j]**2/tau*(1-mp.tanh(2*a[j]*(roots[j]-c[j])/tau)**2)-1),args.digits+10) for j in range(3)]
    # v'=4*a^2/tau*sech^2 - 1? m'=2*a^2/tau*sech^2, hence 2*m'-1.
    C=[exact(v,iv) for v in inverse]
    R=[[iv.mpf(i==j)-C[i]*J[i][j] for j in range(3)] for i in range(3)]
    K=[xi[i]-C[i]*f0[i]+sum(R[i][j]*iv.mpf([-r.b,r.b]) for j in range(3)) for i in range(3)]
    assert all(k.a>b.a and k.b<b.b for k,b in zip(K,box))
    norm=max(sum(abs(v).b for v in row) for row in R);assert norm<1
    # Strict diagonal dominance with positive diagonal proves invertibility
    # throughout this box, independently of the approximate preconditioner.
    margins=[J[i][i]-sum(abs(J[i][j]) for j in range(3) if i!=j) for i in range(3)]
    assert all(v.a>0 for v in margins)
    result=dict(status='passed',engine=args.engine,digits=args.digits,mpmath=mp.__version__,input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),params=raw['params'],points=raw['points'],center=center,radius=args.radius,box=list(map(enc,box)),field_at_center=list(map(enc,f0)),derivative=[[enc(v) for v in row] for row in J],inverse_diagonal=inverse,remainder=[[enc(v) for v in row] for row in R],krawczyk=list(map(enc,K)),contraction=enc(norm),diagonal_dominance=list(map(enc,margins)),cluster_sphere_center=[enc(v/2) for v in box],cluster_sphere_radius=enc(bi),scope='A unique zero s of v=2m-x in the tiny box, with nonsingular Dv. F is undefined at s. For every unit direction u the inverse-function theorem gives points tending to s with F(x) tending to s/2+beta*u. No membership of s in the numerical cloud closure or in the certified invariant subsystem is claimed.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    (args.out/'input.json').write_bytes(args.input.read_bytes())
    for name in [Path(__file__).name,'product_binary_interval.py','verify_product_coverings.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(dict(status='passed',engine=args.engine,center=[float(v) for v in roots],contraction=float(norm),dominance=[float(v.a) for v in margins])))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','out','mpmath-path']:p.add_argument('--'+name,type=Path,required=name!='mpmath-path')
    p.add_argument('--engine',choices=['factorized','full'],default='factorized');p.add_argument('--digits',type=int,default=160)
    p.add_argument('--radius',default='1e-50');main(p.parse_args())
