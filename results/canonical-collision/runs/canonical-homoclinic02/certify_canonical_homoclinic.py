"""Validate a reflection-plane homoclinic return and a transverse secant cone.

The unstable graph is supplied by the collision certificate. Its small vertical
uncertainty is propagated along a parameter bracket, never replaced by a line.
The optional independent mode uses eight Gaussian weights and automatic
derivatives. No entropy or horseshoe theorem is invoked here.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from canonical_cube import parameters, model
from verify_canonical_collision import posterior_reflections, matmul, matvec


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps=mp.iv.dps=300 if args.independent else 250
    iv=mp.iv;start=time.time()
    args.out.mkdir(parents=True,exist_ok=False)
    for name in [Path(__file__).name,'canonical_cube.py','verify_canonical_collision.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    cert=json.loads(args.certificate.read_text(encoding='utf-8'))
    assert cert['status']=='passed' and cert['source_sha256']==sha(args.source)
    pm=parameters(args.source,653,mp.mp);pi=parameters(args.source,653,iv)
    mn=lambda x:model(x,pm,mp.mp,False)
    if args.independent:
        with np.load(args.source) as data:
            sites=data['points'][653,:data['counts'][653]].copy();raw=data['params'][653].copy()
        def exact(v,ctx):
            a,b=float(v).as_integer_ratio();return ctx.mpf(a)/ctx.mpf(b)
        points=[[exact(v,iv) for v in row] for row in sites]
        beta,tau=[exact(v,iv) for v in raw]
        mi=lambda x,derivative=True:posterior_reflections(x,points,beta,tau,iv,derivative)
    else:mi=lambda x,derivative=True:model(x,pi,iv,derivative)
    def point(v):
        p=iv.mpf(v);assert p._mpi_[0]==v._mpf_==p._mpi_[1];return p
    def upper(v):return iv.mpf(v.b)
    def magnitude(v):return upper(abs(v))
    def sym(v):return iv.mpf([-v.b,v.b])
    def hull(a,b):return iv.mpf([a.a,b.b])
    lo=iv.mpf(float(-1.65));hi=iv.mpf(float(-.08))
    Efra=[[Fraction(v) for v in row] for row in cert['unstable_graph']['E_fraction']]
    Ifra=[[Fraction(v) for v in row] for row in cert['unstable_graph']['inverse_fraction']]
    assert matmul(Efra,Ifra)==[[Fraction(i==j) for j in range(3)] for i in range(3)]
    cv=lambda mat:[[iv.mpf(v.numerator)/iv.mpf(v.denominator) for v in row] for row in mat]
    E,I=cv(Efra),cv(Ifra)
    P0=[iv.mpf([mp.mpf(tuple(v[0])),mp.mpf(tuple(v[1]))]) for v in cert['fixed_point']['box_binary']]
    P=[P0[0],P0[0],P0[1]]
    s,t=[mp.mpf(v) for v in cert['fixed_point']['center']]
    L=iv.mpf(cert['unstable_graph']['lipschitz']);r=iv.mpf(cert['unstable_graph']['radius'])
    base=mp.mpf(cert['endpoint']['epsilon']);N=cert['endpoint']['iteration']
    def nominal(c):
        x=[s+c*base,s-c*base,t]
        for _ in range(N):x=mn(x)
        return x
    def difference(c):
        x=nominal(c);return x[0]-x[1]
    root=mp.findroot(difference,('.99','1'),tol=mp.mpf('1e-160'),maxsteps=100)
    middle=mp.nstr(root*base,200)
    eps=iv.mpf(middle);width=iv.mpf('1e-115')
    bracket=hull(eps-width,eps+width)
    assert bracket.a>0 and bracket.b<r.a
    print(json.dumps(dict(stage='located bracket',factor=mp.nstr(root,30))),flush=True)

    def propagate(seed,derivatives=False):
        em=(mp.mpf(seed.a)+mp.mpf(seed.b))/2
        xc=[s+em,s-em,t]
        off=matvec(E,[seed,sym(L*upper(seed)),sym(L*upper(seed))])
        delta=[sym(magnitude(v)) for v in matvec(I,[P[i]+off[i]-point(xc[i]) for i in range(3)])]
        M=[[iv.mpf(i==j) for j in range(3)] for i in range(3)]
        minimum=mp.inf;trail=[]
        for step in range(N):
            center=[point(v) for v in xc];off=matvec(E,delta)
            box=[center[i]+off[i] for i in range(3)]
            assert all(lo.b<v.a and v.b<hi.a for v in box)
            _,J,axis=mi(box);minimum=min(minimum,mp.mpf(axis.a));assert axis.a>0
            G=matmul(matmul(I,J),E)
            if derivatives:M=matmul(G,M)
            fi=mi(center,False);yc=mn(xc)
            roundoff=matvec(I,[fi[i]-point(yc[i]) for i in range(3)])
            advanced=matvec(G,delta)
            delta=[sym(magnitude(roundoff[i]+advanced[i])) for i in range(3)];xc=yc
            if step>=N-40:trail.append(dict(iteration=step+1,center=[mp.nstr(v,70) for v in xc],coordinate_error=[str(v) for v in delta]))
        off=matvec(E,delta);final=[point(xc[i])+off[i] for i in range(3)]
        assert all(lo.b<v.a and v.b<hi.a for v in final)
        return final,M,minimum,trail

    left,_,_,_=propagate(iv.mpf(bracket.a))
    right,_,_,_=propagate(iv.mpf(bracket.b))
    dl=left[0]-left[1];dr=right[0]-right[1]
    print(json.dumps(dict(stage='bracket signs',left=str(dl),right=str(dr))),flush=True)
    assert dl.b<0 and dr.a>0
    whole,M,minradius,trail=propagate(bracket,True)
    # Every secant of the Lipschitz graph has slope (1,v), |v| <= L.
    # Matrix interval products over all intermediate boxes therefore bound
    # the endpoint secant's antisymmetric component by the mean-value theorem.
    slope=2*matvec(M,[iv.mpf(1),sym(L),sym(L)])[0]
    assert slope.a>0
    print(json.dumps(dict(stage='crossing validated',slope=mp.nstr(mp.mpf(slope.a),20))),flush=True)
    # Restrict the crossing enclosure to x=y exactly. Its subsequent orbit
    # stays in the fixed plane by equivariance.
    common=iv.mpf([max(whole[0].a,whole[1].a),min(whole[0].b,whole[1].b)])
    crossing=[common,common,whole[2]]
    # A contracting square in stable coordinates is centered at the exact
    # fixed point. Its derivative bound proves invariance and convergence.
    stable_radius=iv.mpf('1e-4')
    off=matvec(E,[iv.mpf(0),sym(stable_radius),sym(stable_radius)])
    local=[P[i]+off[i] for i in range(3)]
    _,J,_=mi(local);G=matmul(matmul(I,J),E)
    contraction=max(sum(magnitude(v) for v in row[1:]).b for row in G[1:])
    assert contraction<1
    x=crossing;plane_trail=[]
    for step in range(101):
        assert all(lo.b<v.a and v.b<hi.a for v in x)
        coords=matvec(I,[x[i]-P[i] for i in range(3)])
        if max(magnitude(coords[j]).b for j in [1,2])<stable_radius.a:
            break
        y=mi(x,False)
        # Same common input; both components are equal mathematically. The
        # intersection of their interval enclosures remains a valid enclosure.
        common=iv.mpf([max(y[0].a,y[1].a),min(y[0].b,y[1].b)])
        x=[common,common,y[2]]
        plane_trail.append(dict(iteration=step+1,box=[str(v) for v in x]))
    else:raise AssertionError('No entry into the contracting plane square')
    result=dict(status='passed',source_sha256=sha(args.source),certificate_sha256=sha(args.certificate),independent_formula=args.independent,precision_digits=mp.mp.dps,seconds=time.time()-start,
        bracket=dict(center=middle,width='1e-115',binary=[list(v) for v in bracket._mpi_],left_difference=str(dl),right_difference=str(dr),endpoint_secant_slope=str(slope)),
        crossing=dict(iteration=N,box=[str(v) for v in crossing],minimum_radius=mp.nstr(minradius,70),last40=trail),
        return_to_saddle=dict(iterations=step,final_box=[str(v) for v in x],stable_coordinates=[str(v) for v in coords[1:]],local_radius='1e-4',local_contraction=str(contraction),trail=plane_trail),
        scope='A nontrivial homoclinic orbit in K with a validated reflection-plane crossing and positive secant slope across its parameter bracket. Requires the certified local unstable graph. No horseshoe, entropy, transitivity, hyperchaos, equality with the numerical portrait or historical novelty is claimed.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status='passed',return_iterations=step,seconds=result['seconds'],contraction=str(contraction))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','certificate','out','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--independent',action='store_true');main(p.parse_args())
