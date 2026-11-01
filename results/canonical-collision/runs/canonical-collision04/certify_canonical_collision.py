"""Interval box trapping and an unstable-graph reflection crossing.

The finite inequalities support the elementary graph-transform argument in
the accompanying report. This does not identify the entire attracting set
with a sampled portrait, nor certify asymptotic Lyapunov exponents.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from canonical_cube import exact, parameters, mean_one, model


def mm(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def mv(a,b): return [sum(x*y for x,y in zip(row,b)) for row in a]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps=250;mp.iv.dps=250;iv=mp.iv;begin=time.time()
    args.out.mkdir(parents=True,exist_ok=False)
    for name in [Path(__file__).name,'canonical_cube.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    save=lambda p,v:p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    pm=parameters(args.source,653,mp.mp);pi=parameters(args.source,653,iv)
    dec=lambda v:mp.nstr(v,260)
    upper=lambda v:iv.mpf(v.b)
    lower=lambda v:iv.mpf(v.a)
    bound=lambda v:upper(abs(v))
    sym=lambda r:iv.mpf([-r.b,r.b])
    hull=lambda a,b:iv.mpf([a.a,b.b])
    lo=exact(-1.65,iv);hi=exact(-.08,iv)
    # First pass bounds every positive scalar q=x-2m(x).
    cells=[];n=4096
    for k in range(n):
        z=hull(lo+(hi-lo)*k/n,lo+(hi-lo)*(k+1)/n)
        gl,_=mean_one(iv.mpf(z.a),pi,iv)
        gh,_=mean_one(iv.mpf(z.b),pi,iv)
        # The scalar posterior mean is increasing; endpoint evaluation avoids
        # an unnecessarily wide quotient interval in saturated regions.
        g=hull(gl,gh);q=z-2*g
        assert q.a>0
        cells.append((z,g,q))
    qlo=min((q.a for _,_,q in cells));qhi=max((q.b for _,_,q in cells))
    bounds=[];beta=pi[0]
    for z,g,q in cells:
        # -q/sqrt(q^2+s) decreases in q>0 and increases in s>0.
        umin=-upper(q)/iv.sqrt(upper(q)**2+2*iv.mpf(qlo)**2)
        umax=-lower(q)/iv.sqrt(lower(q)**2+2*iv.mpf(qhi)**2)
        a=beta*lower(z)+beta*umin+(1-2*beta)*upper(g)
        b=beta*upper(z)+beta*umax+(1-2*beta)*lower(g)
        bounds.append(hull(a,b))
    image=hull(iv.mpf(min(v.a for v in bounds)),iv.mpf(max(v.b for v in bounds)))
    assert image.a>lo.b and image.b<hi.a
    trap=dict(box=[str(hull(lo,hi))]*3,subintervals=n,q_range=str(hull(iv.mpf(qlo),iv.mpf(qhi))),image_coordinate=str(image),strict_margins=[str(image-lo),str(hi-image)],minimum_radius=str(iv.sqrt(3)*iv.mpf(qlo)),method='Monotonic component bound using all 4096 scalar interval cells; other-coordinate squared radii range from 2*q_min^2 to 2*q_max^2.')
    print(json.dumps(dict(stage='trapping',image=str(image),radius=trap['minimum_radius'])),flush=True)
    # Refine and enclose the reflection-plane fixed point in two variables.
    def equation(s,t):
        y=model([s,s,t],pm,mp.mp,False)
        return y[0]-s,y[2]-t
    s,t=mp.findroot(equation,('-.5744','-.1238'),tol=mp.mpf('1e-235'))
    def point(v):
        result=iv.mpf(v)
        assert result._mpi_[0]==v._mpf_==result._mpi_[1]
        return result
    center=[point(s),point(t)];rad=iv.mpf('1e-220')
    rootbox=[hull(v-rad,v+rad) for v in center]
    Pbox=[rootbox[0],rootbox[0],rootbox[1]]
    y0,J0,_=model([center[0],center[0],center[1]],pi,iv)
    _,Jb,_=model(Pbox,pi,iv)
    H=[[Jb[0][0]+Jb[0][1]-1,Jb[0][2]],[Jb[2][0]+Jb[2][1],Jb[2][2]-1]]
    _,Jn,_=model([s,s,t],pm,mp.mp)
    C=mp.inverse(mp.matrix([[Jn[0][0]+Jn[0][1]-1,Jn[0][2]],[Jn[2][0]+Jn[2][1],Jn[2][2]-1]]))
    Ci=[[point(C[i,j]) for j in range(2)] for i in range(2)]
    CH=mm(Ci,H);R=[[iv.mpf(i==j)-CH[i][j] for j in range(2)] for i in range(2)]
    correction=mv(Ci,[y0[0]-center[0],y0[2]-center[1]])
    remainder=mv(R,[rootbox[i]-center[i] for i in range(2)])
    K=[center[i]-correction[i]+remainder[i] for i in range(2)]
    assert all(k.a>b.a and k.b<b.b for k,b in zip(K,rootbox))
    kr=max(sum(bound(v) for v in row).b for row in R);assert kr<1
    root_record=dict(center=[dec(s),dec(t)],radius='1e-220',box=[str(v) for v in rootbox],box_binary=[[list(v._mpi_[0]),list(v._mpi_[1])] for v in rootbox],krawczyk=[str(v) for v in K],contraction=str(kr),symmetry='The full fixed point is (s,s,t).')
    # Exact dyadic coordinates: antisymmetric u, then two stable eigendirections.
    stable=np.array([[float(Jn[0][0]+Jn[0][1]),float(Jn[0][2])],[float(Jn[2][0]+Jn[2][1]),float(Jn[2][2])]])
    eigen,S=np.linalg.eig(stable);assert np.max(abs(eigen))<1 and np.max(abs(S.imag))==0
    S=S.real;E=np.array([[1.,S[0,0],S[0,1]],[-1.,S[0,0],S[0,1]],[0.,S[1,0],S[1,1]]])
    Ef=[[Fraction.from_float(float(v)) for v in row] for row in E]
    a,b=Ef[0][1:];c,d=Ef[2][1:];det=a*d-b*c
    Ei_f=[[Fraction(1,2),Fraction(-1,2),Fraction(0)],[d/(2*det),d/(2*det),-b/det],[-c/(2*det),-c/(2*det),a/det]]
    assert mm(Ef,Ei_f)==[[Fraction(i==j) for j in range(3)] for i in range(3)]
    cv=lambda q,ctx:ctx.mpf(q.numerator)/ctx.mpf(q.denominator)
    EI=[[cv(q,iv) for q in row] for row in Ef];II=[[cv(q,iv) for q in row] for row in Ei_f]
    EM=[[cv(q,mp.mp) for q in row] for row in Ef]
    r=iv.mpf('2e-80');L=iv.mpf('1e-70')
    B=[sym(r),sym(r*L),sym(r*L)];offset=mv(EI,B)
    physical=[Pbox[i]+offset[i] for i in range(3)]
    _,JB,_=model(physical,pi,iv);G=mm(mm(II,JB),EI)
    aa=lower(G[0][0]);bb=upper(sum(bound(v) for v in G[0][1:]));cc=iv.mpf(max(bound(G[i][0]).b for i in [1,2]));dd=iv.mpf(max(sum(bound(v) for v in G[i][1:]).b for i in [1,2]))
    expansion=aa-bb*L;lip=(cc+dd*L)/expansion;vertical=dd+(cc+dd*L)*bb/expansion
    assert expansion.a>1 and lip.b<L.a and vertical.b<1 and (cc+dd*L).b<L.a
    graph=dict(radius='2e-80',lipschitz='1e-70',E_fraction=[[str(q) for q in row] for row in Ef],inverse_fraction=[[str(q) for q in row] for row in Ei_f],derivative_blocks=dict(a=str(aa),b=str(bb),c=str(cc),d=str(dd)),horizontal_expansion=str(expansion),new_lipschitz_bound=str(lip),graph_contraction=str(vertical),stable_boundary_ratio=str((cc+dd*L)/L),root_in_Q=all(v.a>lo.b and v.b<hi.a for v in Pbox))
    assert graph['root_in_Q']
    print(json.dumps(dict(stage='local graph',expansion=str(expansion),contraction=str(vertical))),flush=True)
    search=json.loads((args.search/'report.json').read_text(encoding='utf-8'));selected=search['selected'];eps_string=selected['epsilon'];eps=iv.mpf(eps_string);em=mp.mpf(eps_string)
    assert eps.a>0 and eps.b<r.a
    # A true graph point has stable coordinates bounded by L*epsilon.
    xc=[s+em,s-em,t];initial=[Pbox[i]+mv(EI,[eps,sym(L*eps),sym(L*eps)])[i] for i in range(3)]
    nominal=[point(v) for v in xc];error=mv(II,[initial[i]-nominal[i] for i in range(3)])
    delta=[sym(bound(v)) for v in error];trail=[];minimum_radius=mp.inf
    for nstep in range(1,selected['crossing_iteration']+1):
        nominal=[point(v) for v in xc];offset=mv(EI,delta)
        xbox=[nominal[i]+offset[i] for i in range(3)]
        assert all(v.a>lo.b and v.b<hi.a for v in xbox)
        _,JX,rx=model(xbox,pi,iv);assert rx.a>0
        minimum_radius=min(minimum_radius,mp.mpf(rx.a))
        image0=model(nominal,pi,iv,False);yc=model(xc,pm,mp.mp,False)
        center_error=mv(II,[image0[i]-point(yc[i]) for i in range(3)])
        propagation=mv(mm(mm(II,JX),EI),delta)
        delta=[sym(bound(center_error[i]+propagation[i])) for i in range(3)];xc=yc
        outbox=[point(xc[i])+mv(EI,delta)[i] for i in range(3)]
        assert all(v.a>lo.b and v.b<hi.a for v in outbox)
        diff=outbox[0]-outbox[1]
        trail.append(dict(iteration=nstep,maximum_coordinate_error=mp.nstr(max(mp.mpf(abs(v).b) for v in delta),20),difference_sign=1 if diff.a>0 else -1 if diff.b<0 else 0))
        if nstep%100==0:print(json.dumps(dict(stage='propagation',iteration=nstep,error=trail[-1]['maximum_coordinate_error'])),flush=True)
    assert diff.b<0
    endpoint=dict(epsilon=eps_string,iteration=selected['crossing_iteration'],final_box=[str(v) for v in outbox],final_x_minus_y=str(diff),maximum_final_coordinate_error=trail[-1]['maximum_coordinate_error'],minimum_radius=dec(minimum_radius),trail=trail)
    result=dict(status='passed',source_sha256=sha(args.source),search_sha256=sha(args.search/'report.json'),precision_digits=250,mpmath_version=mp.__version__,trapping=trap,fixed_point=root_record,unstable_graph=graph,endpoint=endpoint,seconds=time.time()-begin,scope='Strict trapping box and validated reflection-plane crossing of a local unstable graph. The topological argument yields intrinsic noninjectivity on K=intersection F^n(Q). Does not identify K with the numerical portrait or establish hyperchaos, entropy or historical priority.')
    save(args.out/'report.json',result)
    print(json.dumps(dict(status='passed',steps=endpoint['iteration'],final_difference=str(diff),seconds=time.time()-begin)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','search','out','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
