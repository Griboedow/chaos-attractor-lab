"""Ordinary folds and two preimages of a certified normalization singularity.

No membership in an invariant set is inferred from proximity to a cloud.
"""
import argparse, hashlib, json, sys, time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from product_binary_interval import parameters, model, mm, mv
from product_second_derivatives import binary_second,full_second,determinant_gradient,cross
from certify_product_cycles import det


def main(args):
    sys.path.append(str(args.mpmath_path));import mpmath as mp
    mp.mp.dps=mp.iv.dps=args.digits;iv=mp.iv;begin=time.time();args.out.mkdir(parents=True,exist_ok=False)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    exact=lambda x,ctx=iv:ctx.mpf(Q(x).numerator)/Q(x).denominator
    enc=lambda v:dict(decimal=str(v),binary=[list(t) for t in v._mpi_])
    dec=lambda v:mp.nstr(v,args.digits+10)
    save=lambda name,obj:(args.out/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8',newline='\n')
    sources={}
    for name in [Path(__file__).name,'product_binary_interval.py','product_second_derivatives.py','certify_product_cycles.py']:
        p=Path(__file__).with_name(name);(args.out/name).write_bytes(p.read_bytes());sources[name]=sha(p)
    raw=json.loads(args.input.read_text(encoding='utf-8'));singular=json.loads(args.singular.read_text(encoding='utf-8'))
    assert singular['status']=='passed' and singular['input_sha256']==sha(args.input)
    par=parameters(raw['params'],np.array(raw['points']),mp.mp);pari=parameters(raw['params'],np.array(raw['points']),iv)
    field=lambda x:full_second(x,raw,iv) if args.engine=='full' else binary_second(x,pari,iv)
    numeric=lambda x:binary_second(x,par,mp.mp)
    def decode(t):
        sign,mantissa,exponent,bits=t;assert bits==int(mantissa).bit_length()
        return Q((-1 if sign else 1)*mantissa)*Q(2)**exponent
    S=[iv.mpf([exact(decode(v['binary'][0])).a,exact(decode(v['binary'][1])).b]) for v in singular['krawczyk']]
    target=[exact(v,mp.mp) for v in singular['center']]
    proposals=json.loads(args.proposals.read_text(encoding='utf-8'))['singular_preimages']
    # Retain both previously declared proposals inside the rendering coordinate box.
    proposed=[p for p in proposals if p['in_coordinate_box']];assert len(proposed)==2
    preimages=[];radius=exact(args.radius)
    for proposal in proposed:
        root=mp.findroot(lambda *x:[a-b for a,b in zip(numeric(x)[0],target)],proposal['state'],tol=mp.mpf(10)**(-args.digits+25),maxsteps=30)
        center=list(map(dec,root));xi=[exact(v) for v in center];box=[v+iv.mpf([-radius.b,radius.b]) for v in xi]
        y0,_,_,r0=field(xi);_,J,_,r=field(box);assert r.a>0
        Ctext=[[dec(v) for v in row] for row in (mp.matrix(numeric(root)[1])**-1).tolist()];C=[[exact(v) for v in row] for row in Ctext]
        residual=[a-b for a,b in zip(y0,S)];R=[[int(i==j)-v for j,v in enumerate(row)] for i,row in enumerate(mm(C,J))]
        correction=mv(C,residual);rest=mv(R,[iv.mpf([-radius.b,radius.b])]*3);K=[xi[i]-correction[i]+rest[i] for i in range(3)]
        norm=max(sum(abs(v).b for v in row) for row in R)
        assert norm<1 and all(k.a>b.a and k.b<b.b for k,b in zip(K,box)) and abs(det(C)).a>0
        preimages.append(dict(center=center,radius=args.radius,box=list(map(enc,box)),target=list(map(enc,S)),residual=list(map(enc,residual)),
                              derivative=[[enc(v) for v in row] for row in J],preconditioner=Ctext,remainder=[[enc(v) for v in row] for row in R],
                              krawczyk=list(map(enc,K)),contraction=enc(norm),normalization=enc(r),determinant=enc(det(J))))
    with np.load(args.critical) as f:points=f['points'];tags=f['tags'];distances=f['reference_distances'][0]
    centers=np.array([float(v) for v in par[2]]);folds=[];controls=[]
    for axis in range(3):
        for side in [-1,1]:
            mask=(tags==axis)&((points[:,axis]-centers[axis])*side>0)
            index=np.flatnonzero(mask)[np.argmin(distances[mask])];proposal=points[index]
            q=[exact(float(v),mp.mp) for v in proposal]
            def line(t):
                x=list(q);x[axis]=t;return det(numeric(x)[1])
            q[axis]=mp.findroot(line,q[axis],solver='newton',tol=mp.mpf(10)**(-args.digits+25))
            center=[dec(v) if j==axis else str(Q(float(proposal[j]))) for j,v in enumerate(q)]
            xi=[exact(v) for v in center];box=xi.copy();box[axis]+=iv.mpf([-radius.b,radius.b])
            _,J,H,r=field(box);assert r.a>0
            gradient=determinant_gradient(J,H);kernel=cross(J[0],J[1]);fold=sum(a*b for a,b in zip(gradient,kernel))
            endpoints=[]
            for s in [-1,1]:
                x=xi.copy();x[axis]+=s*radius;endpoints.append(det(field(x)[1]))
            assert (endpoints[0].b<0<endpoints[1].a) or (endpoints[1].b<0<endpoints[0].a)
            assert abs(gradient[axis]).a>0 and abs(fold).a>0 and any(abs(v).a>0 for v in kernel)
            # Agreement of independent high-precision expressions at the proposal.
            A=numeric(q);B=full_second(q,raw,mp.mp)
            errors=[max(abs(A[0][i]-B[0][i]) for i in range(3)),max(abs(A[1][i][j]-B[1][i][j]) for i in range(3) for j in range(3)),
                    max(abs(A[2][i][j][k]-B[2][i][j][k]) for i in range(3) for j in range(3) for k in range(3))]
            assert max(errors)<mp.mpf(10)**(-args.digits+20)
            controls.append([dec(v) for v in errors])
            folds.append(dict(axis=axis,side=side,source_index=int(index),source_cloud_distance=float(distances[index]),center=center,radius=args.radius,
                              box=list(map(enc,box)),normalization=enc(r),endpoint_determinants=list(map(enc,endpoints)),
                              derivative=[[enc(v) for v in row] for row in J],hessian=[[[enc(v) for v in row] for row in mat] for mat in H],
                              determinant_gradient=list(map(enc,gradient)),kernel=list(map(enc,kernel)),fold_transversality=enc(fold)))
    report=dict(status='passed',engine=args.engine,digits=args.digits,sources=sources,input_sha256=sha(args.input),singular_sha256=sha(args.singular),
                proposals_sha256=sha(args.proposals),critical_sha256=sha(args.critical),preimages=preimages,folds=folds,independent_expression_errors=controls,seconds=time.time()-begin,
                scope='Two smooth local inverse branches lead exactly to the certified normalization singularity; six ordinary rank-two folds exist on fixed-coordinate lines. These local facts do not prove membership in any invariant set, trapping, or historical novelty.')
    save('report.json',report);print(json.dumps(dict(status=report['status'],preimages=len(preimages),folds=len(folds),seconds=report['seconds'])))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','singular','proposals','critical','out','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--engine',choices=['factorized','full'],default='factorized');p.add_argument('--digits',type=int,default=120);p.add_argument('--radius',default='1e-35')
    main(p.parse_args())
