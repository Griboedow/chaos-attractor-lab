"""Validate an unstable-graph escape from the compact covering atlas.

An escape proves lack of Lyapunov stability for invariant subsets of that
atlas containing cycle 14. It does not disprove a larger attracting set.
"""
import argparse, hashlib, json, sys, time
from fractions import Fraction
from pathlib import Path
import numpy as np
from product_binary_interval import parameters, model, inverse, mm, mv
from certify_product_cycles import iterate, det


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    import verify_product_coverings as full
    full.mp = mp; mp.mp.dps = mp.iv.dps = args.digits; iv = mp.iv
    args.out.mkdir(parents=True, exist_ok=False); begin = time.time()
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    save = lambda p,v: p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')
    enc = lambda v: dict(decimal=str(v),binary=[list(t) for t in v._mpi_])
    exact = lambda x,ctx=iv: ctx.mpf(Fraction(x).numerator)/Fraction(x).denominator
    sym = lambda x: iv.mpf([-abs(x).b,abs(x).b])
    upper = lambda x: iv.mpf(abs(x).b)
    cv = lambda mat,ctx=iv: [[exact(v,ctx) for v in row] for row in mat]
    raw = json.loads(args.windows.read_text(encoding='utf-8'))
    prior = json.loads(args.cycles.read_text(encoding='utf-8'))
    assert prior['status']=='passed' and prior['windows_sha256']==sha(args.windows)
    cycle = next(x for x in prior['cycles'] if x['id']==14)
    assert cycle['period']==4
    par = parameters(raw['params'],np.array(raw['points']),mp.mp)
    pari = parameters(raw['params'],np.array(raw['points']),iv)
    independent = full.Model(raw)
    field = independent.value_jacobian if args.engine=='full' else lambda x:model(x,pari,iv)
    numeric = lambda x:model(x,par,mp.mp)
    center_text = cycle['center_decimal']
    center = [exact(x) for x in center_text]; numeric_center = [exact(x,mp.mp) for x in center_text]
    radius = exact(args.root_radius)
    box = [x+sym(radius) for x in center]
    y0,_,_,_ = iterate(center,4,iv,field)
    _,J,_,_ = iterate(box,4,iv,field)
    H = [[J[i][j]-int(i==j) for j in range(3)] for i in range(3)]
    C = cv(cycle['preconditioner_decimal'])
    R = [[int(i==j)-v for j,v in enumerate(row)] for i,row in enumerate(mm(C,H))]
    correction = mv(C,[a-b for a,b in zip(y0,center)])
    remainder = mv(R,[sym(radius)]*3)
    K = [center[i]-correction[i]+remainder[i] for i in range(3)]
    assert all(k.a>b.a and k.b<b.b for k,b in zip(K,box))
    root_contraction = max(sum(upper(v) for v in row).b for row in R)
    assert root_contraction<1 and (det(C).a>0 or det(C).b<0)
    # An exact rational coordinate chart, shared by both expression engines.
    _,M,_,_ = iterate(numeric_center,4,mp.mp,numeric)
    vals, E = np.linalg.eig(np.array(M,dtype=float))
    order = np.argsort(-np.abs(vals)); E = E[:,order]; vals = vals[order]
    assert np.max(abs(E.imag))==0 and abs(vals[0])>1 and np.max(abs(vals[1:]))<1
    Ef = [[Fraction(float(v)) for v in row] for row in E.real]
    Ei = inverse(Ef); assert mm(Ef,Ei)==[[Fraction(i==j) for j in range(3)] for i in range(3)]
    EI, II, EM = cv(Ef),cv(Ei),cv(Ef,mp.mp)
    r, L = exact(args.graph_radius),exact(args.slope)
    offset = mv(EI,[sym(r),sym(r*L),sym(r*L)])
    graph_box = [K[i]+offset[i] for i in range(3)]
    _,JB,_,radii = iterate(graph_box,4,iv,field)
    chart_J = mm(mm(II,JB),EI)
    assert chart_J[0][0].b<0 or chart_J[0][0].a>0
    a = iv.mpf(abs(chart_J[0][0]).a)
    b = iv.mpf(sum(upper(v) for v in chart_J[0][1:]).b)
    c = iv.mpf(max(upper(chart_J[i][0]).b for i in [1,2]))
    d = iv.mpf(max(sum(upper(v) for v in chart_J[i][1:]).b for i in [1,2]))
    expansion = a-b*L
    slope = (c+d*L)/expansion
    contraction = d+(c+d*L)*b/expansion
    assert expansion.a>1 and slope.b<L.a and contraction.b<1 and (c+d*L).b<L.a
    graph = dict(radius=args.graph_radius,L=args.slope,E=[[str(v) for v in row] for row in Ef],
                 inverse=[[str(v) for v in row] for row in Ei],physical_box=list(map(enc,graph_box)),
                 jacobian=[[enc(v) for v in row] for row in chart_J],
                 a=enc(a),b=enc(b),c=enc(c),d=enc(d),expansion=enc(expansion),
                 slope=enc(slope),contraction=enc(contraction),stable_boundary=enc(c+d*L))
    print(json.dumps(dict(stage='graph',expansion=float(expansion.a),contraction=float(contraction.b),slope=float(slope.b))),flush=True)
    # Enclose X = union_i union_{j=0}^3 F^j(N_i). Each complete affine window
    # is included in its exact coordinate hull; mean-value intersection tightens images.
    atlas=[]
    for ident,w in enumerate(raw['exact_windows']):
        q=[exact(v) for v in raw['centers'][ident]]
        S=cv(w['matrix']); spread=mv(S,[iv.mpf([-1,1])]*3)
        cell=[q[j]+spread[j] for j in range(3)]
        for phase in range(4):
            atlas.append(dict(window=ident,phase=phase,box=list(map(enc,cell))))
            if phase==3:break
            natural,derivative,rr=field(cell); assert rr.a>0
            mid=[iv.mpf((mp.mpf(v.a)+mp.mpf(v.b))/2) for v in cell]
            base,_,_=field(mid); change=mv(derivative,[v-t for v,t in zip(cell,mid)])
            mean=[t+e for t,e in zip(base,change)]
            cell=[iv.mpf([max(t.a,e.a),min(t.b,e.b)]) for t,e in zip(natural,mean)]
    atlas_boxes=[]
    for row in atlas:
        vals=[]
        for v in row['box']:
            def decode(t):
                sign,mantissa,exponent,_=t
                return Fraction((-1 if sign else 1)*mantissa)*(Fraction(2)**exponent)
            lo,hi=[exact(decode(t)) for t in v['binary']]
            vals.append(iv.mpf([lo.a,hi.b]))
        atlas_boxes.append(vals)
    eps=exact(args.epsilon); assert eps.a>0 and eps.b<r.a
    displacement=mv(EI,[eps,sym(L*eps),sym(L*eps)])
    initial=[K[i]+displacement[i] for i in range(3)]
    nominal=[numeric_center[i]+v for i,v in enumerate(mv(EM,[exact(args.epsilon,mp.mp),mp.mpf(0),mp.mpf(0)]))]
    delta=[sym(upper(v)) for v in mv(II,[initial[i]-iv.mpf(nominal[i]) for i in range(3)])]
    trail=[]; separated=None
    for n in range(args.max_steps+1):
        off=mv(EI,delta); current=[iv.mpf(nominal[i])+off[i] for i in range(3)]
        separations=[]
        for abox in atlas_boxes:
            options=[(i,1,iv.mpf(current[i].a)-iv.mpf(abox[i].b)) for i in range(3)]
            options += [(i,-1,iv.mpf(abox[i].a)-iv.mpf(current[i].b)) for i in range(3)]
            axis,sign,gap=max(options,key=lambda t:mp.mpf(t[2].a))
            separations.append(dict(axis=axis,sign=sign,gap=enc(gap)))
        mingap=min(mp.mpf(s['gap']['decimal'].split(',')[0].lstrip('[')) for s in separations)
        # Numeric display above is diagnostic; acceptance uses the actual intervals.
        lower_gaps=[]
        for s,abox in zip(separations,atlas_boxes):
            j=s['axis']; g=current[j]-abox[j] if s['sign']==1 else abox[j]-current[j]
            lower_gaps.append(mp.mpf(g.a))
        mingap=min(lower_gaps)
        trail.append(dict(iteration=n,center=[mp.nstr(v,args.digits+10) for v in nominal],
                          chart_error=list(map(enc,delta)),box=list(map(enc,current)),minimum_atlas_gap=mp.nstr(mingap,30)))
        if mingap>mp.mpf(args.minimum_gap):
            separated=dict(iteration=n,box=list(map(enc,current)),separations=separations,minimum_gap=mp.nstr(mingap,40));break
        if n==args.max_steps:break
        _,jac,rr=field(current); assert rr.a>0
        image0,_,_=field([iv.mpf(v) for v in nominal]); nxt,_,_=numeric(nominal)
        center_error=mv(II,[image0[i]-iv.mpf(nxt[i]) for i in range(3)])
        prop=mv(mm(mm(II,jac),EI),delta)
        delta=[sym(upper(center_error[i]+prop[i])) for i in range(3)]; nominal=nxt
    sources={}
    for name in [Path(__file__).name,'product_binary_interval.py','certify_product_cycles.py','verify_product_coverings.py']:
        p=Path(__file__).with_name(name);(args.out/name).write_bytes(p.read_bytes());sources[name]=sha(p)
    save(args.out/'atlas.json',atlas);save(args.out/'trajectory.json',trail)
    report=dict(status='passed' if separated else 'no_separated_endpoint',engine=args.engine,digits=args.digits,
                windows_sha256=sha(args.windows),cycles_sha256=sha(args.cycles),sources=sources,
                root=dict(center_decimal=center_text,radius=args.root_radius,box=list(map(enc,box)),
                          preconditioner=cycle['preconditioner_decimal'],derivative=[[enc(v) for v in row] for row in H],
                          residual=list(map(enc,[y0[i]-center[i] for i in range(3)])),
                          krawczyk=list(map(enc,K)),remainder=[[enc(v) for v in row] for row in R],contraction=str(root_contraction)),
                graph=graph,epsilon=args.epsilon,atlas_sha256=sha(args.out/'atlas.json'),
                trajectory_sha256=sha(args.out/'trajectory.json'),separated=separated,seconds=time.time()-begin,
                scope='A true local unstable-graph point of period-four cycle 14 has a validated forward image outside the complete covering atlas X. This implies failure of Lyapunov stability for any invariant subset of X containing that cycle, using backward convergence. It does not exclude a larger attractor or identify the numerical cloud.')
    save(args.out/'report.json',report)
    print(json.dumps(dict(status=report['status'],iteration=separated['iteration'] if separated else None,gap=separated['minimum_gap'] if separated else None,seconds=report['seconds'])),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['windows','cycles','out','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--engine',choices=['factorized','full'],default='factorized')
    p.add_argument('--digits',type=int,default=120);p.add_argument('--root-radius',default='1e-50')
    p.add_argument('--graph-radius',default='2e-25');p.add_argument('--slope',default='1e-10');p.add_argument('--epsilon',default='1e-25')
    p.add_argument('--minimum-gap',default='1e-4');p.add_argument('--max-steps',type=int,default=160)
    main(p.parse_args())
