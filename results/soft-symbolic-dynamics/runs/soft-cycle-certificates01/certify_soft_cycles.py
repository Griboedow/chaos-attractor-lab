"""Krawczyk, least-period, real-multiplier and trapping-membership certificates."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np


def model(x,beta,tau,points,period,ctx):
    I=eye(ctx);M=eye(ctx);orbit=[];axes=[]
    for k in range(period):
        orbit.append(x)
        logits=[(2*sum(p[j]*x[j] for j in range(3))-sum(a*a for a in p))/tau for p in points]
        weights=[ctx.exp(a) for a in logits];total=sum(weights);weights=[a/total for a in weights]
        m=[sum(weights[t]*points[t][j] for t in range(len(points))) for j in range(3)]
        D=[[2/tau*sum(weights[t]*(points[t][i]-m[i])*(points[t][j]-m[j]) for t in range(len(points))) for j in range(3)] for i in range(3)]
        rb=[2*m[j]-x[j] for j in range(3)];radius=ctx.sqrt(rb[0]**2+rb[1]**2);axes.append(radius)
        a=[rb[0]/radius,rb[1]/radius,ctx.mpf(0)];tangent=[-a[1],a[0],ctx.mpf(0)]
        J=[[beta*I[i][j]+(1-2*beta)*D[i][j]+beta/radius*tangent[i]*sum(tangent[t]*(2*D[t][j]-I[t][j]) for t in range(3)) for j in range(3)] for i in range(3)]
        M=matmul(J,M);x=[beta/2*(2*a[j]-rb[j]+x[j])+(1-beta)*m[j] for j in range(3)]
    return x,M,orbit,axes


def main(args):
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    global exactfloat,eye,matmul,matvec,det
    from certify_raar_cycles import exactfloat,eye,matmul,matvec,det
    args.out.mkdir(parents=True,exist_ok=False);mp.mp.dps=200;mp.iv.dps=200;iv=mp.iv;begin=time.time();rows=[]
    with np.load(args.trap) as f:trap_params=f['params'];trap_points=f['points'];origin=f['origin'];shape=f['shape'];h=float(f['step']);occupied=set(f['occupied'].tolist())
    for ident in args.ids:
        source=args.source/f'cycle-{ident}.npz'
        with np.load(source) as f:params=f['params'];raw_points=f['points'];seed=f['orbit'][0];period=len(f['orbit'])
        assert np.array_equal(params,trap_params) and np.array_equal(raw_points,trap_points)
        beta,tau=[exactfloat(a,mp.mp) for a in params];points=[[exactfloat(a,mp.mp) for a in p] for p in raw_points];x=[exactfloat(a,mp.mp) for a in seed];I=eye(mp.mp)
        for iteration in range(20):
            y,M,orbit,axes=model(x,beta,tau,points,period,mp.mp);res=[y[i]-x[i] for i in range(3)]
            if max(abs(v) for v in res)<mp.mpf('1e-160'):break
            delta=mp.lu_solve(mp.matrix([[M[i][j]-I[i][j] for j in range(3)] for i in range(3)]),mp.matrix(res));x=[x[i]-delta[i] for i in range(3)]
        else:raise ArithmeticError('Newton refinement failed')
        decimal=lambda v:mp.nstr(v,210)
        center=[iv.mpf(decimal(v)) for v in x];C=mp.inverse(mp.matrix([[M[i][j]-I[i][j] for j in range(3)] for i in range(3)]));ci=[[iv.mpf(decimal(C[i,j])) for j in range(3)] for i in range(3)]
        b,t=[exactfloat(a,iv) for a in params];p=[[exactfloat(a,iv) for a in q] for q in raw_points];rad=iv.mpf('1e-100');box=[iv.mpf([(v-rad).a,(v+rad).b]) for v in center];II=eye(iv)
        y0,_,_,_=model(center,b,t,p,period,iv);fx=[y0[i]-center[i] for i in range(3)];_,MI,orbit_iv,axes_iv=model(box,b,t,p,period,iv)
        H=[[MI[i][j]-II[i][j] for j in range(3)] for i in range(3)];prod=matmul(ci,H);R=[[II[i][j]-prod[i][j] for j in range(3)] for i in range(3)];cf=matvec(ci,fx);remainder=matvec(R,[box[i]-center[i] for i in range(3)]);K=[center[i]-cf[i]+remainder[i] for i in range(3)]
        assert all(v.a>w.a and v.b<w.b for v,w in zip(K,box));contraction=max(sum(abs(v).b for v in row) for row in R);assert contraction<1
        exclusions=[]
        for d in range(1,period):
            if period%d:continue
            difference=[orbit_iv[d][j]-box[j] for j in range(3)];assert any(v.a>0 or v.b<0 for v in difference);exclusions.append(dict(divisor=d,difference=[str(v) for v in difference]))
        coefficients=[-sum(MI[i][i] for i in range(3)),sum(MI[i][i]*MI[j][j]-MI[i][j]*MI[j][i] for i in range(3) for j in range(i+1,3)),-det(MI)]
        def polynomial(v):return ((v+coefficients[0])*v+coefficients[1])*v+coefficients[2]
        eigen=sorted(mp.eig(mp.matrix(M),left=False,right=False),key=lambda v:mp.re(v));brackets=[];unstable=0
        for value in eigen:
            assert abs(mp.im(value))<mp.mpf('1e-150');value=mp.re(value)
            low=iv.mpf(decimal(value-mp.mpf('1e-30')));high=iv.mpf(decimal(value+mp.mpf('1e-30')));left=polynomial(low);right=polynomial(high)
            assert (left.b<0 and right.a>0) or (left.a>0 and right.b<0)
            enclosure=iv.mpf([low.a,high.b]);expanding=bool(enclosure.b<-1 or enclosure.a>1);stable=bool(enclosure.a>-1 and enclosure.b<1);assert expanding or stable
            unstable+=int(expanding);brackets.append(dict(interval=str(enclosure),at_left=str(left),at_right=str(right),expanding=expanding))
        assert all(iv.mpf(brackets[i]['interval']).b<iv.mpf(brackets[i+1]['interval']).a for i in range(2))
        membership=[]
        for phase,q in enumerate(orbit_iv):
            intervals=[(q[j]-exactfloat(origin[j],iv))/exactfloat(h,iv) for j in range(3)];lower=[int(mp.floor(mp.mpf(v.a))) for v in intervals];upper=[int(mp.floor(mp.mpf(v.b))) for v in intervals]
            assert all(0<=lower[j]<=upper[j]<shape[j] for j in range(3));keys=[]
            for ix in range(lower[0],upper[0]+1):
                for iy in range(lower[1],upper[1]+1):
                    for iz in range(lower[2],upper[2]+1):
                        key=int((ix*shape[1]+iy)*shape[2]+iz);assert key in occupied;keys.append(key)
            membership.append(dict(phase=phase,box=[str(v) for v in q],intersected_keys=keys))
        row=dict(id=ident,period=period,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),center_decimal=[decimal(v) for v in x],radius='1e-100',box=[str(v) for v in box],krawczyk=[str(v) for v in K],contraction_bound=str(contraction),newton_residual=decimal(max(abs(v) for v in res)),proper_divisor_exclusions=exclusions,monodromy=[[str(v) for v in r] for r in MI],characteristic_coefficients=[str(v) for v in coefficients],multiplier_brackets=brackets,unstable_dimensions=unstable,stable_dimensions=3-unstable,trapping_membership=membership,minimum_radius=str(min(v.a for v in axes_iv)),numerical_multipliers=[mp.nstr(v,60) for v in eigen])
        rows.append(row);(args.out/f'cycle-{ident}.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(dict(id=ident,period=period,unstable=unstable,multipliers=row['numerical_multipliers'])),flush=True)
    sources={}
    for name in [Path(__file__).name,'certify_raar_cycles.py']:
        data=Path(__file__).with_name(name).read_bytes();(args.out/name).write_bytes(data);sources[name]=hashlib.sha256(data).hexdigest()
    result=dict(status='passed',interval_digits=200,mpmath=mp.__version__,trap_sha256=hashlib.sha256(args.trap.read_bytes()).hexdigest(),params=trap_params.tolist(),points=trap_points.tolist(),sources=sources,cases=rows,seconds=time.time()-begin,scope='Located least-period hyperbolic cycles with exactly classified real multipliers and full orbit membership in the certified trapping region. Their periodicity implies membership in its maximal invariant set. No transitivity, heterodimensional connection, entropy or historical uniqueness claim.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--trap',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',default=[0,5,6]);p.add_argument('--mpmath-path',type=Path);main(p.parse_args())
