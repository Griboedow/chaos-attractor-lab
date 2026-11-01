"""Independent eight-site posterior/automatic-derivative collision verification.

No production factorization or analytic Jacobian is imported. The interval
arithmetic library is shared; formulas, derivative construction, precision and
trapping partition differ from the primary certificate.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time
import numpy as np


class Dual:
    context = None

    def __init__(self, value, derivative=None):
        self.v = value
        self.d = derivative if derivative is not None else [self.context.mpf(0)]*3

    def lift(self, value):
        return value if isinstance(value, Dual) else Dual(self.context.mpf(value))

    def __add__(self, other):
        other = self.lift(other)
        return Dual(self.v+other.v, [a+b for a,b in zip(self.d,other.d)])
    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, [-v for v in self.d])

    def __sub__(self, other): return self+-self.lift(other)
    def __rsub__(self, other): return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        return Dual(self.v*other.v, [a*other.v+self.v*b for a,b in zip(self.d,other.d)])
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.lift(other)
        return Dual(self.v/other.v, [(a*other.v-self.v*b)/(other.v*other.v) for a,b in zip(self.d,other.d)])
    def __rtruediv__(self, other): return self.lift(other)/self

    def exp(self):
        w = self.context.exp(self.v)
        return Dual(w, [w*d for d in self.d])

    def sqrt(self):
        w = self.context.sqrt(self.v)
        return Dual(w, [d/(2*w) for d in self.d])


def posterior_reflections(x, points, beta, tau, ctx, derivatives=False):
    if derivatives:
        Dual.context = ctx
        x = [Dual(v,[ctx.mpf(i==j) for j in range(3)]) for i,v in enumerate(x)]
    weights = []
    for p in points:
        z = -sum((v-b)*(v-b) for v,b in zip(x,p))/tau
        weights.append(z.exp() if derivatives else ctx.exp(z))
    total = sum(weights)
    mean = [sum(w*p[i] for w,p in zip(weights,points))/total for i in range(3)]
    rb = [2*m-v for m,v in zip(mean,x)]
    squared = sum(v*v for v in rb)
    radius = squared.sqrt() if derivatives else ctx.sqrt(squared)
    ra = [2*v/radius-v for v in rb]
    result = [(v+w)*(beta/2)+m*(1-beta) for v,w,m in zip(x,ra,mean)]
    if derivatives: return [z.v for z in result], [z.d for z in result], radius.v
    return result


def matmul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def matvec(a,b): return [sum(row[j]*b[j] for j in range(len(b))) for row in a]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps = mp.iv.dps = 300
    iv = mp.iv
    start = time.time()
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    primary = json.loads(args.certificate.read_text(encoding='utf-8'))
    assert primary['status']=='passed' and primary['source_sha256']==sha(args.source)
    with np.load(args.source) as f:
        raw_points=f['points'][653,:f['counts'][653]].copy()
        raw_params=f['params'][653].copy()
    endpoints=np.unique(raw_points)
    assert len(endpoints)==2 and sorted(map(tuple,raw_points))==sorted(itertools.product(endpoints,repeat=3))
    def exact(v,ctx):
        p,q=float(v).as_integer_ratio()
        return ctx.mpf(p)/ctx.mpf(q)
    points_i=[[exact(v,iv) for v in p] for p in raw_points]
    points_m=[[exact(v,mp.mp) for v in p] for p in raw_points]
    bi,ti=[exact(v,iv) for v in raw_params]
    bm,tm=[exact(v,mp.mp) for v in raw_params]
    def mi(x,derivatives=True): return posterior_reflections(x,points_i,bi,ti,iv,derivatives)
    def mn(x): return posterior_reflections(x,points_m,bm,tm,mp.mp)
    def point(v):
        z=iv.mpf(v)
        assert z._mpi_[0]==v._mpf_==z._mpi_[1]
        return z
    def upper(v): return iv.mpf(v.b)
    def lower(v): return iv.mpf(v.a)
    def magnitude(v): return upper(abs(v))
    def sym(v): return iv.mpf([-v.b,v.b])
    def hull(a,b): return iv.mpf([a.a,b.b])
    def rowsum(a): return max(sum(magnitude(v) for v in row).b for row in a)
    lo,hi=exact(-1.65,iv),exact(-.08,iv)
    bl,bh=[exact(v,iv) for v in endpoints]
    # Two independent scalar Gaussian weights; monotonicity follows from the
    # positive posterior variance. Use twice the primary partition resolution.
    def scalar_mean(z):
        wl=iv.exp(-(z-bl)*(z-bl)/ti)
        wh=iv.exp(-(z-bh)*(z-bh)/ti)
        return (bl*wl+bh*wh)/(wl+wh)
    cells=[]
    for k in range(8192):
        z=hull(lo+(hi-lo)*k/8192,lo+(hi-lo)*(k+1)/8192)
        g=hull(scalar_mean(lower(z)),scalar_mean(upper(z)))
        q=z-2*g
        assert q.a>0
        cells.append((z,g,q))
    qmin=iv.mpf(min(q.a for _,_,q in cells));qmax=iv.mpf(max(q.b for _,_,q in cells))
    images=[]
    for z,g,q in cells:
        ql,qh=lower(q),upper(q)
        ul=-qh/iv.sqrt(qh*qh+2*qmin*qmin)
        uh=-ql/iv.sqrt(ql*ql+2*qmax*qmax)
        images.append(hull(bi*lower(z)+bi*ul+(1-2*bi)*upper(g),bi*upper(z)+bi*uh+(1-2*bi)*lower(g)))
    image=hull(iv.mpf(min(v.a for v in images)),iv.mpf(max(v.b for v in images)))
    assert image.a>lo.b and image.b<hi.a
    print(json.dumps(dict(stage='independent trapping',passed=True)),flush=True)
    def residual(s,t):
        y=mn([s,s,t]);return y[0]-s,y[2]-t
    s,t=mp.findroot(residual,tuple(mp.mpf(v) for v in primary['fixed_point']['center']),tol=mp.mpf('1e-280'))
    center=[point(s),point(t)];rad=iv.mpf('2e-220')
    boxes=[hull(v-rad,v+rad) for v in center];P=[boxes[0],boxes[0],boxes[1]]
    y,Jc,_=mi([center[0],center[0],center[1]])
    _,J,_=mi(P)
    H=[[J[0][0]+J[0][1]-1,J[0][2]],[J[2][0]+J[2][1],J[2][2]-1]]
    Hc=[[Jc[0][0]+Jc[0][1]-1,Jc[0][2]],[Jc[2][0]+Jc[2][1],Jc[2][2]-1]]
    midpoint=lambda z:(mp.mpf(z.a)+mp.mpf(z.b))/2
    C=mp.inverse(mp.matrix([[midpoint(v) for v in row] for row in Hc]))
    Ci=[[point(C[i,j]) for j in range(2)] for i in range(2)]
    CH=matmul(Ci,H);R=[[iv.mpf(i==j)-CH[i][j] for j in range(2)] for i in range(2)]
    correction=matvec(Ci,[y[0]-center[0],y[2]-center[1]])
    rest=matvec(R,[boxes[i]-center[i] for i in range(2)])
    K=[center[i]-correction[i]+rest[i] for i in range(2)]
    assert all(b.a<k.a and k.b<b.b for b,k in zip(boxes,K)) and rowsum(R)<1
    # Include the entire primary root enclosure in the independent uniqueness
    # box. Binary endpoint tuples avoid interpreting rounded display strings.
    old=[iv.mpf([mp.mpf(tuple(v[0])),mp.mpf(tuple(v[1]))]) for v in primary['fixed_point']['box_binary']]
    assert all(b.a<o.a and o.b<b.b for b,o in zip(boxes,old))
    Ef=[[Fraction(v) for v in row] for row in primary['unstable_graph']['E_fraction']]
    If=[[Fraction(v) for v in row] for row in primary['unstable_graph']['inverse_fraction']]
    assert matmul(Ef,If)==[[Fraction(i==j) for j in range(3)] for i in range(3)]
    assert Ef[0][0]==1 and Ef[1][0]==-1 and Ef[2][0]==0 and Ef[0][1:]==Ef[1][1:]
    convert=lambda a:[[iv.mpf(v.numerator)/iv.mpf(v.denominator) for v in row] for row in a]
    E,I=convert(Ef),convert(If)
    radius=iv.mpf(primary['unstable_graph']['radius']);L=iv.mpf(primary['unstable_graph']['lipschitz'])
    offset=matvec(E,[sym(radius),sym(radius*L),sym(radius*L)])
    physical=[P[i]+offset[i] for i in range(3)]
    assert all(v.a>lo.b and v.b<hi.a for v in physical)
    _,J,_=mi(physical);G=matmul(matmul(I,J),E)
    a=lower(G[0][0]);b=iv.mpf(rowsum([G[0][1:]]))
    c=iv.mpf(max(magnitude(G[i][0]).b for i in [1,2]));d=iv.mpf(rowsum([G[i][1:] for i in [1,2]]))
    expansion=a-b*L;newL=(c+d*L)/expansion;contraction=d+(c+d*L)*b/expansion
    assert expansion.a>1 and newL.b<L.a and contraction.b<1 and (c+d*L).b<L.a
    print(json.dumps(dict(stage='independent graph',passed=True)),flush=True)
    eps=iv.mpf(primary['endpoint']['epsilon']);em=mp.mpf(primary['endpoint']['epsilon'])
    assert 0<eps.a and eps.b<radius.a
    xc=[s+em,s-em,t];offset=matvec(E,[eps,sym(L*eps),sym(L*eps)])
    initial=[P[i]+offset[i]-point(xc[i]) for i in range(3)]
    delta=[sym(magnitude(v)) for v in matvec(I,initial)];trail=[]
    for step in range(1,primary['endpoint']['iteration']+1):
        nominal=[point(v) for v in xc];off=matvec(E,delta)
        box=[nominal[i]+off[i] for i in range(3)]
        assert all(lo.b<v.a and v.b<hi.a for v in box)
        _,J,r=mi(box);assert r.a>0
        fi=mi(nominal,False);yc=mn(xc)
        err=matvec(I,[fi[i]-point(yc[i]) for i in range(3)])
        adv=matvec(matmul(matmul(I,J),E),delta)
        delta=[sym(magnitude(err[i]+adv[i])) for i in range(3)];xc=yc
        off=matvec(E,delta);out=[point(xc[i])+off[i] for i in range(3)]
        assert all(lo.b<v.a and v.b<hi.a for v in out)
        diff=out[0]-out[1]
        trail.append(dict(iteration=step,error=mp.nstr(max(mp.mpf(abs(v).b) for v in delta),20),sign=1 if diff.a>0 else -1 if diff.b<0 else 0))
        if step%100==0:print(json.dumps(dict(stage='independent propagation',**trail[-1])),flush=True)
    assert diff.b<0
    report=dict(status='passed',source_sha256=sha(args.source),certificate_sha256=sha(args.certificate),precision_digits=300,mpmath_version=mp.__version__,seconds=time.time()-start,
        trapping=dict(cells=8192,image=str(image),minimum_radius=str(iv.sqrt(3)*qmin)),
        fixed_point=dict(box=[str(v) for v in boxes],krawczyk=[str(v) for v in K],contraction=str(rowsum(R))),
        graph=dict(a=str(a),b=str(b),c=str(c),d=str(d),expansion=str(expansion),new_lipschitz=str(newL),contraction=str(contraction)),
        endpoint=dict(epsilon=primary['endpoint']['epsilon'],iterations=step,box=[str(v) for v in out],x_minus_y=str(diff),trail=trail),
        scope='Independent formulas, automatic differentiation, 300-digit interval precision and 8192-cell trap partition. Uses the same mpmath interval arithmetic library. Revalidates all inequalities needed for the intrinsic collision argument; does not establish equality of K and the sampled portrait, chaos, entropy or historical uniqueness.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status='passed',seconds=report['seconds'],difference=str(diff))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','certificate','out','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
