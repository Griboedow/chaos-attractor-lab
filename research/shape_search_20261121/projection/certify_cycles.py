"""Interval existence, least period and unstable index for selected return cycles."""
import argparse
import json
import os
import sys
from pathlib import Path
import mpmath as mp

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from certify_raar_cycles import exactfloat,eye,matmul,matvec,det


class Jet:
    """First derivatives in three initial coordinates, over either scalar context."""
    def __init__(self,v,d):
        self.v=v
        self.d=d

    def __add__(self,b):
        if not isinstance(b,Jet): return Jet(self.v+b,self.d)
        return Jet(self.v+b.v,[x+y for x,y in zip(self.d,b.d)])
    __radd__=__add__

    def __neg__(self):
        return Jet(-self.v,[-x for x in self.d])

    def __sub__(self,b):
        return self+(-b)

    def __rsub__(self,b):
        return -self+b

    def __mul__(self,b):
        if not isinstance(b,Jet): return Jet(self.v*b,[v*b for v in self.d])
        return Jet(self.v*b.v,[x*b.v+self.v*y for x,y in zip(self.d,b.d)])
    __rmul__=__mul__

    def __truediv__(self,b):
        if not isinstance(b,Jet): return Jet(self.v/b,[v/b for v in self.d])
        value=self.v/b.v
        return Jet(value,[(x-value*y)/b.v for x,y in zip(self.d,b.d)])

    def square(self):
        return Jet(self.v**2,[2*self.v*v for v in self.d])

    def sqrt(self,ctx):
        if not self.v>0: raise ArithmeticError('Unresolved projector singularity')
        root=ctx.sqrt(self.v)
        return Jet(root,[v/(2*root) for v in self.d])

    def exp(self,ctx):
        value=ctx.exp(self.v)
        return Jet(value,[value*v for v in self.d])


def relaxation(raw,ctx):
    beta=exactfloat(raw['beta'],ctx)
    if raw.get('beta_radius') and ctx is mp.iv:
        radius=ctx.mpf(raw['beta_radius'])
        return beta+ctx.mpf([-radius.b,radius.b])
    return beta


def iterate(initial,period,rec,ctx):
    raw=rec['parameters']
    beta=relaxation(raw,ctx);inv=exactfloat(raw['inv_tau'],ctx)
    tube=exactfloat(raw['tube_radius'],ctx)
    sites=[[exactfloat(v,ctx) for v in p] for p in raw['sites']]
    mix=exactfloat(raw.get('mix_temperature',0),ctx)
    ident=eye(ctx)
    x=[Jet(v,ident[i]) for i,v in enumerate(initial)]
    orbit=[]
    for _ in range(period):
        orbit.append([v.v for v in x])
        weights=[(-sum((x[i]-p[i]).square() for i in range(3))*inv).exp(ctx) for p in sites]
        total=sum(weights)
        m=[sum(w*p[i] for w,p in zip(weights,sites))/total for i in range(3)]
        v=[u*2-t for u,t in zip(m,x)]
        def surface(second):
            vx,vy,vz=(v[0]-1,v[2],v[1]) if second else tuple(v)
            rho=(vx.square()+vy.square()).sqrt(ctx)
            length=((rho-1).square()+vz.square()).sqrt(ctx)
            radial=(rho-1)/length*tube+1
            px=vx/rho*radial;py=vy/rho*radial;pz=vz/length*tube
            return ([px+1,pz,py] if second else [px,py,pz]),(length-tube).square()
        p,da=surface(False)
        if raw['second_ring']:
            q,db=surface(True)
            if raw.get('mix_temperature',0):
                wa=(-da/mix).exp(ctx);wb=(-db/mix).exp(ctx)
                p=[(a*wa+b*wb)/(wa+wb) for a,b in zip(p,q)]
            elif db.v<da.v:
                p=q
            elif not da.v<db.v:
                raise ArithmeticError('Interval overlaps the nearest-surface switch')
        # Nested reflector form, differentiated with Jets across the entire cycle.
        x=[(a*2-b+c)*(beta/2)+e*(1-beta) for a,b,c,e in zip(p,v,x,m)]
    return [v.v for v in x],[v.d for v in x],orbit


def refine(seed,period,rec):
    ctx=mp.mp;ident=eye(ctx)
    x=[ctx.mpf(v) for v in seed]
    for _ in range(15):
        y,m,_=iterate(x,period,rec,ctx)
        residual=[y[i]-x[i] for i in range(3)]
        if max(abs(v) for v in residual)<ctx.mpf('1e-80'): break
        h=mp.matrix([[m[i][j]-ident[i][j] for j in range(3)] for i in range(3)])
        change=mp.lu_solve(h,mp.matrix(residual))
        x=[x[i]-change[i] for i in range(3)]
    else: raise ArithmeticError('High-precision Newton refinement did not converge')
    h=mp.matrix([[m[i][j]-ident[i][j] for j in range(3)] for i in range(3)])
    c=mp.inverse(h)
    return [mp.nstr(v,100) for v in x],[[mp.nstr(c[i,j],100) for j in range(3)] for i in range(3)]


def verify(certificate,rec,radius='1e-52',iterate_model=iterate):
    iv=mp.iv;ident=eye(iv)
    p=certificate['period']
    x=[iv.mpf(v) for v in certificate['center']]
    c=[[iv.mpf(v) for v in row] for row in certificate['preconditioner']]
    rad=iv.mpf(radius)
    box=[iv.mpf([(v-rad).a,(v+rad).b]) for v in x]
    y,_,_=iterate_model(x,p,rec,iv)
    _,matrix,orbit=iterate_model(box,p,rec,iv)
    residual=[y[i]-x[i] for i in range(3)]
    derivative=[[matrix[i][j]-ident[i][j] for j in range(3)] for i in range(3)]
    product=matmul(c,derivative)
    rem=[[ident[i][j]-product[i][j] for j in range(3)] for i in range(3)]
    cf=matvec(c,residual)
    remainder=matvec(rem,[box[i]-x[i] for i in range(3)])
    image=[x[i]-cf[i]+remainder[i] for i in range(3)]
    if not all(u.a>v.a and u.b<v.b for u,v in zip(image,box)): raise ArithmeticError('Krawczyk inclusion failed')
    contraction=max(sum(abs(v).b for v in row).b for row in rem)
    if not contraction<1: raise ArithmeticError('Contraction failed')
    excluded=[]
    for d in range(1,p):
        if p%d: continue
        if not any((orbit[d][j]-box[j]).a>0 or (orbit[d][j]-box[j]).b<0 for j in range(3)):
            raise ArithmeticError('Possible lower period')
        excluded.append(d)
    # Transform z=(1+s)/(1-s). Positive real parts of s count multipliers |z|>1.
    a=-sum(matrix[i][i] for i in range(3))
    b=sum(matrix[i][i]*matrix[j][j]-matrix[i][j]*matrix[j][i] for i in range(3) for j in range(i+1,3))
    cc=-det(matrix)
    q0,q1,q2,q3=1+a+b+cc,3+a-b-3*cc,3-a-b+3*cc,1-a+b-cc
    if not (q2>0 or q2<0): raise ArithmeticError('Inconclusive Routh pivot')
    column=[q3,q2,(q2*q1-q3*q0)/q2,q0]
    signs=[1 if v>0 else -1 if v<0 else 0 for v in column]
    if 0 in signs: raise ArithmeticError('Inconclusive spectral signs')
    unstable=sum(a!=b for a,b in zip(signs,signs[1:]))
    return dict(radius=radius,unstable_dimension=unstable,routh_signs=signs,excluded_divisors=excluded,contraction_upper_display=float(contraction.b),max_krawczyk_fraction_display=max(float((abs(image[i]-x[i])/rad).b) for i in range(3)))


def main(args):
    mp.mp.dps=mp.iv.dps=110
    args.out.mkdir(parents=True,exist_ok=True)
    for case in args.cases:
        source=args.source/f'case-{case}.json'
        rec=json.loads(source.read_text())
        dest=args.out/f'cycles-{case}.json'
        if args.verify:
            doc=json.loads(dest.read_text())
            for c in doc['cycles']:
                result=verify(c,rec,c['proof']['radius'])
                if result['routh_signs']!=c['proof']['routh_signs']: raise ArithmeticError('Certificate mismatch')
            print('verified',case,len(doc['cycles']),flush=True)
            continue
        if dest.exists(): continue
        seeds=json.loads((args.seeds/f'seeds-{case}.json').read_text())['cycles']
        order=sorted(seeds,key=lambda r:(r['period'],r['max_orbit_distance_to_sample_rms']))
        # First try one orbit per proposed unstable dimension, then other periods.
        preferred=[]
        for index in sorted({r['numerical_unstable'] for r in order}):
            preferred.append(next(r for r in order if r['numerical_unstable']==index))
        accepted=[];pairs=set();failures=[]
        for row in preferred+[r for r in order if r not in preferred]:
            pair=(row['period'],row['numerical_unstable'])
            if pair in pairs and not args.all_seeds: continue
            try:
                center,pre=refine(row['seed'],row['period'],rec)
                cert=dict(period=row['period'],center=center,preconditioner=pre,numerical_orbit_to_cloud_rms=row['max_orbit_distance_to_sample_rms'])
                if 'cohort_indices' in row:
                    cert['numerical_cohort_indices']=row['cohort_indices']
                cert['proof']=verify(cert,rec)
            except (ArithmeticError,ValueError,ZeroDivisionError) as exc:
                failures.append(dict(period=row['period'],reason=str(exc)))
                continue
            pairs.add(pair);accepted.append(cert)
            if len(accepted)>=args.limit: break
        tmp=dest.with_suffix('.tmp.json')
        tmp.write_text(json.dumps(dict(case=case,precision_decimal=110,arithmetic='mpmath interval with forward automatic differentiation',scope='Unique least-period cycles and unstable indices for the parameter set. Sample proximity is numerical and does not establish membership in a chaotic attracting set.',cycles=accepted,failed_attempts=failures),indent=2)+'\n',encoding='utf8')
        os.replace(tmp,dest)
        print(case,[(c['period'],c['proof']['unstable_dimension']) for c in accepted],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--seeds',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cases',type=int,nargs='+',required=True)
    p.add_argument('--limit',type=int,default=4)
    p.add_argument('--verify',action='store_true')
    p.add_argument('--all-seeds',action='store_true')
    main(p.parse_args())
