"""Moving-frame interval flow bounds for positive-area cells on the reduced Case 180844 section."""
import argparse
import json
from time import monotonic
import numpy as np
from flint import arb,arb_mat,arb_series,ctx
from interval_return import binary,field,jets,record
from native_section import RAW,P,packed
from reduced import reduced
from return_tube import graph

def eye():
    return arb_mat([[int(i==j) for j in range(3)] for i in range(3)])

def enclosure(v):
    return [v[i,0] for i in range(3)]

def column(v):
    return arb_mat([[a] for a in v])

def mid_matrix(a):
    return arb_mat([[a[i,j].mid() for j in range(a.ncols())] for i in range(a.nrows())])

def local(m,B,err,h,p,n):
    e=column([arb(0,r) for r in err]);full=enclosure(m+B*e)
    pads=[2*h*w.abs_upper()+arb('1e-30') for w in field(full,p)]
    for _ in range(20):
        tube=[arb(q.mid(),q.rad()+s) for q,s in zip(full,pads)]
        on=field(tube,p)
        span=arb(h/2,h/2)
        image=[v+span*f for v,f in zip(full,on)]
        if all(v.contains_interior(w) for v,w in zip(tube,image)):break
        pads=[max(arb(5)/4*s,2*(v-q).abs_upper()) for s,v,q in zip(pads,image,full)]
    else:raise ArithmeticError('Picard tube failed')
    deriv=[field([arb_series([v,int(i==k)],2) for i,v in enumerate(tube)],p) for k in range(3)]
    A=arb_mat([[deriv[j][i][1] for j in range(3)] for i in range(3)])
    norm=h*max(sum(A[i,j].abs_upper() for j in range(3)) for i in range(3))
    if not norm<1:raise ArithmeticError('Picard contraction failed')
    # The ordered k-fold product of time-dependent Jacobians is enclosed by A^k
    # computed with intervals. Integration over the simplex contributes h^k/k!.
    term=eye();J=eye();q=arb(1);scalar_sum=q
    for k in range(1,6):
        term=term*A*(h/k);J+=term;q=q*norm/k;scalar_sum+=q
    tail=(norm.exp()-scalar_sum).abs_upper()
    J+=arb_mat([[arb(0,tail) for _ in range(3)] for _ in range(3)])
    coef=jets(enclosure(m),p,n);res=jets(tube,p,n+1);center=[]
    for cs,rs in zip(coef,res):
        v=cs[-1]
        for a in reversed(cs[:-1]):v=v*h+a
        center.append(v+rs[-1]*h**(n+1))
    cm=column(center);newm=mid_matrix(cm)
    jbasis=J*B
    q,_=np.linalg.qr(np.array([[float(jbasis[i,j].mid()) for j in range(3)] for i in range(3)]))
    newB=arb_mat([[binary(f) for f in row] for row in q])
    inverse=newB.inv()
    # Inverting the exact proposed binary matrix removes any reliance on floating QR orthogonality.
    newe=inverse*(cm-newm)+(inverse*jbasis)*e
    return newm,newB,[newe[i,0].abs_upper() for i in range(3)],tube,on,float(norm.upper())


def run(val,offset,rx,rz,width,steps,n):
    ctx.prec=144;ctx.cap=n+3
    p=[binary(x) for x in P];p += [p[10].cos(),p[10].sin()]
    z=arb(val);g=[binary(x) for x in graph().coef];v=z-binary(.99);x=g[-1]
    for q in reversed(g[:-1]):x=x*v+q
    x+=arb(offset);r=reduced(np.array([float(x)]),np.array([float(z)]),.00625,packed())[0]
    ctr=arb(str(round(r[0],8)));w=arb(width);left=ctr-w;right=ctr+w;h=left/steps
    phi=8*arb.pi()/7;c=phi.cos();s=phi.sin()
    B=arb_mat([[c,0,-s],[s,0,c],[0,1,0]])
    m=column([x*c,x*s,z]);merr=enclosure(B.inv()*(m-mid_matrix(m)))
    err=[v.abs_upper()+arb(k) for v,k in zip(merr,[rx,rz,'0'])];m=mid_matrix(m)
    maxrad=0.;ml=0.;clean=True
    def bounds(m,B,err):return enclosure(m+B*column([arb(0,q) for q in err]))
    for k in range(steps):
        try:m,B,err,tube,ft,lip=local(m,B,err,h,p,n)
        except ArithmeticError as exc:raise ArithmeticError((k,str(exc),[v.str(9) for v in err])) from exc
        state=bounds(m,B,err);ml=max(ml,lip);maxrad=max(maxrad,max(float(q.rad()) for q in state))
        if (k+1)*h>=arb('5.7'):
            if not (tube[1]<0 or tube[0]<arb('.77') or tube[0]>arb('.85') or
                    tube[2]<arb('.945') or tube[2]>arb('1.04') or (ft[1]>0 and state[1]<0)):clean=False
    before=state
    m,B,err,tube,ft,lip=local(m,B,err,right-left,p,n);after=bounds(m,B,err)
    passed=before[1]<0 and after[1]>0 and ft[1]>0 and clean
    # Time slab tube is rigorous. Bound graph-relative event output by interval evaluation.
    z1=tube[2];x1=tube[0];j=g[-1]
    for a in reversed(g[:-1]):j=j*(z1-binary(.99))+a
    d=x1-j
    inside=z1>arb('.98071') and z1<arb('.99943') and d.abs_upper()<arb('1e-6')
    return dict(status='certified_return' if passed else 'inconclusive',reduced_tube_membership_certified=bool(passed and inside),
       initial_z=val,initial_x_offset=offset,native_x_halfwidth=rx,native_z_halfwidth=rz,steps=steps,taylor_order=n,precision_bits=144,
       first_chart_event=[record(left),record(right)],leading_gate_cleared=bool(clean),negative_before=record(before[1]),
       positive_after=record(after[1]),normal_on_crossing=record(ft[1]),event_x=record(x1),event_z=record(z1),
       graph_relative_event=record(d),max_state_radius=maxrad,max_picard_contraction=max(ml,lip),
       coefficients=P.tolist(),graph_coefficients=graph().coef.tolist())

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--id',required=True);ap.add_argument('--z',default='.99');ap.add_argument('--delta',default='0')
    ap.add_argument('--rx',default='1e-9');ap.add_argument('--rz',default='1e-9');ap.add_argument('--width',default='5e-8')
    ap.add_argument('--steps',type=int,default=1100);ap.add_argument('--order',type=int,default=12)
    a=ap.parse_args();start=monotonic()
    try:r=run(a.z,a.delta,a.rx,a.rz,a.width,a.steps,a.order)
    except Exception as e:r=dict(status='inconclusive',error=repr(e))
    r['seconds']=monotonic()-start
    print(a.id,r.get('status'),r.get('reduced_tube_membership_certified'),r.get('max_state_radius'),r.get('error'),r['seconds'],flush=True)
    (RAW/f'affine-{a.id}.json').write_text(json.dumps(r,indent=2)+'\n')
