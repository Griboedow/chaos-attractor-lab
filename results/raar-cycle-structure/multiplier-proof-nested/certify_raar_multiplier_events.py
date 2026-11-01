"""Validate parameter brackets containing a -1 Floquet multiplier event.

Uniform parametric Krawczyk inclusion gives a continuous least-period branch.
Opposite characteristic signs at the bracket endpoints force a -1 event.
No nondegenerate flip normal form or uniqueness of the event is asserted.
"""
import argparse,hashlib,json,time
from pathlib import Path
import mpmath as mp
from certify_raar_cycles import iterate,eye,matmul,matvec,det,exactfloat


def lift(x,ctx):
    sign,mantissa,exponent,_=x._mpf_
    return ctx.mpf((-1 if sign else 1)*mantissa)*ctx.mpf(2)**exponent


def characteristic(M):
    a1=-sum(M[i][i] for i in range(3));a2=sum(M[i][i]*M[j][j]-M[i][j]*M[j][i] for i in range(3) for j in range(i+1,3));a3=-det(M)
    return [a1,a2,a3],-1+a1-a2+a3,3-2*a1+a2


def enclosure(center,point_params,interval_params,period,radius):
    iv=mp.iv;I=eye(iv);y,M,_,_=iterate(center,point_params,period,mp.mp)
    C=mp.inverse(mp.matrix([[M[i][j]-int(i==j) for j in range(3)] for i in range(3)]));ci=[[lift(C[i,j],iv) for j in range(3)] for i in range(3)]
    c=[lift(v,iv) for v in center];radius=iv.mpf(radius);box=[iv.mpf([(v-radius).a,(v+radius).b]) for v in c]
    _,mi,orbit,axes=iterate(box,interval_params,period,iv)
    J=[[mi[i][j]-I[i][j] for j in range(3)] for i in range(3)];product=matmul(ci,J);R=[[I[i][j]-product[i][j] for j in range(3)] for i in range(3)]
    yc,_,_,_=iterate(c,interval_params,period,iv);res=[yc[i]-c[i] for i in range(3)];cf=matvec(ci,res);tail=matvec(R,[box[i]-c[i] for i in range(3)]);K=[c[i]-cf[i]+tail[i] for i in range(3)]
    assert all(k.a>b.a and k.b<b.b for k,b in zip(K,box))
    contraction=max(sum(abs(v).b for v in row) for row in R);assert contraction<1
    for d in range(1,period):
        if period%d==0:assert any((orbit[d][i]-box[i]).a>0 or (orbit[d][i]-box[i]).b<0 for i in range(3))
    assert all(v.a>0 for v in axes)
    coefficients,qminus,derivative=characteristic(mi)
    return dict(box=[str(v) for v in box],krawczyk=[str(v) for v in K],contraction=str(contraction),q_minus_one=str(qminus),q_derivative_at_minus_one=str(derivative),characteristic_coefficients=[str(v) for v in coefficients],least_period_certified=True),qminus,derivative


def refine_cycle(x,p,period):
    for _ in range(20):
        y,M,_,_=iterate(x,p,period,mp.mp);res=mp.matrix([y[i]-x[i] for i in range(3)])
        if max(abs(v) for v in res)<mp.mpf('1e-500'):return x
        J=mp.matrix([[M[i][j]-int(i==j) for j in range(3)] for i in range(3)]);d=mp.lu_solve(J,res);x=[x[i]-d[i] for i in range(3)]
    raise ArithmeticError('Endpoint cycle refinement failed')


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);mp.mp.dps=600;mp.iv.dps=600;iv=mp.iv;begin=time.time();rows=[]
    for name in ['certify_raar_multiplier_events.py','certify_raar_cycles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (args.out/'input.json').write_bytes(args.input.read_bytes())
    for case in json.loads(args.input.read_text())['cases']:
        period=case['period'];p=[exactfloat(v,mp.mp) for v in case['params']]
        def event(x,y,z,beta):
            pp=[beta]+p[1:];state=[x,y,z];new,M,_,_=iterate(state,pp,period,mp.mp)
            return tuple([new[i]-state[i] for i in range(3)]+[det([[M[i][j]+int(i==j) for j in range(3)] for i in range(3)])])
        seed=tuple(mp.mpf(v) for v in case['state']+[case['beta']]);sol=mp.findroot(event,seed,tol=mp.mpf('1e-500'),maxsteps=20,verify=True)
        x=list(sol[:3]);beta=sol[3];p=[beta]+p[1:];delta=mp.mpf('1e-300');left=beta-delta;right=beta+delta
        pi=[iv.mpf([lift(left,iv).a,lift(right,iv).b])]+[exactfloat(v,iv) for v in case['params'][1:]]
        uniform,_,derivative=enclosure(x,p,pi,period,'1e-180')
        assert derivative.a>0 or derivative.b<0
        endpoints=[];signs=[]
        for endpoint in [left,right]:
            pp=[endpoint]+p[1:];xx=refine_cycle(x.copy(),pp,period);piv=[lift(endpoint,iv)]+pi[1:]
            # Tie each endpoint root to the unique uniform branch, rather than
            # merely proving existence of an unrelated periodic point there.
            assert all(abs(lift(xx[i],iv)-lift(x[i],iv)).b+iv.mpf('1e-420') < iv.mpf('0.5e-180') for i in range(3))
            record,qminus,_=enclosure(xx,pp,piv,period,'1e-420')
            assert qminus.a>0 or qminus.b<0;signs.append(1 if qminus.a>0 else -1)
            endpoints.append(dict(beta_decimal=mp.nstr(endpoint,610),center_decimal=[mp.nstr(v,610) for v in xx],inside_uniform_branch_box=True,**record))
        assert signs[0]*signs[1]<0
        record=dict(period=period,base_parameters=case['params'],beta_display=mp.nstr(beta,30),beta_bracket=[mp.nstr(left,610),mp.nstr(right,610)],beta_bracket_width=mp.nstr(right-left,30),center_decimal=[mp.nstr(v,610) for v in x],uniform_branch_certificate=uniform,endpoints=endpoints,endpoint_characteristic_signs=signs,simple_multiplier_certified=True,minus_one_event_exists=True)
        (args.out/f'period-{period}.json').write_text(json.dumps(record,indent=2),encoding='utf-8');rows.append(record);print(json.dumps(dict(period=period,beta_display=record['beta_display'],bracket_width=record['beta_bracket_width'],endpoint_signs=signs)),flush=True)
    report=dict(status='passed',interval_digits=600,mpmath=mp.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),seconds=time.time()-begin,cases=rows,scope='For each bracket, a continuous least-period cycle branch exists and its characteristic polynomial changes sign at -1, forcing an event with a simple -1 multiplier. Neither unique event, nondegenerate flip normal form, connecting doubled branch, infinite cascade, asymptotic chaos nor historical novelty is proved.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
