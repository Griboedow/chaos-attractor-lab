"""Independent interval Krawczyk existence and stability certificates for RAAR cycles."""
import argparse,hashlib,json,time
from pathlib import Path
import mpmath as mp


def eye(ctx):return [[ctx.mpf(int(i==j)) for j in range(3)] for i in range(3)]
def matmul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def matvec(A,x):return [sum(A[i][j]*x[j] for j in range(3)) for i in range(3)]
def det(A):return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])


def projection(x,c,n,r,ctx):
    dx=[x[i]-c[i] for i in range(3)];dot=sum(dx[i]*n[i] for i in range(3));v=[dx[i]-dot*n[i] for i in range(3)]
    length=ctx.sqrt(sum(z**2 for z in v));u=[z/length for z in v]
    P=[c[i]+r*u[i] for i in range(3)]
    J=[[r/length*(int(i==j)-n[i]*n[j]-u[i]*u[j]) for j in range(3)] for i in range(3)]
    return P,J,length


def iterate(x,p,period,ctx):
    beta,r,theta,*center=p;normal=[ctx.mpf(0),-ctx.sin(theta),ctx.cos(theta)];zero=[ctx.mpf(0)]*3;vertical=[ctx.mpf(0),ctx.mpf(0),ctx.mpf(1)]
    I=eye(ctx);M=eye(ctx);orbit=[];axes=[]
    for _ in range(period):
        orbit.append(x);b,jb,lb=projection(x,center,normal,r,ctx);rb=[2*b[i]-x[i] for i in range(3)];a,ja,la=projection(rb,zero,vertical,ctx.mpf(1),ctx)
        ra=[2*a[i]-rb[i] for i in range(3)];new=[beta*(ra[i]+x[i])/2+(1-beta)*b[i] for i in range(3)]
        product=matmul(ja,[[2*jb[i][j]-I[i][j] for j in range(3)] for i in range(3)])
        J=[[beta*I[i][j]+beta*product[i][j]+(1-2*beta)*jb[i][j] for j in range(3)] for i in range(3)]
        M=matmul(J,M);x=new;axes.extend([la,lb])
    return x,M,orbit,axes


def exactfloat(x,ctx):
    n,d=float(x).as_integer_ratio();return ctx.mpf(n)/d


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);mp.mp.dps=200;iv=mp.iv;iv.dps=200;begin=time.time();rows=[]
    source=Path(__file__);(args.out/source.name).write_bytes(source.read_bytes());(args.out/'input.json').write_bytes(args.input.read_bytes())
    for case in json.loads(args.input.read_text())['cases']:
        p=[exactfloat(v,mp.mp) for v in case['params']];x=[mp.mpf(v) for v in case['state']];period=case['period'];I=eye(mp.mp)
        for iteration in range(20):
            y,M,_,_=iterate(x,p,period,mp.mp);res=[y[i]-x[i] for i in range(3)]
            if max(abs(v) for v in res)<mp.mpf('1e-160'):break
            J=mp.matrix([[M[i][j]-I[i][j] for j in range(3)] for i in range(3)]);delta=mp.lu_solve(J,mp.matrix(res));x=[x[i]-delta[i] for i in range(3)]
        else:raise ArithmeticError('Newton failed to refine '+case['label'])
        C=mp.inverse(mp.matrix([[M[i][j]-I[i][j] for j in range(3)] for i in range(3)]))
        decimal=lambda v:mp.nstr(v,210)
        center=[iv.mpf(decimal(v)) for v in x];ci=[[iv.mpf(decimal(C[i,j])) for j in range(3)] for i in range(3)];params=[exactfloat(v,iv) for v in case['params']]
        radius=iv.mpf('1e-120');box=[iv.mpf([(v-radius).a,(v+radius).b]) for v in center];II=eye(iv)
        y0,_,_,_=iterate(center,params,period,iv);fx=[y0[i]-center[i] for i in range(3)]
        y,mi,orbit,axes=iterate(box,params,period,iv);J=[[mi[i][j]-II[i][j] for j in range(3)] for i in range(3)]
        prod=matmul(ci,J);R=[[II[i][j]-prod[i][j] for j in range(3)] for i in range(3)]
        cf=matvec(ci,fx);remaining=matvec(R,[box[i]-center[i] for i in range(3)]);K=[center[i]-cf[i]+remaining[i] for i in range(3)]
        assert all(k.a>b.a and k.b<b.b for k,b in zip(K,box)),case['label']
        contraction=max(sum(abs(v).b for v in row) for row in R);assert contraction<1
        divisors=[]
        for d in range(1,period):
            if period%d:continue
            residual=[orbit[d][j]-box[j] for j in range(3)];assert any(v.a>0 or v.b<0 for v in residual)
            divisors.append(dict(divisor=d,coordinate_difference=[str(v) for v in residual]))
        a1=-sum(mi[i][i] for i in range(3));a2=sum(mi[i][i]*mi[j][j]-mi[i][j]*mi[j][i] for i in range(3) for j in range(i+1,3));a3=-det(mi)
        jury=[1-abs(a3),1+a1+a2+a3,1-a1+a2-a3,1-a2+a1*a3-a3*a3]
        stable=all(v.a>0 for v in jury)
        # A real root below -1 is guaranteed when q(-1)>0 for a monic cubic.
        qm1=-1+a1-a2+a3;qp1=1+a1+a2+a3
        unstable=bool(qm1.a>0 or qp1.b<0 or abs(a3).a>1)
        expected=case['label'].startswith('stable-');assert stable==expected and (stable or unstable),case['label']
        record=dict(label=case['label'],period=period,parameters=case['params'],center_decimal=[decimal(v) for v in x],box_radius='1e-120',box=[str(v) for v in box],krawczyk=[str(v) for v in K],contraction_bound=str(contraction),minimal_period_certified=True,proper_divisor_exclusions=divisors,axis_distances=[str(v) for v in axes],monodromy=[[str(v) for v in row] for row in mi],characteristic_coefficients=[str(v) for v in [a1,a2,a3]],schur_jury_margins=[str(v) for v in jury],characteristic_at_minus_one=str(qm1),characteristic_at_plus_one=str(qp1),stable_certified=stable,unstable_certified=unstable,newton_residual=decimal(max(abs(v) for v in res)))
        (args.out/(case['label']+'.json')).write_text(json.dumps(record,indent=2),encoding='utf-8');rows.append(record);print(json.dumps({k:record[k] for k in ['label','period','minimal_period_certified','stable_certified','unstable_certified','newton_residual']}),flush=True)
    report=dict(status='passed',interval_digits=200,mpmath=mp.__version__,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),seconds=time.time()-begin,cases=rows,scope='Existence, uniqueness within each stated box, least period and local Floquet stability at isolated exact settings. No proof of connecting bifurcations, cascade accumulation, attracting chaos or historical uniqueness.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
