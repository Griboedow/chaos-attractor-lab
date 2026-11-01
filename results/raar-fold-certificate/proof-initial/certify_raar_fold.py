"""Interval certificate of a regular Whitney fold near a finite reference cloud."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import mpmath as mp
from certify_raar_cycles import exactfloat,eye,matmul,det,iterate
from certify_raar_multiplier_events import lift
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def projection2(x,c,n,r,ctx):
    Q=[[ctx.mpf(int(i==j))-n[i]*n[j] for j in range(3)] for i in range(3)];v=[sum(Q[i][j]*(x[j]-c[j]) for j in range(3)) for i in range(3)];L=ctx.sqrt(sum(a*a for a in v));P=[c[i]+r*v[i]/L for i in range(3)]
    J=[[r*(Q[i][j]/L-v[i]*v[j]/L**3) for j in range(3)] for i in range(3)]
    H=[[[r*(-(Q[i][j]*v[k]+Q[i][k]*v[j]+Q[j][k]*v[i])/L**3+3*v[i]*v[j]*v[k]/L**5) for k in range(3)] for j in range(3)] for i in range(3)];return P,J,H,L
def jet(x,p,ctx):
    beta,r,theta,*c=p;I=eye(ctx);zero=[ctx.mpf(0)]*3;n=[ctx.mpf(0),-ctx.sin(theta),ctx.cos(theta)];vertical=[ctx.mpf(0),ctx.mpf(0),ctx.mpf(1)]
    b,B,HB,lb=projection2(x,c,n,r,ctx);a,A,HA,la=projection2([2*b[i]-x[i] for i in range(3)],zero,vertical,ctx.mpf(1),ctx);U=[[2*B[i][j]-I[i][j] for j in range(3)] for i in range(3)];AU=matmul(A,U)
    J=[[beta*I[i][j]+beta*AU[i][j]+(1-2*beta)*B[i][j] for j in range(3)] for i in range(3)]
    H=[[[beta*sum(HA[i][m][n]*U[m][j]*U[n][k] for m in range(3) for n in range(3))+2*beta*sum(A[i][m]*HB[m][j][k] for m in range(3))+(1-2*beta)*HB[i][j][k] for k in range(3)] for j in range(3)] for i in range(3)]
    co=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rr=[q for q in range(3) if q!=i];cc=[q for q in range(3) if q!=j];row.append((-1)**(i+j)*(J[rr[0]][cc[0]]*J[rr[1]][cc[1]]-J[rr[0]][cc[1]]*J[rr[1]][cc[0]]))
        co.append(row)
    grad=[sum(co[i][j]*H[i][j][k] for i in range(3) for j in range(3)) for k in range(3)]
    return det(J),grad,J,H,[la,lb]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();mp.mp.dps=200;mp.iv.dps=200;iv=mp.iv;source=json.loads(args.input.read_text());p=[exactfloat(v,mp.mp) for v in source['params']];a,b=[[exactfloat(v,mp.mp) for v in q] for q in source['segment_endpoints']];delta=[b[i]-a[i] for i in range(3)];t=mp.mpf(source['t_proposal'])
    for _ in range(20):
        x=[a[i]+t*delta[i] for i in range(3)];g,grad,J,H,axes=jet(x,p,mp.mp)
        if abs(g)<mp.mpf('1e-160'):break
        t-=g/sum(grad[i]*delta[i] for i in range(3))
    else:raise ArithmeticError('Newton failed')
    assert 0<t<1
    # Independently differentiate the nested-reflector map (not this Hessian formula).
    derivative_errors=[]
    for i,j,k in itertools.product(range(3),repeat=3):
        order=[0,0,0];order[j]+=1;order[k]+=1
        q=mp.diff(lambda *v:iterate(list(v),p,1,mp.mp)[0][i],tuple(x),tuple(order));derivative_errors.append(abs(q-H[i][j][k]))
    assert max(derivative_errors)<mp.mpf('1e-150')
    center=lift(t,iv);radius=iv.mpf('1e-100');T=iv.mpf([(center-radius).a,(center+radius).b]);ai,bi=[[exactfloat(v,iv) for v in q] for q in source['segment_endpoints']];di=[bi[i]-ai[i] for i in range(3)];pi=[exactfloat(v,iv) for v in source['params']];X=[ai[i]+T*di[i] for i in range(3)];gi,gradient,Ji,Hi,axes=jet(X,pi,iv);derivative=sum(gradient[i]*di[i] for i in range(3));assert derivative.a>0 or derivative.b<0
    gc=jet([ai[i]+center*di[i] for i in range(3)],pi,iv)[0];N=center-gc/derivative;assert T.a<N.a and N.b<T.b
    pairs=list(itertools.combinations(range(3),2));rows=max(pairs,key=lambda ij:sum(v*v for v in cross(J[ij[0]],J[ij[1]])));kernel=cross(Ji[rows[0]],Ji[rows[1]]);component=max(range(3),key=lambda k:abs(cross(J[rows[0]],J[rows[1]])[k]));assert kernel[component].a>0 or kernel[component].b<0
    transversality=sum(gradient[i]*kernel[i] for i in range(3));assert transversality.a>0 or transversality.b<0;assert all(v.a>0 for v in axes)
    image=iterate(X,pi,1,iv)[0];report=dict(status='passed',interval_digits=200,mpmath=mp.__version__,input_sha256=sha(args.input),source_sha256=sha(Path(__file__)),t_center=mp.nstr(t,210),t_box=str(T),interval_newton=str(N),state_box=[str(v) for v in X],image_box=[str(v) for v in image],determinant_box=str(gi),line_derivative=str(derivative),null_rows=list(rows),nonzero_minor_component=component,nonzero_minor=str(kernel[component]),kernel=[str(v) for v in kernel],determinant_gradient=[str(v) for v in gradient],kernel_transversality=str(transversality),axis_distances=[str(v) for v in axes],independent_hessian_max_residual=mp.nstr(max(derivative_errors),210),state_center=[mp.nstr(v,210) for v in x],seconds=time.time()-begin,scope='One certified regular corank-one critical point with kernel transverse to det(DF)=0: the standard fold criterion. Exact binary64 segment endpoints/parameters. Near a finite reference sample, not certified membership in an attracting set, causal explanation of the arch, new singularity type, or historical-priority proof.')
    for name in ['certify_raar_fold.py','certify_raar_cycles.py','certify_raar_multiplier_events.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (args.out/'input.json').write_bytes(args.input.read_bytes());(args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({k:report[k] for k in ['status','t_center','nonzero_minor','kernel_transversality','axis_distances','independent_hessian_max_residual','seconds']}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
