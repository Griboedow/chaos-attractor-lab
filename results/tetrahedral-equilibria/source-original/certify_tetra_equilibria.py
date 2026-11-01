"""Local interval root certificates and Routh counts; no global completeness claim."""
from pathlib import Path
import argparse,json,time
import numpy as np
import mpmath as mp

parser=argparse.ArgumentParser();parser.add_argument('--precision',type=int,default=60);args=parser.parse_args()
base=Path(__file__).parent;d=np.load(base/'tetra-v1/data.npz');inventory=json.loads((base/'tetra-equilibria/inventory.json').read_text());out=base/'tetra-equilibria';mp.iv.dps=args.precision
iv=mp.iv;one=iv.mpf(1);zero=iv.mpf(0);rad=iv.mpf('1e-8')
def exact(a):
    n,q=float(a).as_integer_ratio();return iv.mpf(n)/q
def matmul(a,b):return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def matvec(a,b):return [sum(a[i][k]*b[k] for k in range(3)) for i in range(3)]
def determinant(a):return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
def absupper(a):return max(abs(float(a.a)),abs(float(a.b)))
def sign(a):
    if a.a>0:return 1
    if a.b<0:return -1
    raise ValueError('Interval sign is unresolved')
rows=[];begin=time.time()
for case in inventory:
    ident=case['id'];centers=[[exact(a) for a in c] for c in d['centers'][ident]];axes=[[exact(a) for a in c] for c in d['axes'][ident]];sigma=exact(d['sigma'][ident]);gamma=exact(d['gamma'][ident])
    def evaluate(x):
        f=[-gamma*a for a in x];jac=[[-gamma if i==j else zero for j in range(3)] for i in range(3)]
        for c,u in zip(centers,axes):
            r=[x[j]-c[j] for j in range(3)];cross=[u[1]*r[2]-u[2]*r[1],u[2]*r[0]-u[0]*r[2],u[0]*r[1]-u[1]*r[0]];g=iv.exp(-sum(a**2 for a in r)/(2*sigma**2))
            for j in range(3):
                f[j]+=g*cross[j]
                for k in range(3):jac[j][k]-=g*cross[j]*r[k]/sigma**2
            jac[0][1]-=g*u[2];jac[0][2]+=g*u[1];jac[1][0]+=g*u[2];jac[1][2]-=g*u[0];jac[2][0]-=g*u[1];jac[2][1]+=g*u[0]
        return f,jac
    certified=[]
    for root in case['roots']:
        x=[exact(a) for a in root['point']];box=[iv.mpf([a.a-rad,a.b+rad]) for a in x]
        f,j0=evaluate(x);_,j=evaluate(box);mid=np.array([[float(v.mid) for v in row] for row in j0]);C=[[exact(v) for v in row] for row in np.linalg.inv(mid)]
        cj=matmul(C,j);err=[[(one if a==b else zero)-cj[a][b] for b in range(3)] for a in range(3)];correction=matvec(C,f)
        delta=matvec(err,[box[i]-x[i] for i in range(3)]);K=[x[i]-correction[i]+delta[i] for i in range(3)]
        assert all(K[i].a>box[i].a and K[i].b<box[i].b for i in range(3))
        # The strict norm bound is evaluated as intervals, not rounded floats.
        normrows=[sum(abs(v) for v in row) for row in err];assert all(v.b<1 for v in normrows)
        a1=-sum(j[i][i] for i in range(3));a2=j[0][0]*j[1][1]-j[0][1]*j[1][0]+j[0][0]*j[2][2]-j[0][2]*j[2][0]+j[1][1]*j[2][2]-j[1][2]*j[2][1];a3=-determinant(j)
        first=[one,a1,(a1*a2-a3)/a1,a3];signs=[sign(v) for v in first];positive=sum(signs[i]!=signs[i+1] for i in range(3));assert positive==root['positive_real_count']
        certified.append({'center':root['point'],'radius':'1e-8','Krawczyk_coordinates':[str(a) for a in K],'contraction_row_sums':[str(a) for a in normrows],'Routh_first_column':[str(a) for a in first],'Routh_signs':signs,'positive_real_eigenvalue_count':positive})
    centers_float=np.array([r['center'] for r in certified]);dist=np.linalg.norm(centers_float[:,None]-centers_float[None,:],axis=2);np.fill_diagonal(dist,np.inf);assert dist.min()>1e-6
    rows.append({'id':ident,'certified_distinct_roots':len(certified),'roots':certified,'not_exhaustive':True})
report={'method':f'Outward mpmath interval arithmetic, {args.precision} decimal digits. Exact binary64 parameters converted by integer ratios. Strict Krawczyk inclusion with an interval contraction norm below one certifies one root per disjoint box. Routh first-column signs certify unstable half-plane eigenvalue counts throughout each box. No claim of exhaustive root enumeration or chaos certification.','seconds':time.time()-begin,'cases':rows}
(out/f'certificate-{args.precision}.json').write_text(json.dumps(report,indent=2));print(json.dumps({'seconds':report['seconds'],'cases':[{'id':r['id'],'certified_distinct_roots':r['certified_distinct_roots']} for r in rows]}),flush=True)
