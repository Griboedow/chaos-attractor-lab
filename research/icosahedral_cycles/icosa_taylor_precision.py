"""High-precision midpoint Taylor integration; no remainder enclosure claim.

The Gaussian coefficients are expanded recursively. Ball arithmetic bounds local
arithmetic operations, but step-end radii are discarded and no Taylor remainder
is bounded. This is an independent numerical precision/order check, not a proof.
"""
from pathlib import Path
import sys
from cycle_paths import ROOT,RESEARCH,INPUT,STRUCTURE,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import argparse,json,time,math
import numpy as np
from flint import arb,ctx

ap=argparse.ArgumentParser();ap.add_argument('--index',type=int,default=0);ap.add_argument('--order',type=int,default=24);ap.add_argument('--step',type=float,default=.04);ap.add_argument('--digits',type=int,default=80);ap.add_argument('--full',action='store_true');ap.add_argument('--benchmark',action='store_true');args=ap.parse_args()
ctx.dps=args.digits;base=OUTPUT;d=np.load(INPUT/'icosa-v1/input.npz');rows=json.loads((base/'icosa-cycle-refinement/report.json').read_text())['rows'];row=next(r for r in rows if r['index']==args.index)
def exact(a):
    n,q=float(a).as_integer_ratio();return arb(n)/q
c=[[exact(v) for v in a] for a in d['centers'][177]];u=[[exact(v) for v in a] for a in d['axes'][177]];s=exact(d['sigma'][177]);gamma=exact(d['gamma'][177]);den=2*s*s
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
uc=[cross(a,b) for a,b in zip(u,c)];c2=[sum(v*v for v in a) for a in c]
def advance(start,h,order):
    x=[start];q=[[] for _ in c];e=[[] for _ in c];A=[];B=[]
    for n in range(order):
        S=sum(sum(x[k][a]*x[n-k][a] for k in range(n+1)) for a in range(3))
        for j in range(len(c)):
            qj=(-S+2*sum(c[j][a]*x[n][a] for a in range(3))-(c2[j] if n==0 else 0))/den;q[j].append(qj)
            e[j].append(qj.exp() if n==0 else sum(k*q[j][k]*e[j][n-k] for k in range(1,n+1))/n)
        A.append([sum(u[j][a]*e[j][n] for j in range(len(c))) for a in range(3)]);B.append([sum(uc[j][a]*e[j][n] for j in range(len(c))) for a in range(3)])
        terms=[cross(A[k],x[n-k]) for k in range(n+1)]
        x.append([(sum(t[a] for t in terms)-B[n][a]-gamma*x[n][a])/(n+1) for a in range(3)])
    end=[];rounding=0.
    for a in range(3):
        v=x[-1][a]
        for n in range(order-1,-1,-1):v=v*h+x[n][a]
        rounding=max(rounding,float(v.rad()));end.append(v.mid())
    return end,rounding
start=[exact(v) for v in row['point']];T=row['period'] if args.full else row['segment_time']
if args.benchmark:T=.04
n=math.ceil(T/args.step);h=exact(T)/n;x=start;begin=time.time();maxround=0.
for k in range(n):
    x,r=advance(x,h,args.order);maxround=max(maxround,r)
    if (k+1)%100==0:print(json.dumps({'index':args.index,'step':k+1,'of':n,'seconds':time.time()-begin}),flush=True)
G=[[exact(v) for v in a] for a in d['group'][row['rotation']]]
target=start if args.full else [sum(G[a][b]*start[b] for b in range(3)) for a in range(3)]
error=sum((x[a]-target[a])**2 for a in range(3)).sqrt()
report={'index':args.index,'full_period':args.full,'benchmark':args.benchmark,'digits':args.digits,'order':args.order,'nominal_step':args.step,'steps':n,'time':T,'endpoint':[v.str(args.digits+5,radius=False) for v in x],'closure_error':error.str(args.digits,radius=False),'largest_discarded_rounding_radius':maxround,'seconds':time.time()-begin,'limits':'Midpoint Taylor approximation; step-end arithmetic radii are discarded and truncation remainder is not enclosed. This is not validated ODE integration or a periodic-orbit certificate.'}
if not args.benchmark:
    out=base/'icosa-taylor-precision';out.mkdir(exist_ok=True);name=f"{args.index}-{'full' if args.full else 'relative'}-d{args.digits}-o{args.order}-h{args.step}.json";assert not (out/name).exists();(out/name).write_text(json.dumps(report,indent=2))
print(json.dumps(report),flush=True)
