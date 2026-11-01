"""Interval covering checks for exact binary64 windows of the spatial RAAR map.

The output concerns a compact symbolic invariant set, not attracting chaos or
historical uniqueness. An incomplete run is retained as incomplete.
"""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import mpmath as mp
from certify_raar_cycles import exactfloat,projection,matmul,matvec,det,eye

def inverse(A):
    denominator=det(A);assert denominator.a>0 or denominator.b<0
    result=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rr=[k for k in range(3) if k!=j];cc=[k for k in range(3) if k!=i]
            row.append((-1)**(i+j)*(A[rr[0]][cc[0]]*A[rr[1]][cc[1]]-A[rr[0]][cc[1]]*A[rr[1]][cc[0]])/denominator)
        result.append(row)
    return result

class Model:
    def __init__(self,p):
        iv=mp.iv;self.iv=iv;self.p=[exactfloat(v,iv) for v in p];self.beta,self.r,self.theta,*self.c=self.p;self.n=[iv.mpf(0),-iv.sin(self.theta),iv.cos(self.theta)];self.zero=[iv.mpf(0)]*3;self.vertical=[iv.mpf(0),iv.mpf(0),iv.mpf(1)];self.I=eye(iv)
    def step(self,x,jacobian=True):
        iv=self.iv;b,jb,lb=projection(x,self.c,self.n,self.r,iv);rb=[2*b[i]-x[i] for i in range(3)];a,ja,la=projection(rb,self.zero,self.vertical,iv.mpf(1),iv)
        if not(la.a>0 and lb.a>0):raise ArithmeticError('Projection-axis enclosure touches zero')
        y=[self.beta*x[i]+self.beta*a[i]+(1-2*self.beta)*b[i] for i in range(3)]
        if not jacobian:return y,None,min(la.a,lb.a)
        product=matmul(ja,[[2*jb[i][j]-self.I[i][j] for j in range(3)] for i in range(3)]);J=[[self.beta*self.I[i][j]+self.beta*product[i][j]+(1-2*self.beta)*jb[i][j] for j in range(3)] for i in range(3)]
        return y,J,min(la.a,lb.a)

def enclose(model,source,target,mid,rad):
    iv=mp.iv;S=source['matrix'];center=source['center'];q=[exactfloat(v,iv) for v in mid];delta=[iv.mpf([-exactfloat(v,iv).b,exactfloat(v,iv).b]) for v in rad];offset=matvec(S,q);C=[center[i]+offset[i] for i in range(3)];T=S;offset=matvec(T,delta);X=[C[i]+offset[i] for i in range(3)];minimum=iv.mpf('inf')
    for _ in range(7):
        _,J,axis=model.step(X);minimum=min(minimum,axis);C,_,_=model.step(C,False);T=matmul(J,T);offset=matvec(T,delta);X=[C[i]+offset[i] for i in range(3)]
    transformed=matvec(target['inverse'],[C[i]-target['center'][i] for i in range(3)]);D=matmul(target['inverse'],T);offset=matvec(D,delta);Y=[transformed[i]+offset[i] for i in range(3)]
    return Y,minimum

def successful(Y,kind,orientation):
    if kind=='volume':return bool(Y[0].b<-1 or Y[0].a>1 or all(v.a>-1 and v.b<1 for v in Y[1:]))
    signed=(-orientation if kind=='left' else orientation)*Y[0]
    return bool(signed.a>1)

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);mp.mp.dps=args.digits;mp.iv.dps=args.digits;iv=mp.iv;raw=json.loads(args.input.read_text());model=Model(raw['params']);windows=[];begin=time.time()
    for c,B,r in zip(raw['centers'],raw['bases'],raw['radii']):
        center=[exactfloat(v,iv) for v in c];matrix=[[exactfloat(B[i][j],iv)*exactfloat(r[j],iv) for j in range(3)] for i in range(3)];windows.append(dict(center=center,matrix=matrix,inverse=inverse(matrix)))
    transformed=matvec(windows[0]['inverse'],[windows[1]['center'][i]-windows[0]['center'][i] for i in range(3)]);D=matmul(windows[0]['inverse'],windows[1]['matrix']);offset=matvec(D,[iv.mpf([-1,1])]*3);separation=transformed[0]+offset[0];assert separation.b<-1 or separation.a>1
    records=[];edge_indices=range(len(raw['edges'])) if args.edge is None else [args.edge]
    for edgeid in edge_indices:
        i,j=raw['edges'][edgeid];orientation=int(raw['orientation'][edgeid]);source,target=windows[i],windows[j];scales=[sum(abs(source['matrix'][k][l]).b for k in range(3)) for l in range(3)]
        for kind in ['volume','left','right']:
            mid=[0.,0.,0.];rad=[1.,1.,1.]
            if kind!='volume':mid[0]=-1. if kind=='left' else 1.;rad[0]=0.
            stack=[(mid,rad,'')];accepted=0;examined=0;maxdepth=0;last=time.time();certificate=args.out/f'edge-{edgeid}-{kind}.jsonl.gz';minimum_axis=iv.mpf('inf')
            with gzip.open(certificate,'wt',encoding='utf-8') as f:
                while stack and examined<args.max_cells:
                    mid,rad,path=stack.pop();examined+=1;maxdepth=max(maxdepth,len(path)//2)
                    try:
                        Y,gap=enclose(model,source,target,mid,rad);good=successful(Y,kind,orientation)
                    except (ArithmeticError,ZeroDivisionError,ValueError):good=False;Y=None
                    if good:
                        accepted+=1;minimum_axis=min(minimum_axis,gap);f.write(json.dumps(dict(path=path,image=[str(v) for v in Y],minimum_axis=str(gap)),separators=(',',':'))+'\n')
                    else:
                        split=max(range(3),key=lambda k:float(scales[k])*rad[k]);assert rad[split]>1e-16,(edgeid,kind,path)
                        childrad=rad.copy();childrad[split]*=.5
                        for sign in [-1,1]:
                            childmid=mid.copy();childmid[split]+=sign*childrad[split];stack.append((childmid,childrad,path+str(split)+('0' if sign<0 else '1')))
                    if time.time()-last>15:
                        print(json.dumps(dict(edge=edgeid,kind=kind,examined=examined,accepted=accepted,pending=len(stack),depth=maxdepth,seconds=time.time()-begin)),flush=True);last=time.time()
            row=dict(edge=edgeid,source=i,target=j,kind=kind,status='passed' if not stack else 'incomplete',examined=examined,accepted=accepted,pending=len(stack),maximum_depth=maxdepth,minimum_axis=str(minimum_axis),certificate_sha256=hashlib.sha256(certificate.read_bytes()).hexdigest());records.append(row);print(json.dumps(row),flush=True)
            if stack:(args.out/f'edge-{edgeid}-{kind}-pending.json').write_text(json.dumps(stack),encoding='utf-8')
            (args.out/'progress.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    report=dict(status='passed' if all(r['status']=='passed' for r in records) else 'incomplete',edges=list(edge_indices),checks=records,interval_digits=args.digits,mpmath=mp.__version__,input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),window_0_1_separation=str(separation),seconds=time.time()-begin,scope='Strict interval covering conditions for the listed G=F^7 window edges. Full branching graph requires all edges; attracting chaos, manifold smoothness and historical novelty are not certified.')
    (args.out/'input.json').write_bytes(args.input.read_bytes());(args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for name in ['certify_raar_coverings.py','certify_raar_cycles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--digits',type=int,default=70);p.add_argument('--edge',type=int);p.add_argument('--max-cells',type=int,default=50000);main(p.parse_args())
