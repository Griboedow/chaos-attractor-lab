"""Independent full-eight-site automatic-derivative verification for G=F^4.

Revisits every primary leaf in reverse order and refines where necessary.
Shares mpmath interval arithmetic; does not use the Cartesian factorization
or its analytic Jacobian. Exact windows and the mean-value principle are shared.
"""
import argparse,gzip,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def exact(x):
    f=Fraction(x);return mp.iv.mpf(f.numerator)/f.denominator

class Dual:
    def __init__(self,v,d=None):self.v=v;self.d=[mp.iv.mpf(0)]*3 if d is None else d
    def __add__(self,b):
        if not isinstance(b,Dual):b=Dual(b)
        return Dual(self.v+b.v,[a+c for a,c in zip(self.d,b.d)])
    __radd__=__add__
    def __neg__(self):return Dual(-self.v,[-a for a in self.d])
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return b+-self
    def __mul__(self,b):
        if not isinstance(b,Dual):b=Dual(b)
        return Dual(self.v*b.v,[a*b.v+self.v*c for a,c in zip(self.d,b.d)])
    __rmul__=__mul__
    def __truediv__(self,b):
        if not isinstance(b,Dual):b=Dual(b)
        y=self.v/b.v;return Dual(y,[(a-y*c)/b.v for a,c in zip(self.d,b.d)])
    def __rtruediv__(self,b):return Dual(b)/self
    def square(self):return Dual(self.v**2,[2*self.v*a for a in self.d])
    def sqrt(self):
        y=mp.iv.sqrt(self.v);return Dual(y,[a/(2*y) for a in self.d])
    def exp(self):
        y=mp.iv.exp(self.v);return Dual(y,[y*a for a in self.d])

def square(x):return x.square() if isinstance(x,Dual) else x**2
def sqrt(x):return x.sqrt() if isinstance(x,Dual) else mp.iv.sqrt(x)
def exp(x):return x.exp() if isinstance(x,Dual) else mp.iv.exp(x)
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def mv(A,x):return [sum(A[i][k]*x[k] for k in range(3)) for i in range(3)]

class Model:
    def __init__(self,raw):self.beta,self.tau=map(exact,raw['params']);self.points=[[exact(x) for x in p] for p in raw['points']]
    def field(self,x):
        weights=[exp(-sum(square(x[j]-p[j]) for j in range(3))/self.tau) for p in self.points];total=sum(weights)
        # Subtract one fixed site inside the weighted sum. This exact identity
        # reduces dependency width when a posterior weight nearly equals one.
        anchor=max(range(len(weights)),key=lambda k:(weights[k].v if isinstance(weights[k],Dual) else weights[k]).b)
        m=[sum(weights[k]*(p[j]-self.points[anchor][j]) for k,p in enumerate(self.points))/total+self.points[anchor][j] for j in range(3)]
        reflected=[2*m[j]-x[j] for j in range(3)];r=sqrt(sum(square(v) for v in reflected));axis=r.v if isinstance(r,Dual) else r
        if not axis.a>0:raise ArithmeticError('Sphere normalization contains zero')
        projection=[v/r for v in reflected];second=[2*projection[j]-reflected[j] for j in range(3)]
        return [(second[j]+x[j])*self.beta/2+m[j]*(1-self.beta) for j in range(3)],axis
    def value_jacobian(self,box):
        y,axis=self.field([Dual(v,[mp.iv.mpf(int(i==j)) for j in range(3)]) for i,v in enumerate(box)]);return [v.v for v in y],[v.d for v in y],axis


    def return_value(self, x):
        for _ in range(4): x, _ = self.field(x)
        return x
    def return_jacobian(self, box):
        total = [[mp.iv.mpf(i == j) for j in range(3)] for i in range(3)]
        minimum = mp.iv.mpf('inf')
        for _ in range(4):
            _, jac, axis = self.value_jacobian(box)
            minimum = min(minimum, axis.a)
            center = [mp.iv.mpf((mp.mpf(v.a)+mp.mpf(v.b))/2) for v in box]
            delta = [v-c for v,c in zip(box,center)]
            value, _ = self.field(center)
            spread = mv(jac,delta)
            box = [a+b for a,b in zip(value,spread)]
            total = mm(jac,total)
        return total, minimum

def cell(address,kind):
    mid=[Fraction(0)]*3;rad=[Fraction(1)]*3
    if kind!='volume':mid[0]=Fraction(-1 if kind=='left' else 1);rad[0]=Fraction(0)
    for k in range(0,len(address),2):
        axis=int(address[k]);sign=-1 if address[k+1]=='0' else 1;rad[axis]/=2;mid[axis]+=sign*rad[axis]
    return mid,rad

def good(Y,kind,orientation):
    if kind=='volume':return Y[0].b<-1 or Y[0].a>1 or all(v.a>-1 and v.b<1 for v in Y[1:])
    sign=-orientation if kind=='left' else orientation;return (Y[0].a if sign==1 else -Y[0].b)>1

def enclosure(model,ci,S,cj,T,mid,rad):
    q=list(map(exact,mid));delta=[mp.iv.mpf([-exact(r).b,exact(r).b]) for r in rad];offset=mv(S,q);center=[a+b for a,b in zip(ci,offset)];spread=mv(S,delta);box=[a+b for a,b in zip(center,spread)];J,axis=model.return_jacobian(box);f=model.return_value(center);base=mv(T,[a-b for a,b in zip(f,cj)]);change=mv(mm(mm(T,J),S),delta);return [a+b for a,b in zip(base,change)],axis

def main(args):
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    global mp
    import mpmath as mp
    mp.mp.dps=args.digits;mp.iv.dps=args.digits;primary=json.loads((args.source/'report.json').read_text());assert primary['status']=='passed';args.out.mkdir(parents=True,exist_ok=False);begin=time.time();raw=json.loads((args.source/'input.json').read_text());model=Model(raw);centers=[[exact(v) for v in row] for row in raw['centers']];S=[[[exact(v) for v in row] for row in w['matrix']] for w in raw['exact_windows']];T=[[[exact(v) for v in row] for row in w['inverse']] for w in raw['exact_windows']]
    (args.out/'input.json').write_bytes((args.source/'input.json').read_bytes());(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());checks=[]
    for e,(i,j) in enumerate(raw['edges']):
        for kind in ['volume','left','right']:
            original=args.source/f'edge-{e}-{kind}.jsonl.gz'
            with gzip.open(original,'rt',encoding='utf-8') as handle:paths=[json.loads(line)['path'] for line in handle]
            stack=list(paths);examined=accepted=refinements=0;depth=0;minimum=mp.iv.mpf('inf');outpath=args.out/original.name;last=time.time()
            with gzip.open(outpath,'wt',encoding='utf-8',newline='\n') as handle:
                while stack:
                    address=stack.pop();mid,rad=cell(address,kind);Y,axis=enclosure(model,centers[i],S[i],centers[j],T[j],mid,rad);examined+=1;depth=max(depth,len(address)//2)
                    if good(Y,kind,raw['orientation'][e]):
                        accepted+=1;minimum=min(minimum,axis.a);handle.write(json.dumps(dict(path=address,image=[str(v) for v in Y],image_binary=[[list(t) for t in v._mpi_] for v in Y],minimum_axis=str(axis.a),radius_binary=list(axis._mpi_[0])),separators=(',',':'))+'\n')
                    else:
                        physical=[float(rad[a])*sum(abs(raw['bases'][i][k][a]*raw['radii'][i][a]) for k in range(3)) for a in range(3)];split=max(range(3),key=lambda a:physical[a]);assert rad[split]>Fraction(1,2**30);stack.extend([address+str(split)+'0',address+str(split)+'1']);refinements+=1
                    if time.time()-last>20:print(json.dumps(dict(edge=e,kind=kind,examined=examined,pending=len(stack),seconds=time.time()-begin)),flush=True);last=time.time()
            row=dict(edge=e,source=i,target=j,kind=kind,status='passed',pending=0,primary_leaves=len(paths),recomputation_calls=examined,refinements=refinements,accepted=accepted,examined=2*accepted-1,maximum_depth=depth,minimum_axis=str(minimum),certificate_sha256=sha(outpath));checks.append(row);print(json.dumps(row),flush=True);(args.out/'progress.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8',newline='\n')
    result=dict(status='passed',interval_digits=args.digits,mpmath=mp.__version__,input_sha256=sha(args.out/'input.json'),primary_report_sha256=sha(args.source/'report.json'),sources={Path(__file__).name:sha(args.out/Path(__file__).name)},checks=checks,seconds=time.time()-begin,scope='Independent full eight-site squared-distance formula, nested spherical reflections, forward automatic derivatives and higher precision. Every primary leaf is recomputed; extra subdivisions retained. Shared mpmath interval library, windows and mean-value principle. No trapping, attraction or historical novelty claim.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--digits',type=int,default=80);p.add_argument('--mpmath-path',type=Path);main(p.parse_args())
