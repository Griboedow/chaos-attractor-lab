"""Second interval realization: scalar circle coordinates and automatic derivatives.

Uses the same subdivision/certificate driver and mpmath interval library, but
replaces the map, derivative propagation and matrix inversion implementations.
"""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
import certify_raar_coverings as driver

def exact(v):
    n,d=float(v).as_integer_ratio();return mp.iv.mpf(n)/d

class Dual:
    def __init__(self,v,d=None):self.v=v;self.d=[mp.iv.mpf(0)]*3 if d is None else d
    def __add__(self,other):
        if not isinstance(other,Dual):other=Dual(other)
        return Dual(self.v+other.v,[a+b for a,b in zip(self.d,other.d)])
    __radd__=__add__
    def __neg__(self):return Dual(-self.v,[-v for v in self.d])
    def __sub__(self,other):return self+-other
    def __rsub__(self,other):return other+-self
    def __mul__(self,other):
        if not isinstance(other,Dual):other=Dual(other)
        return Dual(self.v*other.v,[a*other.v+self.v*b for a,b in zip(self.d,other.d)])
    __rmul__=__mul__
    def __truediv__(self,other):
        if not isinstance(other,Dual):other=Dual(other)
        v=self.v/other.v;return Dual(v,[(a-v*b)/other.v for a,b in zip(self.d,other.d)])
    def __rtruediv__(self,other):return Dual(other)/self
    def sqrt(self):
        v=mp.iv.sqrt(self.v);return Dual(v,[a/(2*v) for a in self.d])

def square(x):
    if isinstance(x,Dual):return Dual(x.v**2,[2*x.v*a for a in x.d])
    return x**2

def sqrt(x):return x.sqrt() if isinstance(x,Dual) else mp.iv.sqrt(x)

class Model:
    def __init__(self,parameters):
        self.beta,self.radius,self.theta,*self.center=[exact(v) for v in parameters];self.ct=mp.iv.cos(self.theta);self.st=mp.iv.sin(self.theta)
    def field(self,x):
        cx,cy,cz=self.center;dx=x[0]-cx;dy=x[1]-cy;dz=x[2]-cz;along=dy*self.ct+dz*self.st;lb=sqrt(square(dx)+square(along));b=[dx/lb*self.radius+cx,along/lb*self.radius*self.ct+cy,along/lb*self.radius*self.st+cz]
        rb=[2*b[i]-x[i] for i in range(3)];la=sqrt(square(rb[0])+square(rb[1]));a=[rb[0]/la,rb[1]/la,rb[2]*0]
        ra=[2*a[i]-rb[i] for i in range(3)];y=[(ra[i]+x[i])*self.beta/2+b[i]*(1-self.beta) for i in range(3)]
        av=la.v if isinstance(la,Dual) else la;bv=lb.v if isinstance(lb,Dual) else lb
        if not(av.a>0 and bv.a>0):raise ArithmeticError('Axis enclosure contains zero')
        return y,min(av.a,bv.a)

def exact_inverse(A):
    def rational(v):
        lo,hi=v._mpi_;assert lo==hi
        sign,mantissa,exponent,_=lo
        return Fraction((-1 if sign else 1)*mantissa)*(Fraction(2)**exponent)
    a=[[rational(v) for v in row]+[Fraction(i==j) for j in range(3)] for i,row in enumerate(A)]
    for k in range(3):
        pivot=next(i for i in range(k,3) if a[i][k]);a[k],a[pivot]=a[pivot],a[k];v=a[k][k];a[k]=[x/v for x in a[k]]
        for i in range(3):
            if i!=k:v=a[i][k];a[i]=[x-v*y for x,y in zip(a[i],a[k])]
    return [[mp.iv.mpf(v.numerator)/v.denominator for v in row[3:]] for row in a]

def enclose(model,source,target,mid,rad):
    iv=mp.iv;S=source['matrix'];q=[exact(v) for v in mid];delta=[iv.mpf([-exact(v).b,exact(v).b]) for v in rad];center=[source['center'][i]+sum(S[i][j]*q[j] for j in range(3)) for i in range(3)];derivative=S;minimum=iv.mpf('inf')
    for _ in range(7):
        box=[Dual(center[i]+sum(derivative[i][j]*delta[j] for j in range(3)),derivative[i]) for i in range(3)];images,axis=model.field(box);minimum=min(minimum,axis);center,_=model.field(center);derivative=[v.d for v in images]
    C=target['inverse'];y=[]
    for i in range(3):
        base=sum(C[i][k]*(center[k]-target['center'][k]) for k in range(3));gradient=[sum(C[i][k]*derivative[k][j] for k in range(3)) for j in range(3)];y.append(base+sum(gradient[j]*delta[j] for j in range(3)))
    return y,minimum

def main(args):
    driver.Model=Model;driver.inverse=exact_inverse;driver.enclose=enclose;driver.main(args)
    path=args.out/'report.json';report=json.loads(path.read_text());report['second_implementation']='Scalar orthonormal circle coordinates; nested reflectors; forward interval automatic differentiation; exact rational Gaussian matrix inversion. Shared subdivision driver and interval library.';report['second_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--digits',type=int,default=100);p.add_argument('--edge',type=int);p.add_argument('--max-cells',type=int,default=200000);main(p.parse_args())
