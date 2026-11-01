"""Second fold-region calculation, with scalar forward first/second derivatives.

Independent scalar projector and nested-reflector expressions replace the
vector projector formulas. The doubleton enclosure, subdivision driver and
mpmath interval library are shared. No finite-difference error is used.
"""
import argparse
import json
from pathlib import Path
import mpmath as mp
import certify_raar_critical_fold_region as driver
from check_raar_critical_crossing_ad import Model
from certify_raar_cycles import det


class Jet2:
    def __init__(self,v,d=None,h=None):
        self.v=v;self.d=[mp.iv.mpf(0)]*3 if d is None else d
        self.h=[[mp.iv.mpf(0) for _ in range(3)] for _ in range(3)] if h is None else h

    def __add__(self,y):
        if not isinstance(y,Jet2):y=Jet2(y)
        return Jet2(self.v+y.v,[self.d[i]+y.d[i] for i in range(3)],[[self.h[i][j]+y.h[i][j] for j in range(3)] for i in range(3)])
    __radd__=__add__

    def __neg__(self):return Jet2(-self.v,[-v for v in self.d],[[-v for v in row] for row in self.h])
    def __sub__(self,y):return self+-y
    def __rsub__(self,y):return y+-self

    def __mul__(self,y):
        if not isinstance(y,Jet2):y=Jet2(y)
        return Jet2(self.v*y.v,[self.d[i]*y.v+self.v*y.d[i] for i in range(3)],[[self.h[i][j]*y.v+self.d[i]*y.d[j]+self.d[j]*y.d[i]+self.v*y.h[i][j] for j in range(3)] for i in range(3)])
    __rmul__=__mul__

    def inverse(self):
        v=1/self.v
        return Jet2(v,[-self.d[i]*v**2 for i in range(3)],[[2*self.d[i]*self.d[j]*v**3-self.h[i][j]*v**2 for j in range(3)] for i in range(3)])
    def __truediv__(self,y):return self*(y.inverse() if isinstance(y,Jet2) else 1/y)
    def __rtruediv__(self,y):return self.inverse()*y

    def square(self):
        return Jet2(self.v**2,[2*self.v*a for a in self.d],[[2*(self.d[i]*self.d[j]+self.v*self.h[i][j]) for j in range(3)] for i in range(3)])

    def sqrt(self):
        v=mp.iv.sqrt(self.v)
        if not v.a>0:raise ArithmeticError('Square-root derivative requires a positive interval')
        return Jet2(v,[a/(2*v) for a in self.d],[[self.h[i][j]/(2*v)-self.d[i]*self.d[j]/(4*v**3) for j in range(3)] for i in range(3)])


def jet(x,p,ctx):
    assert ctx is mp.iv
    beta,radius,theta,cx,cy,cz=p;ct,st=ctx.cos(theta),ctx.sin(theta)
    variables=[Jet2(x[i],[ctx.mpf(int(i==j)) for j in range(3)]) for i in range(3)]
    dx=variables[0]-cx;along=(variables[1]-cy)*ct+(variables[2]-cz)*st
    lb=(dx.square()+along.square()).sqrt()
    b=[dx/lb*radius+cx,along/lb*radius*ct+cy,along/lb*radius*st+cz]
    reflected=[2*b[i]-variables[i] for i in range(3)]
    la=(reflected[0].square()+reflected[1].square()).sqrt()
    a=[reflected[0]/la,reflected[1]/la,Jet2(ctx.mpf(0))]
    ra=[2*a[i]-reflected[i] for i in range(3)]
    values=[(ra[i]+variables[i])*beta/2+b[i]*(1-beta) for i in range(3)]
    J=[v.d for v in values];H=[v.h for v in values]
    # Differentiate the determinant by its six signed permutation products.
    permutations=[((0,1,2),1),((1,2,0),1),((2,0,1),1),((0,2,1),-1),((1,0,2),-1),((2,1,0),-1)]
    gradient=[]
    for k in range(3):
        g=ctx.mpf(0)
        for permutation,sign in permutations:
            for row in range(3):
                term=sign*H[row][permutation[row]][k]
                for other in range(3):
                    if other!=row:term*=J[other][permutation[other]]
                g+=term
        gradient.append(g)
    return det(J),gradient,J,H,[la.v,lb.v]


def main(args):
    driver.Model=Model;driver.jet=jet
    extra=['check_raar_critical_fold_ad.py','check_raar_critical_crossing_ad.py','check_raar_coverings_ad.py']
    try:
        driver.main(args)
    finally:
        if args.out.exists():
            for name in extra:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    path=args.out/'report.json';report=json.loads(path.read_text())
    report['second_implementation']='Scalar circle coordinates and nested reflectors with forward first/second automatic derivatives; determinant gradient differentiated by permutation products. Shared enclosure, partition driver and interval library.'
    report['second_sources']={name:driver.sha(args.out/name) for name in extra}
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'second_implementation':report['second_implementation']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--window',type=int,default=0);p.add_argument('--steps',type=int,default=42);p.add_argument('--digits',type=int,default=100);p.add_argument('--max-cells',type=int,default=20000);main(p.parse_args())
