"""Isolate transversal crossings of the multiplier -1 along the certified diagonal."""
import argparse
import json
from pathlib import Path

from flint import arb,arb_mat,ctx
import mpmath as mp
import numpy as np

from period_windows import parameters
from search import save


class Jet:
    """A scalar, its gradient and Hessian with respect to (x, y, z, t)."""
    def __init__(self,value,g,h):
        self.v,self.g,self.h=value,g,h

    def __add__(self,b):
        if isinstance(b,Jet):
            return Jet(self.v+b.v,self.g+b.g,self.h+b.h)
        return Jet(self.v+b,self.g,self.h)

    __radd__=__add__

    def __neg__(self):
        return Jet(-self.v,-self.g,-self.h)

    def __sub__(self,b):
        return self+-b

    def __rsub__(self,b):
        return -self+b

    def __mul__(self,b):
        if not isinstance(b,Jet):
            return Jet(self.v*b,self.g*b,self.h*b)
        cross=np.outer(self.g,b.g)
        return Jet(self.v*b.v,self.g*b.v+b.g*self.v,self.h*b.v+b.h*self.v+cross+cross.T)

    __rmul__=__mul__

    def __truediv__(self,n):
        return Jet(self.v/n,self.g/n,self.h/n)

    def exp(self,op):
        e=op(self.v)
        return Jet(e,self.g*e,(self.h+np.outer(self.g,self.g))*e)


def evaluate(values,period,exp):
    unit=np.identity(4,dtype=int)
    variables=[Jet(v,np.array([v*0+int(i==k) for k in range(4)],dtype=object),
                   np.full((4,4),v*0,dtype=object)) for i,v in enumerate(values)]
    x,y,z,t=variables
    r,a,b,c,s=parameters(t)
    for _ in range(period):
        n=s*x+(1-s)*z*(r*(1-a*x-b*y-c*z)).exp(exp)
        x,y,z=n,x,y
    target=[x,y,z]
    j=np.array([v.g[:3] for v in target],dtype=object)+unit[:3,:3]
    cof=[]
    for i in range(3):
        u,v=(i+1)%3,(i+2)%3
        cof.append([j[u,(k+1)%3]*j[v,(k+2)%3]-j[u,(k+2)%3]*j[v,(k+1)%3] for k in range(3)])
    det=sum(j[0,k]*cof[0][k] for k in range(3))
    row=[sum(cof[i][k]*target[i].h[k,l] for i in range(3) for k in range(3)) for l in range(4)]
    f=[v.v-values[i] for i,v in enumerate(target)]+[det]
    jac=[list(v.g-unit[i]) for i,v in enumerate(target)]+[row]
    return f,jac


def propose(row,guess,period):
    # Connect the stored stable phase to the tested parameter. The final
    # interval verifier uses only the rational center and the period.
    x=np.array([float(q) for q in row["center"]])
    start=row.get("t_numerator",1)/row.get("t_denominator",1)
    for t in np.linspace(start,guess,11)[1:]:
        for _ in range(12):
            f,j=evaluate([*x,t],period,np.exp)
            dx=np.linalg.solve(np.array(j,dtype=float)[:3,:3],np.array(f[:3],dtype=float))
            x-=dx
            if max(abs(dx))<2e-12:
                break
    mp.mp.dps=period+64
    z=mp.matrix([mp.mpf(float(q)) for q in [*x,guess]])
    for _ in range(12):
        f,j=evaluate(z,period,mp.exp)
        delta=mp.lu_solve(mp.matrix(j),mp.matrix(f))
        z-=delta
        if max(abs(v) for v in delta)<mp.power(10,-period-40):
            break
    f,j=evaluate(z,period,mp.exp)
    assert max(abs(v) for v in f)<mp.power(10,-period-32)
    R=mp.matrix(j)**-1
    print("root",period,mp.nstr(z[3],16),flush=True)
    return {"period":period,"center":[mp.nstr(q,period+46) for q in z],
            "inverse":[[mp.nstr(R[i,j],period+46) for j in range(4)] for i in range(4)],
            "radius":"1e-"+str(period+24)}


def verify(row):
    n=row["period"]
    ctx.prec=4*(n+76)
    mid=list(map(arb,row["center"]))
    rad=arb(row["radius"])
    cube=[q+arb(0,row["radius"]) for q in mid]
    assert cube[3]>0 and cube[3]<1
    R=arb_mat([[arb(q) for q in v] for v in row["inverse"]])
    assert not R.det().contains(0)
    f,_=evaluate(mid,n,lambda x:x.exp())
    _,j=evaluate(cube,n,lambda x:x.exp())
    raw_j=arb_mat(j)
    delta=arb_mat([[int(a==b) for b in range(4)] for a in range(4)])-R*raw_j
    residual=R*arb_mat([[q] for q in f])
    margins=[abs(residual[i,0])/rad+sum(abs(delta[i,k]) for k in range(4)) for i in range(4)]
    assert all(q<1 for q in margins)
    # J[:3,:3]=D(F^p)-I. Its nonsingularity yields the local periodic branch.
    ret=arb_mat([[j[i][k] for k in range(3)] for i in range(3)])
    assert not ret.det().contains(0)
    tangent=-(ret.inv()*arb_mat([[j[i][3]] for i in range(3)]))
    slope=j[3][3]+sum(j[3][k]*tangent[k,0] for k in range(3))
    assert not slope.contains(0)
    # At an exact root, the cubic characteristic polynomial factors as
    # (lambda+1)(lambda^2+u*lambda+v). The strict real Jury inequalities
    # place both remaining multipliers inside the unit circle.
    ident=arb_mat([[int(i==k) for k in range(3)] for i in range(3)])
    deriv=ret+ident
    u=-1-sum(deriv[k,k] for k in range(3))
    v=-deriv.det()
    assert 1+u+v>0 and 1-u+v>0 and 1-v>0
    x,y,z,t=cube
    r,a,b,c,s=parameters(t)
    for h in range(1,n):
        x,y,z=s*x+(1-s)*z*(r*(1-a*x-b*y-c*z)).exp(),x,y
        if n%h==0:
            assert any(not q.overlaps(cube[k]) for k,q in enumerate((x,y,z)))
    return {"precision_bits":ctx.prec,"maximum_inclusion_ratio":str(max(margins)),
            "det_plus_I_derivative_along_branch":str(slope),
            "remaining_multipliers_in_open_unit_disk":True,"least_period":n}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--periods",type=int,nargs="+",default=[39,78,156])
    ap.add_argument("--verify",action="store_true")
    args=ap.parse_args()
    if args.verify:
        for row in json.loads(args.out.read_text()):
            verify(row)
    else:
        data=[json.loads((args.from_dir/"period-39.json").read_text())]
        data+=json.loads((args.from_dir/"period-windows.json").read_text())
        data+=json.loads((args.from_dir/"high-period-windows.json").read_text())
        guesses={39:.9941,78:.98582,156:.98424,312:.9839,624:.98382}
        results=[]
        for n in args.periods:
            item=next(s for s in data if s["period"]==n)
            row=propose(item,guesses[n],n)
            row["verified"]=verify(row)
            results.append(row)
            save(args.out,results)
            print("certified transversal multiplier crossing",n,flush=True)
