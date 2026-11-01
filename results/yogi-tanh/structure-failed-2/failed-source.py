"""Exact and numerical controls for the credited tanh-Yogi structure note."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import sympy as s


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    b,r,q,z=s.symbols('b r q z',real=True);A=(1-b)/(1+b)
    L=s.Matrix([[1+2*b-2*r*q,2*b/A,-2*(1-q)/A],[-(1-b),-b,0],[q*r*A,0,1-q]])
    coeff=[1,-b-2+q*(1+2*r),1+2*b+q*(2*b*r-b-1),b*(q-1)]
    poly=sum(c*z**(3-j) for j,c in enumerate(coeff))
    assert all(s.factor(left-right)==0 for left,right in zip(L.charpoly(z).all_coeffs(),coeff))
    a1,a2,a3=coeff[1:];lower=(4*b*r-(1-b)**2)/(b*(1-b+2*r));upper=2*(1+b)/(1+b+r*(1-b))
    identities={
        'p(1)':(poly.subs(z,1),2*q*r*(1+b)),
        '-p(-1)':(-poly.subs(z,-1),4*(1+b)-2*q*((1+b)+r*(1-b))),
        'Jury residual':(1-a2+a1*a3-a3*a3,q*((1-b)**2-4*b*r+b*(1-b+2*r)*q)),
        'nonempty interval':(upper-lower,(1-b)*(1+r)*((1+b)**2-4*b*r)/(b*(1-b+2*r)*(1+b+r*(1-b)))),
        'lower minus one':(lower-1,(b*(1+2*r)-1)/(b*(1-b+2*r)))
    }
    for left,right in identities.values():assert s.factor(left-right)==0
    # The actual derivative, evaluated at the exact cycle, equals L after reflection.
    x,m,v,c,e,k=s.symbols('x m v c e k',positive=True)
    T=s.tanh(k*(v-x*x));M=b*m+(1-b)*x;V=v-(1-c)*T*x*x
    F=s.Matrix([x-M/(s.sqrt(V)+e),M,V]);J=F.jacobian([x,m,v])
    a=s.symbols('a',positive=True)
    at={x:a,m:-A*a,v:a*a}
    J0=s.simplify(J.subs(at));J0=s.simplify(J0.subs(e,A/2-a).subs(k,q/((1-c)*a*a)).subs(a,r*A/2))
    assert all(s.factor(w)==0 for w in s.diag(-1,-1,1)*J0-L)
    # Adam with epsilon outside sqrt and c_eff=1-q has exactly the same local J.
    VA=(1-q)*v+q*x*x;FA=s.Matrix([x-M/(s.sqrt(VA)+e),M,VA])
    JA=s.simplify(FA.jacobian([x,m,v]).subs(at).subs(e,A/2-a).subs(a,r*A/2))
    assert all(s.factor(w)==0 for w in J0-JA)
    rng=np.random.default_rng(2026110165);records=[]
    for i in range(1024):
        bb,rr=rng.uniform(.002,.998,2);qq=rng.uniform(.0001,4.)
        numeric=np.array([float(w.subs({b:bb,r:rr,q:qq})) for w in coeff]);roots=np.roots(numeric)
        lo=float(lower.subs({b:bb,r:rr}));hi=float(upper.subs({b:bb,r:rr}))
        predicted=max(0.,lo)<qq<hi;radius=float(max(abs(roots)))
        assert abs(radius-1)>1e-10 and predicted==(radius<1)
        cc=float(rng.uniform(.001,.999));AA=(1-bb)/(1+bb);aa=rr*AA/2;ee=AA/2-aa;kk=qq/((1-cc)*aa*aa)
        def advance(w):
            xx,mm,vv=w;MM=bb*mm+(1-bb)*xx;VV=vv-(1-cc)*np.tanh(kk*(vv-xx*xx))*xx*xx
            return np.array([xx-MM/(np.sqrt(VV)+ee),MM,VV])
        plus=np.array([aa,-AA*aa,aa*aa]);minus=plus*np.array([-1.,-1.,1.])
        error=float(max(np.max(abs(advance(plus)-minus)),np.max(abs(advance(minus)-plus))))
        assert error<1e-13
        records.append(dict(b=bb,r=rr,q=qq,c=cc,e=ee,k=kk,lower=lo,upper=hi,spectral_radius=radius,stable=predicted,cycle_residual=error))
    result=dict(symbolic_identities={name:s.sstr(right) for name,(_,right) in identities.items()},
        derivative_equals_relative_matrix=True,adam_epsilon_outside_local_derivative_identity=True,
        declared=len(records),completed=len(records),runs=records,
        scope='Controls check algebra and implementation, not historical priority or nonlinear behavior at unit-modulus boundaries. The symmetric Adam cycle and its local Jacobian are credited, not claimed as a new mechanism.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print('PASS: symbolic identities, Adam local equivalence, 1024 cycle and Schur controls',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
