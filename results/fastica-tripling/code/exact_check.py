"""Exact moment identities and rational bounds for the FastICA construction."""
from pathlib import Path
import hashlib,json
import sympy as s
import numpy as np
from model import EXP,coefficients

A,B,C,D,e=s.symbols('A B C D e',real=True)
w=[A,B,C,D];base=[]
for r in [s.Integer(1)]*4+[s.sqrt(6)]:
    base += [(r,0),(-r,0),(0,r),(0,-r)]
vectors=[s.Matrix([x,y,e*x*(x*x+y*y-4),e*y*(x*x+y*y-4)]) for x,y in base]
mean=sum(vectors,s.zeros(4,1))/20
cov=sum((v*v.T for v in vectors),s.zeros(4))/20+s.diag(0,0,1-6*e**2,1-6*e**2)
assert mean==s.zeros(4,1) and s.simplify(cov-s.eye(4))==s.zeros(4)
V=[(v.T*s.Matrix(w))[0] for v in vectors]
V2=s.expand(sum(v**2 for v in V)/20)
V4=s.expand(sum(v**4 for v in V)/20)
# The independent six-point quadratures have moments (1,0,1,0,3).
noise=[s.sqrt(3),-s.sqrt(3),s.Integer(0),s.Integer(0),s.Integer(0),s.Integer(0)]
moments=[s.simplify(sum(x**j for x in noise)/6) for j in range(5)]
assert moments==[1,0,1,0,3]
N2=(1-6*e**2)*(C*C+D*D)
P=s.expand(V4+6*V2*N2+3*N2*N2-3*sum(x*x for x in w)**2)
coef=[1,1,-6,24*e,24*e,72*e**2,72*e**2,-36*e**2,-36*e**2,72*e**3,72*e**3,-18*e**4,-18*e**4,-216*e**4]
expected=sum(c*s.prod(x**int(n) for x,n in zip(w,p)) for c,p in zip(coef,EXP))
assert s.expand(P-expected)==0
for ep in [s.Rational(1,25),s.Rational(3,20)]:
    assert np.allclose(coefficients(float(ep)),[float(s.sympify(c).subs(e,ep)) for c in coef],rtol=1e-15,atol=0)
h=[s.diff(P,x)/4 for x in w]
assert s.expand(h[0].subs(e,0)-(A**3-3*A*B**2))==0
assert s.expand(h[1].subs(e,0)-(B**3-3*B*A**2))==0
assert all(s.expand(v.subs(e,0))==0 for v in h[2:])
assert s.expand((A**3-3*A*B**2)**2+(B**3-3*B*A**2)**2-(A*A+B*B)**3)==0
ep=s.Rational(1,25);k=s.Rational(2,5)
L=1-18*ep*k-36*ep**2*k**2-18*ep**3*k**3
U=6*ep+36*ep**2*k+54*ep**3*k**2+108*ep**4*k**3
assert L>0 and U<k*L
result={'status':'passed','scope':'Exact symbolic moment identities and rational endpoint inequalities; topological argument is written in the report, not machine formalized.',
 'rows':720,'mean':str(mean),'covariance':str(s.simplify(cov)),'quadrature_moments':list(map(str,moments)),
 'cumulant_polynomial':str(P),'gradient_quarter':list(map(str,h)),
 'epsilon_interval':['0','1/25'],'tube_ratio':'2/5','L':str(L),'U':str(U),'ratio_bound':str(U/L),
 'strict_margin':str(k*L-U),'decimal_ratio_bound':float(U/L),'sympy_version':s.__version__,
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'model_sha256':hashlib.sha256(Path(__file__).with_name('model.py').read_bytes()).hexdigest()}
Path(__file__).with_name('exact-check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','L','U','ratio_bound','strict_margin']}))
