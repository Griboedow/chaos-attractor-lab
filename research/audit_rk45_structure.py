"""Exact existing-pair polynomials and accepted-solution bound certificates.

This audits algebra supporting a written proof, not chaos or historical novelty.
The rectangular Bernstein certificate covers all 0<=zeta<=1 and 0<=h<=1/2.
Selected binary64 parameter values are treated as exact rationals.
"""
import argparse,hashlib,json
from pathlib import Path
import sympy as s


def tableau():
    cv=lambda row:list(map(s.Rational,row.split()))
    A=[[],cv('1/5'),cv('3/40 9/40'),cv('44/45 -56/15 32/9'),cv('19372/6561 -25360/2187 64448/6561 -212/729'),cv('9017/3168 -355/33 46732/5247 49/176 -5103/18656')]
    return A,cv('35/384 0 500/1113 125/192 -2187/6784 11/84'),cv('-71/57600 0 71/16695 -71/1920 17253/339200 -22/525 1/40')


def modulus_squared(coefficients,h,zeta):
    n=len(coefficients)
    return s.expand(sum(coefficients[i]**2*h**(2*i) for i in range(n))+2*sum(coefficients[i]*coefficients[j]*h**(i+j)*s.chebyshevt(i-j,-zeta) for i in range(n) for j in range(i)))


def main(args):
    x,h,z=s.symbols('x h z');A,B,E=tableau();stages=[s.Integer(1)]
    for row in A[1:]:stages.append(s.expand(1+x*sum(a*k for a,k in zip(row,stages))))
    R=s.expand(1+x*sum(b*k for b,k in zip(B,stages)));stages.append(R)
    error=s.expand(x*sum(e*k for e,k in zip(E,stages)))
    assert R==sum(x**i/s.factorial(i) for i in range(6))+x**6/600
    assert s.factor(error-x**5*(5*x*x-39*x+97)/120000)==0
    coefficients=[R.coeff(x,k) for k in range(7)];r2=modulus_squared(coefficients,h,z)
    q2=modulus_squared([s.Integer(97),s.Integer(-39),s.Integer(5)],h,z)
    assert s.expand(q2-(25*h**4+390*z*h**3+(551+1940*z*z)*h*h+7566*z*h+9409))==0
    p=s.Poly(s.expand((1-r2)/h),h,z);m=p.degree(h);n=p.degree(z);H=s.Rational(1,2)
    bernstein=[]
    for i in range(m+1):
        row=[]
        for j in range(n+1):
            value=sum(a*H**k*s.binomial(i,k)/s.binomial(m,k)*s.binomial(j,l)/s.binomial(n,l) for (k,l),a in p.terms() if k<=i and l<=j)
            assert value>=0
            row.append(value)
        bernstein.append(row)
    assert bernstein[0][1]>0 and bernstein[-1][-1]>0
    # Reconstruct the exact polynomial, so coefficient indexing is independently checked.
    u,v=s.symbols('u v')
    reconstructed=sum(bernstein[i][j]*s.binomial(m,i)*u**i*(1-u)**(m-i)*s.binomial(n,j)*v**j*(1-v)**(n-j) for i in range(m+1) for j in range(n+1))
    assert s.Poly(s.expand(reconstructed-p.as_expr().subs({h:H*u,z:v})),u,v).is_zero
    broadK=s.factor(24000*sum(a*H**(k-7) for k,a in enumerate(coefficients)))
    selected=[]
    cases=json.loads(args.cases.read_text())['cases'] if args.cases else []
    for case in cases:
        z0,rho=[s.Rational(float(value)) for value in case['params']]
        polynomial=s.Poly(s.expand((r2-1)/h).subs(z,z0),h)
        # Choose the largest tenth-grid endpoint in [0.5,6] certified wholly contractive.
        valid=[]
        for numerator in range(5,61):
            end=s.Rational(numerator,10)
            if polynomial.eval(0)<0 and polynomial.eval(end)<0 and polynomial.count_roots(0,end)==0:valid.append(end)
        assert valid
        end=max(valid);triangleK=s.factor(24000*sum(a*end**(k-7) for k,a in enumerate(coefficients)))
        K=triangleK;shifted_coefficients=None
        # A positive-coefficient polynomial in t=h-H certifies the ratio on
        # the entire unbounded interval, without a sampled-grid inference.
        for candidate in [s.Rational(1,4),s.Rational(1,2),s.Integer(1),s.Integer(2),s.Integer(4),s.Integer(8),s.Integer(16),s.Integer(32),s.Integer(64),s.Integer(128),s.Integer(256)]:
            t=s.symbols('t')
            diff=s.Poly(s.expand((candidate**2*h**10*q2/120000**2-r2).subs({z:z0,h:end+t})),t)
            coeffs=diff.all_coeffs()
            if diff.eval(0)>0 and all(v>=0 for v in coeffs):
                K=candidate;shifted_coefficients=coeffs;break
        b2=2*(1+z0)/(1-z0);gain2=s.factor(K*K*b2*rho*rho)
        bound=K*s.sqrt(2*(1+z0))/(1-s.sqrt(gain2)) if gain2<1 else None
        selected.append(dict(id=case['id'],zeta=str(z0),rtol=str(rho),H=str(end),sturm_roots=0,polynomial_descending_coefficients=[str(v) for v in polynomial.all_coeffs()],endpoint_values=[str(polynomial.eval(0)),str(polynomial.eval(end))],triangle_K=str(triangleK),K=str(K),ratio_certificate='Every coefficient of K^2 |E((H+t)lambda)|^2-|R((H+t)lambda)|^2 is nonnegative, with positive constant.' if shifted_coefficients is not None else 'Triangle bound',ratio_polynomial_descending_coefficients=[str(v) for v in shifted_coefficients] if shifted_coefficients is not None else None,squared_gain=str(gain2),squared_gain_decimal=float(gain2),gain_condition=bool(gain2<1),normalized_P_bound_constant=float(bound) if bound is not None else None))
    report=dict(status='passed',sympy=s.__version__,R=str(R),E=str(error),R_modulus_squared=str(r2),error_quadratic_modulus_squared=str(q2),
                bernstein_certificate=dict(polynomial='(1-|R(h lambda)|^2)/h',h_interval=['0','1/2'],zeta_interval=['0','1'],degree=[m,n],coefficients=[[str(v) for v in row] for row in bernstein],negative=0,zero=sum(v==0 for row in bernstein for v in row),reconstruction_identity=True),
                broad_ratio_bound=str(broadK),selected_cases=selected,
                scope='Exact algebra and accepted-solution sufficient bounds only. No bound on the proposed step, trapping region or chaos proof.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',broad_ratio_bound=str(broadK),selected=[{k:r[k] for k in ['id','H','squared_gain_decimal','gain_condition','normalized_P_bound_constant']} for r in selected])),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--cases',type=Path);main(p.parse_args())
