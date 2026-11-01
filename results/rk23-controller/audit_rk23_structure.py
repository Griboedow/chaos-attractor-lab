"""Exact symbolic checks for the stated accepted-solution bound (SymPy)."""
import argparse,json,hashlib
from pathlib import Path
import sympy as s

def main(out):
    z,h,w=s.symbols('z h w',real=True)
    A=s.Matrix([[0,1],[-1,-2*z]]);P=s.Matrix([[1,z],[z,1]])
    assert A.T*P+P*A==-2*z*P
    R=s.eye(2)+h*A+h*h*A*A/2+h**3*A**3/6
    E=(h**3*A**3+h**4*A**4)/48
    lam=-z+s.I*w
    def reduce(expr):
        return s.rem(s.Poly(s.expand(expr),w),s.Poly(w*w+z*z-1,w)).as_expr().expand()
    def modulus(poly):return reduce(poly*s.conjugate(poly))
    r2=modulus(1+h*lam+(h*lam)**2/2+(h*lam)**3/6)
    e2=modulus((h*lam)**3*(1+h*lam)/48)
    assert (R.T*P*R-r2*P).applyfunc(s.expand)==s.zeros(2)
    assert (E.T*P*E-e2*P).applyfunc(s.expand)==s.zeros(2)
    expansion=-2*z*h*(1-z*h+s.Rational(2,3)*z*z*h*h)+(s.Rational(2,3)*z*z-s.Rational(1,12))*h**4-z*h**5/6+h**6/36
    assert s.expand(r2-1-expansion)==0
    # Nonnegative decomposition of U - (|R|^2-1), 0<z<1, 0<h<=1.
    upper=-s.Rational(7,12)*z*h-h**4/18
    certificate=s.Rational(4,3)*z*h*(z*h-s.Rational(3,4))**2+s.Rational(2,3)*z*h*(1-z*h**3)+z*h**5/6+(h**4-h**6)/36
    assert s.expand(upper-(r2-1)-certificate)==0
    assert s.expand(e2-h**6*((h-z)**2+1-z*z)/2304)==0
    # 8/3 h^3 dominates the coefficient triangle bound for h>=1.
    assert s.expand(s.Rational(8,3)*h**3-(1+h+h*h/2+h**3/6)-(h-1)*(5*h*h+4*h+2)/2)==0
    c=128*s.sqrt(2);a=s.sqrt(2*(1+z));b=s.sqrt(2*(1+z)/(1-z));K=128/s.sqrt(1-z*z)
    # Squared identities avoid unstated square-root branch manipulations.
    assert s.factor((K*a)**2-c*c/(1-z))==0
    assert s.factor((K*b)**2-c*c/(1-z)**2)==0
    # A sharper accepted-solution bound for the explicitly declared O526 pair.
    # Rational(float) is the exact binary64 value, not a decimal approximation.
    z0=s.Rational(.1717106344755739);rho0=s.Rational(.03929059404001495);H=s.Rational(23,10)
    poly=s.Poly(s.cancel((r2-1)/h).subs(z,z0),h)
    root_count=poly.count_roots(0,H)
    assert root_count==0 and poly.eval(0)<0 and poly.eval(H)<0
    B=H**-3+H**-2+1/(2*H)+s.Rational(1,6)
    K2=48**2*B**2/(1-2*z0*H+H**2)
    gain2=s.factor(K2*2*(1+z0)/(1-z0)*rho0**2)
    assert 0<gain2<1
    C=s.sqrt(K2)*s.sqrt(2*(1+z0))/(1-s.sqrt(gain2))
    special=dict(id=526,zeta=str(z0),rtol=str(rho0),H=str(H),polynomial_descending_coefficients=[str(v) for v in poly.all_coeffs()],
                 sturm_roots_on_closed_interval=int(root_count),endpoint_values=[str(poly.eval(0)),str(poly.eval(H))],
                 squared_controller_gain=str(gain2),squared_controller_gain_decimal=float(gain2),normalized_P_bound_constant=float(C),
                 rationale='The exact polynomial has no root on [0,23/10] and is negative at 0. Above H, the triangle-ratio numerator decreases and its denominator increases. The same accepted-solution proof applies with this sharper K.')
    report=dict(status='passed',sympy=s.__version__,checks=['Lyapunov metric identity','matrix R and E norm identities','stability polynomial expansion','nonnegative small-step bound certificate','large-step error denominator','large-step solution triangle bound','positive square-root constant identities','O526 exact Sturm small-step certificate and rational gain inequality'],
                special_case=special,
                small_step_certificate=str(certificate),scope='Exact symbolic identities supporting the written proof. Positivity uses its explicit 0<zeta<1 and h domains. No invariant-set or chaos proof.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)
