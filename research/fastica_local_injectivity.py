"""Exact algebra for local noninjectivity at the structured FastICA attractor.

The global topological argument is in docs/fastica-tripling-structure.md.
This program checks symbolic identities and rational uniform bounds; numerical
periodic-point coordinates and finite differences are supplementary controls.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys
import sympy as sp
import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def positive_polynomial_certificate(expr,e,v,vlo,vhi):
    """Exact Bernstein coefficients on e^2 in [0,1/625], v in [vlo,vhi]."""
    z,u,b=sp.symbols('z u b')
    powers=sp.Poly(sp.expand(sp.cancel(expr)),e,v)
    assert all(i%2==0 for (i,j),c in powers.terms())
    poly=sum(c*z**(i//2)*v**j for (i,j),c in powers.terms())
    unit=sp.Poly(sp.expand(poly.subs({z:u/sp.Integer(625),v:vlo+(vhi-vlo)*b})),u,b)
    n,m=unit.degree(u),unit.degree(b);coeff=[];reconstructed=0
    for k in range(n+1):
        for l in range(m+1):
            value=sum(c*sp.binomial(k,i)/sp.binomial(n,i)*sp.binomial(l,j)/sp.binomial(m,j)
                      for (i,j),c in unit.terms() if i<=k and j<=l)
            coeff.append(value)
            reconstructed+=value*sp.binomial(n,k)*u**k*(1-u)**(n-k)*sp.binomial(m,l)*b**l*(1-b)**(m-l)
    assert sp.expand(reconstructed-unit.as_expr())==0
    assert min(coeff)>0
    return {'polynomial_in_e_squared_and_v':str(poly),'v_interval':[vlo,vhi],
      'bernstein_degrees':[n,m],'coefficients':list(map(str,coeff)),'strict_lower_bound':str(min(coeff))}

def check():
    A,B,C,D,e,t=sp.symbols('A B C D e t',real=True)
    variables=[A,B,C,D]
    P=(A**4+B**4-6*A*A*B*B+24*e*(A**3*C+B**3*D)
       +e**2*(72*(A*A*C*C+B*B*D*D)-36*(A*A*D*D+B*B*C*C))
       +72*e**3*(A*C**3+B*D**3)-18*e**4*(C**4+D**4)-216*e**4*C*C*D*D)
    prior=ROOT/'results/fastica-tripling/exact-check.json'
    frozen=json.loads(prior.read_text())
    assert sp.expand(P-sp.sympify(frozen['cumulant_polynomial'],locals=dict(zip(['A','B','C','D','e'],[A,B,C,D,e]))))==0
    h=sp.Matrix([sp.diff(P,x)/4 for x in variables]);H=h.jacobian(variables)
    assert sp.simplify(H*sp.Matrix(variables)-3*h)==sp.zeros(4,1)
    q_axis=18*e**3*t**4+(36*e**2+18*e**4)*t**3+(18*e-54*e**3)*t**2+(1-36*e**2)*t-6*e
    q_diag=18*e**3*t**4+(18*e**2+126*e**4)*t**3+(18*e-54*e**3)*t**2-(2+18*e**2)*t-6*e
    d_axis=-2916*e**4*(1+6*e**2*t*t)**3*(4+18*e*t+15*e**2*t*t)
    d_diag=46656*e**4*(1+6*e*t+3*e**2*t*t)**3*(2-6*e*t+15*e**2*t*t)
    parts=[]
    for name,subs,q,det,bracket in [
      ('axis',{A:1,B:0,C:t,D:0},q_axis,d_axis,[4,6]),
      ('diagonal',{A:1,B:1,C:t,D:t},q_diag,d_diag,[-4,-2])]:
        hs=h.subs(subs);assert sp.expand(t*hs[0]-hs[2]-q)==0
        assert sp.expand(H.subs(subs).det()-det)==0
        if name=='axis':assert hs[1]==0 and hs[3]==0
        else:assert sp.expand(hs[0]-hs[1])==0 and sp.expand(hs[2]-hs[3])==0
        parts.append({'family':name,'q':str(q),'h_first':str(hs[0]),'det_raw_jacobian':str(det),
                      't_over_e_bracket':bracket,'endpoint_q_over_e':[str(sp.factor(q.subs(t,k*e)/e)) for k in bracket]})
    # Each bound below follows by dropping terms of a known sign and using e<=E.
    E=sp.Rational(1,25)
    bounds={
      'axis_left_negative_margin':1-72*E**2-720*E**4-2880*E**6,
      'axis_right_q_over_e3_lower':sp.Integer(432),
      'diagonal_left_positive_half_margin':1-1008*E**4-1728*E**6,
      'diagonal_right_negative_half_margin':1-54*E**2,
      'axis_q_derivative_lower':1-36*E**2-648*E**4,
      'diagonal_q_derivative_negative_margin':2-1296*E**4-6048*E**6,
      'diagonal_h_first_negative_margin':2-288*E**4,
      'diagonal_determinant_factor_lower':1-24*E**2,
      'axis_tube_margin':sp.Rational(2,5)-6*E,
      'diagonal_tube_margin':sp.Rational(2,5)-4*E}
    assert all(v>0 for v in bounds.values())
    # Verify the endpoint identities used to make the uniform sign deductions.
    assert sp.expand(q_axis.subs(t,4*e)/e-2*(2880*e**6+720*e**4+72*e**2-1))==0
    assert sp.expand(q_axis.subs(t,6*e)/e-216*e**2*(126*e**4+27*e**2+2))==0
    assert sp.expand(q_diag.subs(t,-4*e)/e-2*(1+180*e**2-1008*e**4-1728*e**6))==0
    assert sp.expand(q_diag.subs(t,-2*e)/e+2*(1-54*e**2+180*e**4+360*e**6))==0
    v=sp.symbols('v')
    positive_expressions={
      'negative_axis_left':(-q_axis.subs(t,4*e)/(2*e),0,1),
      'positive_axis_right':(q_axis.subs(t,6*e)/e**3,0,1),
      'positive_diagonal_left':(q_diag.subs(t,-4*e)/(2*e),0,1),
      'negative_diagonal_right':(-q_diag.subs(t,-2*e)/(2*e),0,1),
      'axis_unique_root':(sp.diff(q_axis,t).subs(t,e*v),4,6),
      'diagonal_unique_root':(-sp.diff(q_diag,t).subs(t,e*v),-4,-2),
      'negative_diagonal_multiplier':(-h[0].subs({A:1,B:1,C:e*v,D:e*v}),-4,-2),
      'positive_diagonal_det_factor':(1+6*e**2*v+3*e**4*v*v,-4,-2),
      'positive_axis_multiplier':(h[0].subs({A:1,B:0,C:e*v,D:0}),4,6)}
    certificates={name:positive_polynomial_certificate(expr,e,v,lo,hi) for name,(expr,lo,hi) in positive_expressions.items()}
    return {'scope':'Exact symbolic identities and rational uniform margins; the topological argument is not machine formalized.',
      'epsilon_interval':'0 < epsilon <= 1/25','symbolic_cases':parts,
      'positive_bounds':{k:{'exact':str(v),'decimal':float(v)} for k,v in bounds.items()},
      'uniform_positivity_certificates':certificates,
      'euler_identity':'Dh(w) w = 3 h(w)',
      'oriented_tangent_determinant':'det(Dh(w)) / (3 ||h(w)||^4)',
      'previous_algebra_sha256':hashlib.sha256(prior.read_bytes()).hexdigest()},(e,t,q_axis,q_diag)

def oriented_frame(w):
    Q=np.linalg.qr(np.column_stack([w,np.eye(4)]))[0]
    if Q[:,0]@w<0:Q[:,0]*=-1
    if np.linalg.det(Q)<0:Q[:,1]*=-1
    assert np.max(abs(Q[:,0]-w))<1e-14 and np.linalg.det(Q)>.999999999
    return Q[:,1:]

def numerical_controls(polynomials):
    e,t,qa,qd=polynomials;mp.mp.dps=90
    spec=importlib.util.spec_from_file_location('frozen_precision',ROOT/'results/fastica-tripling/code/precision_check.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=[]
    for ep in ['0.01','0.02','0.03','0.04']:
        E=mp.mpf(ep)
        for family,q,limits in [('axis',qa,[4,6]),('diagonal',qd,[-4,-2])]:
            fn=sp.lambdify((e,t),q,modules='mpmath');lo,hi=[k*E for k in limits]
            sign=mp.sign(fn(E,lo))
            for _ in range(260):
                mid=(lo+hi)/2
                if mp.sign(fn(E,mid))==sign:lo=mid
                else:hi=mid
            tau=(lo+hi)/2;w=mp.matrix([1,0,tau,0] if family=='axis' else [1,1,tau,tau]);w/=mp.norm(w)
            h,H=module.raw(w,E);hn=mp.norm(h);f=h/hn;df=(mp.eye(4)-f*f.T)*H/hn
            target=w if family=='axis' else -w;res=mp.norm(f-target)
            wb=np.array(list(w),float);fb=np.array(list(f),float);V=oriented_frame(wb);U=oriented_frame(fb)
            J=np.array(df.tolist(),float);oriented=np.linalg.det(U.T@J@V);formula=mp.det(H)/(3*hn**4)
            assert np.sign(oriented)==(-1 if family=='axis' else 1)
            assert abs(oriented-float(formula))<1e-12 and res<mp.mpf('1e-70')
            # Centered sphere finite differences use the independent explicit cubic.
            deltas=[]
            for step in [1e-5,1e-6]:
                columns=[]
                for v in V.T:
                    p=mp.matrix((wb+step*v).tolist());p/=mp.norm(p);m=mp.matrix((wb-step*v).tolist());m/=mp.norm(m)
                    hp,_=module.raw(p,E);hm,_=module.raw(m,E);column=(hp/mp.norm(hp)-hm/mp.norm(hm))/(2*step)
                    columns.append(np.array(list(column),float))
                finite=U.T@np.array(columns).T;error=float(abs(finite-U.T@J@V).max())
                assert error<1e-8 and np.sign(np.linalg.det(finite))==np.sign(oriented)
                deltas.append({'step':step,'determinant':float(np.linalg.det(finite)),'max_entry_error':error})
            rows.append({'epsilon':ep,'family':family,'t':str(tau),'weight':list(map(str,w)),'periodic_residual_90_digits':str(res),
              'oriented_tangent_determinant':str(formula),'determinant_via_oriented_frames':oriented,'finite_differences':deltas})
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
    if out.is_relative_to(ROOT):p.error('Write new results outside the repository, normally on D:.')
    out.mkdir(parents=True,exist_ok=False)
    result,polys=check();result['numerical_controls']=numerical_controls(polys)
    result.update(status='passed',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=sp.__version__,mpmath_version=mp.__version__,numpy_version=np.__version__)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':'passed','bounds':len(result['positive_bounds']),'periodic_point_controls':len(result['numerical_controls']),
      'epsilon_004':[{'family':r['family'],'t':r['t'][:20],'tangent_determinant':r['determinant_via_oriented_frames']} for r in result['numerical_controls'][-2:]]}))

if __name__=='__main__':main()
