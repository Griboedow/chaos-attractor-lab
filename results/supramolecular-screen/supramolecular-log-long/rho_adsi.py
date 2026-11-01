"""Literal homogeneous ADSI kinetics, Chomchai et al. (2025), Eq. 7/Table S2.

Concentrations are (RT,RD,F), in micromolar. Computational time tau=.025*t,
where t is seconds. Parameters are GEF, GAP, inhibitor speed, exchange speed.
The two speeds multiply complete production/removal pairs, preserving their
nullclines. Source values have both speeds equal to one. No diffusion is added.
"""
import numpy as np

BASE=np.array([.00625,.3125,1.,.0625,.05625,.0625,.0208333,.001875,.1406,.25,.025])
REFERENCE=np.array([1.,1.05,1.,1.])
TIME_RATE=.025
SOURCE='https://biology.ed.ac.uk/sites/default/files/2025-04/70-Goryachev-Fcurbiol0325.pdf'


def field_jac(state,p):
    x,y,z=state;gef,gap,speed,exchange=p
    k=BASE.copy();k[5:7]*=exchange;k[[7,8,10]]*=speed
    a=k[0]+gef*k[1]*x**3/(1+k[2]*x*x)
    da=gef*k[1]*x*x*(3+k[2]*x*x)/(1+k[2]*x*x)**2
    b=k[3]+k[4]*gap*z
    reaction=a*y-b*x
    field=np.array([reaction,-reaction+k[5]-k[6]*y,k[7]+k[8]*x*x/(1+k[9]*x*x)-k[10]*z])
    j=np.array([[da*y-b,a,-k[4]*gap*x],[-da*y+b,-a-k[6],k[4]*gap*x],
                [2*k[8]*x/(1+k[9]*x*x)**2,0.,-k[10]]])
    return field/TIME_RATE,j/TIME_RATE


def independent_field(state,p):
    """Independent Eq. 8 total-Rho coordinates, then transform back to Eq. 7."""
    active,inactive,actin=state;gef,gap,speed,exchange=p
    total=active+inactive
    active_dot=(.00625+gef*.3125*active**3/(1+active**2))*(total-active)-(.0625+.05625*gap*actin)*active
    total_dot=exchange*(.0625-.0208333*(total-active))
    actin_dot=speed*(.001875+.1406*active**2/(1+.25*active**2)-.025*actin)
    return np.array([active_dot,total_dot-active_dot,actin_dot])/.025


def equilibria(p):
    """Positive stationary states from the exact fifth-degree numerator.

    Numerical root proposals, not validated roots or a global boundedness proof.
    The two paired speed factors do not affect the stationary states.
    """
    gef,gap,_,_=p;k=BASE;q=k[5]/k[6]
    a=np.polynomial.Polynomial([1.,0.,k[2]])
    b=np.polynomial.Polynomial([1.,0.,k[9]])
    x=np.polynomial.Polynomial([0.,1.])
    numerator=q*(k[0]*a+gef*k[1]*x**3)*b-x*(k[3]*a*b+k[4]*gap/k[10]*a*(k[7]*b+k[8]*x*x))
    roots=numerator.roots();answer=[]
    for root in roots:
        if abs(root.imag)<1e-8 and root.real>0:
            active=root.real;actin=(k[7]+k[8]*active**2/(1+k[9]*active**2))/k[10]
            state=np.array([active,q,actin]);f,j=field_jac(state,p)
            answer.append(dict(state=state,eigenvalues=np.linalg.eigvals(j),residual=float(max(abs(f)))))
    return answer
