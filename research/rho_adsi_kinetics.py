"""The same ADSI source equations with independently variable source constants.

The first eleven entries are k0..k10, with GEF and GAP folded into k1 and k4.
Five padding entries are ignored by the GPU double16 representation. No extra
reactions, clipping, diffusion, noise or forcing are introduced.
"""
import numpy as np
from decimal import Decimal, localcontext
from rho_adsi import BASE


def field_jac(state,k):
    x,y,z=state
    den=1+k[2]*x*x
    a=k[0]+k[1]*x**3/den
    da=k[1]*x*x*(3+k[2]*x*x)/den**2
    b=k[3]+k[4]*z;rx=a*y-b*x
    field=np.array([rx,-rx+k[5]-k[6]*y,k[7]+k[8]*x*x/(1+k[9]*x*x)-k[10]*z])/.025
    jac=np.array([[da*y-b,a,-k[4]*x],[-da*y+b,-a-k[6],k[4]*x],
                  [2*k[8]*x/(1+k[9]*x*x)**2,0.,-k[10]]])/.025
    return field,jac


def independent_field(state,k):
    active,inactive,actin=state;total=active+inactive
    active_dot=(k[0]+k[1]*active**3/(1+k[2]*active**2))*(total-active)-(k[3]+k[4]*actin)*active
    total_dot=k[5]-k[6]*(total-active)
    actin_dot=k[7]+k[8]*active**2/(1+k[9]*active**2)-k[10]*actin
    return np.array([active_dot,total_dot-active_dot,actin_dot])/.025


def decimal_source_check(state,k):
    """70-digit source-form field and absolute reaction-flux scale.

    The total-minus-active source coordinates can lose digits when active Rho
    is large. Normalize differences by the physical terms being subtracted,
    while retaining the raw residual separately. Decimal uses the exact input
    binary floats; no extra integration precision is claimed.
    """
    with localcontext() as ctx:
        ctx.prec=70
        x,y,z=[Decimal.from_float(float(v)) for v in state]
        q=[Decimal.from_float(float(v)) for v in k]
        one=Decimal(1);tau=Decimal.from_float(.025);total=x+y
        activate=(q[0]+q[1]*x**3/(one+q[2]*x*x))*(total-x)
        inactivate=(q[3]+q[4]*z)*x
        inflow=q[5];outflow=q[6]*(total-x)
        basal=q[7];stimulus=q[8]*x*x/(one+q[9]*x*x);decay=q[10]*z
        f=[activate-inactivate,inflow-outflow-activate+inactivate,basal+stimulus-decay]
        scale=[abs(activate)+abs(inactivate),abs(inflow)+abs(outflow)+abs(activate)+abs(inactivate),abs(basal)+abs(stimulus)+abs(decay)]
        return np.array([float(v/tau) for v in f]),np.array([float(v/tau) for v in scale])


def equilibria(k):
    poly=np.polynomial.Polynomial
    a=poly([1.,0.,k[2]]);b=poly([1.,0.,k[9]]);x=poly([0.,1.]);q=k[5]/k[6]
    numerator=q*(k[0]*a+k[1]*x**3)*b-x*(k[3]*a*b+k[4]/k[10]*a*(k[7]*b+k[8]*x*x))
    answer=[]
    for root in numerator.roots():
        if abs(root.imag)<1e-8 and root.real>0:
            active=root.real;actin=(k[7]+k[8]*active**2/(1+k[9]*active**2))/k[10]
            state=np.array([active,q,actin]);f,j=field_jac(state,k)
            answer.append(dict(state=state,eigenvalues=np.linalg.eigvals(j),residual=float(max(abs(f)))))
    return answer


def parameterize(rows):
    """GEF,GAP/GEF,inhibitor speed,exchange speed,k0 factor,pool factor,k2 factor,k9 factor."""
    k=np.zeros((len(rows),16));k[:,:11]=BASE
    k[:,1]*=rows[:,0];k[:,4]*=rows[:,0]*rows[:,1]
    k[:,[7,8,10]]*=rows[:,2,None];k[:,[5,6]]*=rows[:,3,None]
    k[:,0]*=rows[:,4];k[:,5]*=rows[:,5];k[:,2]*=rows[:,6];k[:,9]*=rows[:,7]
    return k
