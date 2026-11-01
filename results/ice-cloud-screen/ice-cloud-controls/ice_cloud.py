"""Literal positive-interior ice-cloud model, Bergner/Spichtinger v1 Eq.55--57.

The authors freeze environmental pressure/temperature while retaining the
updraft source in saturation. This is their autonomous approximation, not a
vertically evolving atmospheric parcel. Nucleation mass is retained here.
"""
import numpy as np
from scipy.optimize import brentq, root

SOURCE='https://arxiv.org/abs/2507.03475v1; final article DOI10.1063/5.0297531'


def coefficients(T=230.,pressure=30000.,updraft=.1,sedimentation=1.,aerosols=3e8,cp=1005.):
    """SI quantities; cp is a declared dry-air convention, not recovered source code."""
    Ra=287.058;Rv=461.553;eps=.622;r0=3.;g=9.81
    latent=(46782.5+35.8925*T-.07414*T*T+541.5*np.exp(-(T/123.75)**2))/.01801528
    psi=np.exp(9.550426-5723.265/T+3.53068*np.log(T)-.00728332*T)
    conductivity=.002646*T**1.5/(T+245.*10**(-12./T))
    diffusivity=2.1422e-5*(T/273.15)**2*101325./pressure
    howell=1./((latent/(Rv*T)-1)*latent/T*diffusivity/conductivity+Rv*T/psi)
    radius_factor=(3./(4*np.pi*810.))**(1./3)
    Bq=4*np.pi*howell*diffusivity*radius_factor*r0**(-1./9)
    Cn=600000.*r0**(-1./9)/100.
    Cq=600000.*r0**(5./9)/100.
    density=pressure/(Ra*T)
    Vsol=4*np.pi/3*(75e-9)**3*np.exp(4.5*np.log(1.5)**2)
    An=aerosols/density*Vsol*1e16
    As=(latent/(cp*Rv*T*T)-1./(Ra*T))*g*updraft
    Bs=pressure/(eps*psi)
    p1=np.log(10)*(-50.4085+.9368*T)
    p2=-1.36989e-5*T*T+.00228*T+1.67469
    return np.array([An,1e-16*An,Bq,sedimentation*Cn,sedimentation*Cq,As,Bs,p1,p2])


def physical_field(x,c):
    """Independent literal dimensional Eq.55--57 evaluation."""
    n,q,s=x;An,Aq,Bq,Cn,Cq,As,Bs,k,p2=c
    if n<=0 or q<=0:raise ValueError('The unregularized model requires n,q>0')
    nuc=np.exp(k*(s-p2))
    growth=Bq*n**(2./3)*q**(1./3)*(s-1)
    return np.array([An*nuc+(n/q)*min(growth,0.)-Cn*n**(1./3)*q**(2./3),
                     Aq*nuc+growth-Cq*n**(-2./3)*q**(5./3),
                     As*s-Bs*(Aq*nuc+growth)])


def make_normalized(c):
    """u=log(n/n0),v=log(q/q0),s; tau=As*t. Scale uses full equilibrium."""
    An,Aq,Bq,Cn,Cq,As,Bs,k,p2=c
    log_b=np.log(Bs*Bq*An/Cn)+.25*np.log(Cq/Bq)
    s0=brentq(lambda s:np.log(As*s)-log_b-.75*np.log(s-1)-k*(s-p2),1+1e-12,3.)
    n0=An/Cn*np.exp(k*(s0-p2))*(Bq/Cq*(s0-1))**(-.5)
    q0=n0*(Bq/Cq*(s0-1))**(.75)
    def build(n,q,s):
        ne=np.exp(k*(s-p2));growth=Bq*(n/q)**(2./3)/As
        return np.array([An*ne/n/As,Aq*ne/q/As,growth,Cn*(q/n)**(2./3)/As,
            Cq*(q/n)**(2./3)/As,1.,Bs*Aq*ne/As,Bs*Bq*n**(2./3)*q**(1./3)/As,k,s])
    preliminary=build(n0,q0,s0)
    sol=root(lambda x:field_jac(x,preliminary)[0],[0.,0.,s0],jac=lambda x:field_jac(x,preliminary)[1],tol=1e-11)
    if not sol.success and np.max(abs(sol.fun))>1e-9:raise RuntimeError(sol.message)
    n0*=np.exp(sol.x[0]);q0*=np.exp(sol.x[1]);s0=sol.x[2]
    p=build(n0,q0,s0)
    assert np.max(abs(field_jac(np.array([0.,0.,s0]),p)[0]))<1e-8
    return p,np.array([n0,q0,1.]),As


def field_jac(z,p):
    u,v,s=z;a,d,e,c,f,g,h,i,k,sref=p
    en=np.exp(k*(s-sref));N=a*en*np.exp(-u);M=d*en*np.exp(-v);H=h*en
    G=e*np.exp(2*(u-v)/3);C=c*np.exp(2*(v-u)/3);F=f*np.exp(2*(v-u)/3)
    I=i*np.exp((2*u+v)/3);r=s-1;rm=min(r,0.)
    field=np.array([N+G*rm-C,M+G*r-F,g*s-H-I*r])
    jac=np.array([[-N+2*(G*rm+C)/3,-2*(G*rm+C)/3,k*N+G*(s<1)],
                  [2*(G*r+F)/3,-M-2*(G*r+F)/3,k*M+G],
                  [-2*I*r/3,-I*r/3,g-k*H-I]])
    return field,jac


def to_physical(z,scale):
    a=np.asarray(z);out=np.array(a,copy=True);out[...,:2]=np.exp(a[...,:2])*scale[:2]
    return out
