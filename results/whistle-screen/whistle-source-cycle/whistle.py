"""Literal six-state empirical whistle equations of da Silva et al. (2026).

Eq2 pressure sign is retained as printed. Ambient numerical constants are
declared assumptions at the source's20C/1atm, not recovered fitted values.
"""
import numpy as np
from numba import njit

SOURCE='https://doi.org/10.1016/j.jsv.2026.119952'


def reference(name='B'):
    length,volume,pressure,velocity,frequency,helmholtz={
        'A':(.038,80.4e-6,123.2,9.99,505.1,463.3),
        'B':(.0635,134e-6,93.,8.16,392.9,358.1),
        'C':(.0889,188e-6,70.,7.10,331.6,302.)}[name]
    return dict(name=name,length=length,volume=volume,pressure=pressure,velocity=velocity,
        observed_frequency=frequency,reported_helmholtz=helmholtz,diameter=.00635,thickness=.0127,
        density=1.2041,sound_speed=343.2,kinematic_viscosity=1.535e-5,
        end_correction=.72,discharge=.98,resistance=1.,Cv=.0013,alpha=.0127,tau0=.005,gamma=.09,mu=3.1,kappa=1.,beta=.3)


def normalized(p):
    area=np.pi*p['diameter']**2/4;le=p['thickness']/2+p['end_correction']*p['diameter']
    mass=p['density']*le/area;K=p['density']*p['sound_speed']**2/p['volume']
    loss=p['density']/(2*p['discharge']**2*area**2)
    omega=np.sqrt(2*K/mass);Q=area*p['velocity'];P=p['pressure'];G=p['alpha']*p['velocity']
    rem=p['velocity']*np.sqrt(p['length']*p['diameter'])/p['kinematic_viscosity']
    st=697.4*rem**(-.832);wv=2*np.pi*p['velocity']/p['diameter']*st
    coefficients=np.array([P/(mass*Q*omega),p['resistance']/(mass*omega),loss*Q/(mass*omega),
        K*Q/(P*omega),K*p['Cv']*G/(P*omega),p['gamma']*P,p['beta']*p['velocity'],
        1/(p['tau0']*omega),p['mu']/omega,(wv/omega)**2,p['kappa']*P/omega**2])
    scale=np.array([Q,Q,P,G,1.,omega])
    return coefficients,scale,omega,dict(modified_Reynolds=rem,Strouhal=st,vortex_frequency=wv/(2*np.pi),
        computed_helmholtz=omega/(2*np.pi),area=area,effective_length=le)


@njit(cache=True)
def field_jac(state,p):
    x,y,pressure,g,q,v=state
    A,R,B,C,D,G,H,J,U,W,F=p
    e=np.exp(G*abs(pressure));decay=J*e*(1+H*abs(x))
    f=np.array([A*(1-pressure)-R*x-B*x*abs(x),A*pressure-R*y-B*y*abs(y),
        C*(x-y)-D*g,v-decay*g,v,U*(1-q*q)*v-W*q+F*pressure])
    jac=np.zeros((6,6))
    jac[0,0]=-R-2*B*abs(x);jac[0,2]=-A
    jac[1,1]=-R-2*B*abs(y);jac[1,2]=A
    jac[2,0]=C;jac[2,1]=-C;jac[2,3]=-D
    jac[3,0]=-g*J*e*H*np.sign(x);jac[3,2]=-g*decay*G*np.sign(pressure);jac[3,3]=-decay;jac[3,5]=1
    jac[4,5]=1;jac[5,2]=F;jac[5,4]=-2*U*q*v-W;jac[5,5]=U*(1-q*q)
    return f,jac


def dimensional_field(state,p):
    Q1,Q2,pc,circulation,q,velocity=state
    area=np.pi*p['diameter']**2/4;le=p['thickness']/2+p['end_correction']*p['diameter']
    inertia=p['density']*le/area;compliance=p['density']*p['sound_speed']**2/p['volume']
    loss=p['density']/(2*p['discharge']**2*area**2)
    rem=p['velocity']*np.sqrt(p['length']*p['diameter'])/p['kinematic_viscosity']
    frequency=2*np.pi*p['velocity']/p['diameter']*697.4*rem**(-.832)
    inverse_lifetime=np.exp(p['gamma']*abs(pc))*(1+p['beta']*abs(Q1/area))/p['tau0']
    return np.array([(p['pressure']-pc-p['resistance']*Q1-loss*Q1*abs(Q1))/inertia,
        (pc-p['resistance']*Q2-loss*Q2*abs(Q2))/inertia,
        -compliance*(Q2-Q1+p['Cv']*circulation),p['alpha']*p['velocity']*velocity-inverse_lifetime*circulation,
        velocity,p['mu']*(1-q*q)*velocity-frequency**2*q+p['kappa']*pc])
