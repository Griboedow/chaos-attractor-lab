"""Dimensionless delayed-stretch-activation model of Lynch et al. (2022).

Source equations (4), (9), not an invented oscillator or a chaos claim.
tau = sqrt(k/I)*t, theta_scale=I/Gamma, h=mu*f/(k*theta_scale).
State (theta/theta_scale, dtheta/dtau/theta_scale, h, dh/dtau).
Parameters (r3/sqrt(k/I), kappa=r4/r3, mu/k, time_scale).
The last parameter is a computational clock factor, normally one.
"""
import numpy as np

SOURCE='https://sponberg.gatech.edu/wp-content/uploads/2022/08/Gau-et-al-2022_02.pdf'


def coefficients(p):
    r,ratio,gain,scale=p
    return r*(1+ratio),ratio*r*r,gain*r*(1-ratio),scale


def field_jac(q,p):
    x,v,h,w=q;a,b,c,scale=coefficients(p)
    f=np.array([v,-x-v*abs(v)+h,w,c*v-b*h-a*w])
    j=np.array([[0,1,0,0],[-1,-2*abs(v),1,0],[0,0,0,1],[0,c,-b,-a]])
    return scale*f,scale*j


def dimensional_field(q,k,inertia,drag,gain,r3,ratio):
    theta,velocity,force,force_rate=q
    return np.array([velocity,(-k*theta-drag*velocity*abs(velocity)+gain*force)/inertia,
        force_rate,r3*(1-ratio)*velocity-ratio*r3*r3*force-r3*(1+ratio)*force_rate])


def source_equivalence(q,p,k=.416,inertia=.00196,drag=.003):
    r,ratio,gain,scale=p;omega=np.sqrt(k/inertia);theta_scale=inertia/drag;mu=gain*k
    transform=np.array([theta_scale,theta_scale*omega,k*theta_scale/mu,k*theta_scale*omega/mu])
    original=dimensional_field(q*transform,k,inertia,drag,mu,r*omega,ratio)
    return scale*original/(omega*transform)


def equilibrium_polynomial(p):
    a,b,c,_=coefficients(p)
    return np.array([1,a,1+b,a-c,b])
