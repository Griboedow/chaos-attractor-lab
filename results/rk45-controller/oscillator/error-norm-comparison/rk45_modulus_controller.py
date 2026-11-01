"""SciPy-style RK45 accepted-step amplitude map on stable complex linear ODEs.

The reduction removes the S1 phase symmetry; any added phase exponent is zero.
It transcribes the installed SciPy acceptance/factor rule, omitting finite-time
endpoint clipping and floating-point t+h subtraction. An upstream check is
required before describing a candidate as reproduced by SciPy itself.
"""
import argparse,hashlib,inspect,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.integrate._ivp.rk import RK45,rk_step,RungeKutta
from scipy.stats import qmc
import scipy

@njit(cache=True)
def rk45_step(state,params):
    # state=(log proposed h, log |y|/atol); params=(angle,rtol).
    x,y=state;angle,rtol=params
    dx=np.array([1.,0.]);dy=np.array([0.,1.])
    lam=complex(-np.cos(angle),np.sin(angle))
    rejected=False
    for attempt in range(100):
        h=np.exp(x);z=h*lam
        R=1+z+z*z/2+z**3/6+z**4/24+z**5/120+z**6/600
        E=z**5*(5*z*z-39*z+97)/120000
        modulus=abs(R)
        if modulus==0:
            # Exact arrival at y=0 leaves this positive-amplitude chart.
            return state,np.eye(2),0.,attempt,2
        mu=(z*(1+z+z*z/2+z**3/6+z**4/24+z**5/100)/R).real
        scale_term=rtol*np.exp(y)*max(1.,modulus)
        dscale=dy.copy()
        if modulus>1:dscale+=mu*dx
        if E==0:
            # The capped-factor branch is constant in a neighborhood of this
            # zero. Do not form E'/E before selecting that branch.
            logerror=-np.inf;de=np.zeros(2)
        else:
            nu=(5+z*(10*z-39)/(5*z*z-39*z+97)).real
            logerror=np.log(abs(E))+y-np.log1p(scale_term)
            de=nu*dx+dy-scale_term/(1+scale_term)*dscale
        logfactor=np.log(.9)-logerror/5
        df=-de/5
        if logerror<0:
            ceiling=0. if rejected else np.log(10.)
            if logfactor>ceiling:logfactor=ceiling;df*=0
            result=np.array([x+logfactor,y+np.log(modulus)])
            jac=np.empty((2,2));jac[0]=dx+df;jac[1]=dy+mu*dx
            return result,jac,np.angle(R),attempt,0
        if logfactor<np.log(.2):logfactor=np.log(.2);df*=0
        x+=logfactor;dx+=df;rejected=True
    return state,np.eye(2),0.,100,1

def verify():
    rng=np.random.default_rng(2026110172);errors=[];polynomial_errors=[]
    for _ in range(512):
        angle=rng.uniform(0,np.pi/2-.001);rtol=10**rng.uniform(-12,-1);p=np.array([angle,rtol]);s=rng.uniform(-1.5,1.5,2)
        _,j,_,_,status=rk45_step(s,p);assert status==0
        finite=np.column_stack([(rk45_step(s+np.eye(2)[i]*1e-6,p)[0]-rk45_step(s-np.eye(2)[i]*1e-6,p)[0])/2e-6 for i in range(2)])
        errors.append(float(np.max(abs(j-finite))))
        lam=-np.cos(angle)+1j*np.sin(angle);h=np.exp(s[0]);z=h*lam;y=np.array([np.exp(s[1])+0j]);K=np.zeros((7,1),dtype=complex)
        new,_=rk_step(lambda t,u:lam*u,0.,y,lam*y,h,RK45.A,RK45.B,RK45.C,K)
        expected=(1+z+z*z/2+z**3/6+z**4/24+z**5/120+z**6/600)*y
        polynomial_errors.append(float(max(abs(new-expected).max(),abs(h*(K.T@RK45.E)-z**5*(5*z*z-39*z+97)/120000*y).max())))
    assert max(errors)<1e-6,max(errors)
    assert max(polynomial_errors)<1e-12,max(polynomial_errors)
    return dict(finite_difference_count=512,max_jacobian_error=max(errors),max_installed_rk_stage_polynomial_error=max(polynomial_errors))
