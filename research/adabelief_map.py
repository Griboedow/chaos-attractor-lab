"""Scalar quadratic AdaBelief, with explicit epsilon placement and state scaling.

Author implementation: Adabelief-Optimizer 2855178eff6776762793f915ce2c615f776a1e4b,
pypi_packages/adabelief_pytorch0.2.1/adabelief_pytorch/AdaBelief.py.
This is an independent NumPy implementation of the nonrectified, non-AMSGrad,
zero-weight-decay branch. The autonomous map is its bias-correction limit.
For f(theta)=theta**2/2 use x=theta/alpha, m=moment/alpha, s=variance/alpha**2.
"""
import numpy as np


def step(state, b, c, eps_s=1e-12, eps_out=1e-14, iteration=None):
    """State has last axis (x,m,s); iteration is 1-based, None is autonomous."""
    x, m, s = np.moveaxis(np.asarray(state), -1, 0)
    mn = b*m + (1-b)*x
    residual = x-mn
    sn = c*s + (1-c)*residual**2 + eps_s
    correction1 = 1 if iteration is None else 1-b**iteration
    correction2 = 1 if iteration is None else 1-c**iteration
    xn = x - (mn/correction1)/(np.sqrt(sn/correction2)+eps_out)
    return np.stack((xn, mn, sn), axis=-1)


def jacobian(state, b, c, eps_s=1e-12, eps_out=1e-14):
    """Analytic derivative of the autonomous map; supports batched states."""
    x, m, s = np.moveaxis(np.asarray(state), -1, 0)
    mn = b*m + (1-b)*x
    residual = x-mn
    sn = c*s + (1-c)*residual**2 + eps_s
    root = np.sqrt(sn)
    denom = root+eps_out
    j = np.zeros(np.shape(x)+(3,3), dtype=np.result_type(state, float))
    j[...,1,0] = 1-b
    j[...,1,1] = b
    j[...,2,0] = 2*(1-c)*b*residual
    j[...,2,1] = -j[...,2,0]
    j[...,2,2] = c
    j[...,0,:] = -j[...,1,:]/denom[...,None]
    j[...,0,:] += (mn/(2*root*denom**2))[...,None]*j[...,2,:]
    j[...,0,0] += 1
    return j


def inverse_step(state, b, c, eps_s=1e-12, eps_out=1e-14):
    """Inverse on the image of the positive-variance domain."""
    xn, mn, sn = np.moveaxis(np.asarray(state), -1, 0)
    x = xn + mn/(np.sqrt(sn)+eps_out)
    m = (mn-(1-b)*x)/b
    s = (sn-(1-c)*(x-mn)**2-eps_s)/c
    return np.stack((x,m,s),axis=-1)


def spectrum(states, b, c, *, warmup=10000, steps=100000, window=10000,
             eps_s=1e-12, eps_out=1e-14, keep=20000, bias_warmup=False):
    """Every-step QR, finite-time windows, unthinned tail, all-state extrema.

    Tangents align during warmup. Bias-corrected warmup aligns to the limiting
    Jacobian; only subsequent autonomous iterations enter reported spectra.
    """
    if steps <= 0 or steps % window or not 0 <= keep <= steps:
        raise ValueError('steps must be a positive multiple of window; 0<=keep<=steps')
    z = np.atleast_2d(np.array(states, dtype=float))
    b, c = np.broadcast_arrays(np.asarray(b), np.asarray(c))
    q = np.broadcast_to(np.eye(3), (len(z),3,3)).copy()
    acc = np.zeros((len(z),3)); windows=[]; tail=[]
    lo = np.full_like(z, np.inf); hi=-lo
    for t in range(warmup+steps):
        a = jacobian(z,b,c,eps_s,eps_out) @ q
        q,r = np.linalg.qr(a)
        z = step(z,b,c,eps_s,eps_out,
                 t+1 if bias_warmup and t < warmup else None)
        if not np.isfinite(z).all() or (z[:,2] <= 0).any():
            raise FloatingPointError(f'invalid state at iteration {t+1}')
        if t >= warmup:
            acc += np.log(np.abs(np.diagonal(r,axis1=-2,axis2=-1)))
            lo=np.minimum(lo,z); hi=np.maximum(hi,z)
            if (t-warmup+1)%window == 0:
                windows.append(acc/window); acc[:]=0
            if t >= warmup+steps-keep:
                tail.append(z.copy())
    return dict(end=z, windows=np.array(windows), spectrum=np.mean(windows,axis=0),
                tail=np.array(tail), minimum=lo, maximum=hi,
                determinant_sum=np.log(b*c))
