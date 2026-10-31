"""Planar bead-drag follower-force chain, a declared extension of De Canio 2017.

Angles are unwrapped elastic coordinates, not reduced modulo 2*pi. There is no
contact law or long-range hydrodynamic interaction. The N=2, unit-parameter case
recovers Eqs. (30)-(31) of arXiv:1707.02346 (the second equation has a missing
minus sign in the preprint typesetting; virtual work and Eq. 34 fix the sign).
"""
import numpy as np
from numba import njit


def matrices(theta, force, stiffness, drag, length):
    """Independent general-N virtual-work implementation (also complex-safe)."""
    n = len(theta)
    tail = np.cumsum(np.asarray(drag)[::-1])[::-1]
    base = tail[np.maximum.outer(np.arange(n), np.arange(n))]
    base = base * np.outer(length, length)
    resistance = base * np.cos(theta[:, None] - theta[None, :])
    difference = np.eye(n) - np.eye(n, k=-1)
    elastic = difference.T @ np.diag(stiffness) @ difference
    torque = force * length * np.sin(theta - theta[-1]) - elastic @ theta
    return resistance, torque


def field(theta, force, stiffness, drag, length):
    m, g = matrices(theta, force, stiffness, drag, length)
    return np.linalg.solve(m, g)


@njit(cache=True)
def solve3(m, v):
    a, b, c = m[0, 0], m[0, 1], m[0, 2]
    d, e, f = m[1, 1], m[1, 2], m[2, 2]
    aa, bb, cc = d*f-e*e, c*e-b*f, b*e-c*d
    dd, ee, ff = a*f-c*c, b*c-a*e, a*d-b*b
    det = a*aa+b*bb+c*cc
    return np.array([aa*v[0]+bb*v[1]+cc*v[2],
                     bb*v[0]+dd*v[1]+ee*v[2],
                     cc*v[0]+ee*v[1]+ff*v[2]])/det


@njit(cache=True)
def augmented3(z, force, elastic, base, length):
    x, v = z[:3], z[3:]
    m = np.empty((3, 3))
    dm = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            d = x[i]-x[j]
            m[i, j] = base[i, j]*np.cos(d)
            dm[i, j] = -base[i, j]*np.sin(d)*(v[i]-v[j])
    g = force*length*np.sin(x-x[2])-elastic@x
    f = solve3(m, g)
    dg = force*length*np.cos(x-x[2])*(v-v[2])-elastic@v
    dv = solve3(m, dg-dm@f)
    result = np.empty(6)
    result[:3], result[3:] = f, dv
    return result


@njit(cache=True)
def pilot(force, stiffness, drag, length, initial, step, warm, duration, stride):
    """Fixed-step screen only; candidates require independent adaptive checks."""
    difference = np.eye(3)-np.eye(3, k=-1)
    elastic = difference.T@np.diag(stiffness)@difference
    tail = np.cumsum(drag[::-1])[::-1]
    base = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            base[i, j] = tail[max(i, j)]*length[i]*length[j]
    nw, no = int(round(warm/step)), int(round(duration/step))
    trace = np.empty((no//stride, 3))
    z = np.empty(6)
    z[:3], z[3:] = initial, np.array([1., 0., 0.])
    logs, blocklogs, min_norm, max_angle = 0., np.zeros(4), 1e300, 0.
    qr = max(1, int(round(.1/step)))
    assert nw % qr == 0 and no % (4*qr) == 0 and no % stride == 0
    for i in range(nw+no):
        a = augmented3(z, force, elastic, base, length)
        b = augmented3(z+.5*step*a, force, elastic, base, length)
        c = augmented3(z+.5*step*b, force, elastic, base, length)
        d = augmented3(z+step*c, force, elastic, base, length)
        z += (step/6)*(a+2*b+2*c+d)
        if not np.isfinite(z).all():
            return trace, np.nan, blocklogs, z, np.nan, np.nan
        if (i+1) % qr == 0:
            norm = np.linalg.norm(z[3:])
            z[3:] /= norm
            if i >= nw:
                inc = np.log(norm)
                logs += inc
                blocklogs[min(3, (i-nw)//(no//4))] += inc
        if i >= nw:
            min_norm = min(min_norm, np.linalg.norm(z[:3]))
            max_angle = max(max_angle, np.max(np.abs(z[:3])))
            if (i-nw+1) % stride == 0:
                trace[(i-nw)//stride] = z[:3]
    return trace, logs/duration, blocklogs/(duration/4), z, min_norm, max_angle
