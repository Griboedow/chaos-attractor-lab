"""Exact binary Cartesian posterior and its spherical reflected iteration.

Input floats denote exact dyadic real numbers. Half-sums and half-differences
are formed after lifting each stored endpoint into the supplied context.
"""
from fractions import Fraction
import itertools
import numpy as np


def exact(value, ctx):
    n, d = float(value).as_integer_ratio()
    return ctx.mpf(n)/ctx.mpf(d)


def parameters(raw, points, ctx):
    levels = [np.unique(points[:, j]) for j in range(3)]
    assert all(len(v) == 2 for v in levels)
    assert sorted(map(tuple, points)) == sorted(itertools.product(*levels))
    pairs = [[exact(v, ctx) for v in row] for row in levels]
    beta, tau = [exact(v, ctx) for v in raw]
    return beta, tau, [(lo+hi)/2 for lo, hi in pairs], [(hi-lo)/2 for lo, hi in pairs]


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def mv(a, b):
    return [sum(x*y for x, y in zip(row, b)) for row in a]


def inverse(a):
    a = [list(row)+[Fraction(int(i == j)) for j in range(3)] for i, row in enumerate(a)]
    for j in range(3):
        pivot = next(i for i in range(j, 3) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        value = a[j][j]
        a[j] = [v/value for v in a[j]]
        for i in range(3):
            if i != j:
                value = a[i][j]
                a[i] = [a[i][k]-value*a[j][k] for k in range(6)]
    return [row[3:] for row in a]


def model(x, par, ctx):
    beta, tau, c, a = par
    m, d = [], []
    for j in range(3):
        k = 2*a[j]/tau
        w = ctx.exp(2*k*(x[j]-c[j]))
        # Equivalent monotone logistic form avoids artificial width in the
        # saturated branches of (w-1)/(w+1).
        response = 1-2/(w+1) if w >= 1 else -1+2*w/(w+1) if w <= 1 else (w-1)/(w+1)
        m.append(c[j]+a[j]*response)
        d.append(4*a[j]*k*w/(w+1)**2)
    v = [2*s-t for s, t in zip(m, x)]
    r = ctx.sqrt(sum(t**2 for t in v))
    u = [t/r for t in v]
    y = [beta*t+beta*s+(1-2*beta)*q for t, s, q in zip(x, u, m)]
    jac = [[beta*(i == j)+(1-2*beta)*d[i]*(i == j)
            +beta/r*((i == j)-u[i]*u[j])*(2*d[j]-1)
            for j in range(3)] for i in range(3)]
    return y, jac, r


def fourth(x, par, ctx, step=model):
    jac = [[ctx.mpf(i == j) for j in range(3)] for i in range(3)]
    radii = []
    for _ in range(4):
        x, j, r = step(x, par, ctx)
        jac = mm(j, jac)
        radii.append(r)
    return x, jac, radii
