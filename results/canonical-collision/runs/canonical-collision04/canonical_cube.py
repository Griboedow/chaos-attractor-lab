"""Exact factorization of the stored equal-center cubic posterior.

The endpoints, beta, and tau are exact binary64 real inputs. Their half-sum
and half-difference are formed in the supplied arbitrary-precision context.
This module does not identify rounded table parameters with those endpoints.
"""
import itertools
import numpy as np


def exact(value, ctx):
    n, d = float(value).as_integer_ratio()
    return ctx.mpf(n)/ctx.mpf(d)


def parameters(source, ident, ctx):
    with np.load(source) as f:
        raw = f['params'][ident].copy()
        points = f['points'][ident, :f['counts'][ident]].copy()
    assert len(points) == 8
    endpoints = np.unique(points)
    assert len(endpoints) == 2
    expected = np.array(list(itertools.product(endpoints, repeat=3)))
    assert sorted(map(tuple, points)) == sorted(map(tuple, expected))
    low, high = [exact(v, ctx) for v in endpoints]
    beta, tau = [exact(v, ctx) for v in raw]
    return beta, tau, (low+high)/2, (high-low)/2


def mean_one(x, par, ctx):
    beta, tau, c, a = par
    k = 2*a/tau
    w = ctx.exp(2*k*(x-c))
    return c+a*(w-1)/(w+1), 4*a*k*w/(w+1)**2


def model(x, par, ctx, derivative=True):
    beta, tau, c, a = par
    pairs = [mean_one(t, par, ctx) for t in x]
    m = [p[0] for p in pairs]; d = [p[1] for p in pairs]
    v = [2*s-t for s, t in zip(m, x)]
    r = ctx.sqrt(sum(t*t for t in v)); u = [t/r for t in v]
    y = [beta*t+beta*s+(1-2*beta)*q for t, s, q in zip(x, u, m)]
    if not derivative: return y
    J = [[beta*(i==j)+(1-2*beta)*d[i]*(i==j)
          +beta/r*((i==j)-u[i]*u[j])*(2*d[j]-1)
          for j in range(3)] for i in range(3)]
    return y, J, r
