"""Outward binary64 scalar interval primitives for the spatial-circle RAAR formula.

Every elementary arithmetic result is expanded by nextafter. No fastmath or
parallel reductions are used. Square roots require nonnegative inputs and use
the platform binary64 sqrt, expanded by one representable number. Parameters
are exact binary64 lifts; sin(theta), cos(theta) must be supplied as outward
intervals generated separately with arbitrary-precision interval arithmetic.
This module alone proves no trapping or attraction statement.
"""
import numpy as np
from numba import njit


@njit(cache=False, inline="always")
def down(x):
    return np.nextafter(x, -np.inf)


@njit(cache=False, inline="always")
def up(x):
    return np.nextafter(x, np.inf)


@njit(cache=False, inline="always")
def add(a, b):
    return (down(a[0] + b[0]), up(a[1] + b[1]))


@njit(cache=False, inline="always")
def neg(a):
    return (-a[1], -a[0])


@njit(cache=False, inline="always")
def sub(a, b):
    return add(a, neg(b))


@njit(cache=False, inline="always")
def mul(a, b):
    x0, x1, x2, x3 = a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]
    return (down(min(x0, x1, x2, x3)), up(max(x0, x1, x2, x3)))


@njit(cache=False, inline="always")
def square(a):
    if a[0] <= 0 <= a[1]:
        return (0., up(max(a[0]*a[0], a[1]*a[1])))
    return (down(min(a[0]*a[0], a[1]*a[1])), up(max(a[0]*a[0], a[1]*a[1])))


@njit(cache=False, inline="always")
def reciprocal(a):
    assert a[0] > 0 or a[1] < 0
    return (down(1./a[1]), up(1./a[0]))


@njit(cache=False, inline="always")
def root(a):
    assert a[0] >= 0
    return (max(0., down(np.sqrt(a[0]))), up(np.sqrt(a[1])))


