"""Outward binary64 interval jets for the spatial-circle RAAR formula.

Every elementary arithmetic result is expanded by nextafter. No fastmath or
parallel reductions are used. Square roots require nonnegative inputs and use
the platform binary64 sqrt, expanded by one representable number. Parameters
are exact binary64 lifts; sin(theta), cos(theta) must be supplied as outward
intervals generated separately with arbitrary-precision interval arithmetic.
This module alone proves no trapping or attraction statement.
"""
import numpy as np
from numba import njit


@njit(cache=True, inline="always")
def down(x):
    return np.nextafter(x, -np.inf)


@njit(cache=True, inline="always")
def up(x):
    return np.nextafter(x, np.inf)


@njit(cache=True, inline="always")
def add(a, b):
    return np.array([down(a[0] + b[0]), up(a[1] + b[1])])


@njit(cache=True, inline="always")
def neg(a):
    return np.array([-a[1], -a[0]])


@njit(cache=True, inline="always")
def sub(a, b):
    return add(a, neg(b))


@njit(cache=True, inline="always")
def mul(a, b):
    x0, x1, x2, x3 = a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]
    return np.array([down(min(x0, x1, x2, x3)), up(max(x0, x1, x2, x3))])


@njit(cache=True, inline="always")
def square(a):
    if a[0] <= 0 <= a[1]:
        return np.array([0., up(max(a[0]*a[0], a[1]*a[1]))])
    return np.array([down(min(a[0]*a[0], a[1]*a[1])), up(max(a[0]*a[0], a[1]*a[1]))])


@njit(cache=True, inline="always")
def reciprocal(a):
    assert a[0] > 0 or a[1] < 0
    return np.array([down(1./a[1]), up(1./a[0])])


@njit(cache=True, inline="always")
def root(a):
    assert a[0] >= 0
    return np.array([max(0., down(np.sqrt(a[0]))), up(np.sqrt(a[1]))])


@njit(cache=True, inline="always")
def constant(a):
    j = np.zeros((4, 2))
    j[0] = a
    return j


@njit(cache=True, inline="always")
def scalar(x):
    return constant(np.array([x, x]))


@njit(cache=True, inline="always")
def jadd(a, b):
    out = np.empty((4, 2))
    for k in range(4):
        out[k] = add(a[k], b[k])
    return out


@njit(cache=True, inline="always")
def jneg(a):
    out = np.empty((4, 2))
    for k in range(4):
        out[k] = neg(a[k])
    return out


@njit(cache=True, inline="always")
def jsub(a, b):
    return jadd(a, jneg(b))


@njit(cache=True, inline="always")
def jmul(a, b):
    out = np.empty((4, 2))
    out[0] = mul(a[0], b[0])
    for k in range(1, 4):
        out[k] = add(mul(a[k], b[0]), mul(a[0], b[k]))
    return out


@njit(cache=True, inline="always")
def jsquare(a):
    out = np.empty((4, 2))
    out[0] = square(a[0])
    two_a = mul(np.array([2., 2.]), a[0])
    for k in range(1, 4):
        out[k] = mul(two_a, a[k])
    return out


@njit(cache=True, inline="always")
def jreciprocal(a):
    out = np.empty((4, 2))
    out[0] = reciprocal(a[0])
    derivative = neg(square(out[0]))
    for k in range(1, 4):
        out[k] = mul(derivative, a[k])
    return out


@njit(cache=True, inline="always")
def jroot(a):
    out = np.empty((4, 2))
    out[0] = root(a[0])
    derivative = reciprocal(mul(np.array([2., 2.]), out[0]))
    for k in range(1, 4):
        out[k] = mul(derivative, a[k])
    return out


@njit(cache=True)
def evaluate_box(box, parameters, trig):
    """Return F(box), DF(box), lower axis distance and failure flag.

    box has shape (3,2), trig is [[sin lower, sin upper],[cos lower,cos upper]].
    A failure means this extension could not establish regularity, not that
    the box necessarily contains a true singularity.
    """
    x = np.empty((3, 4, 2))
    for i in range(3):
        x[i] = constant(box[i])
        x[i, i+1] = np.array([1., 1.])
    beta, radius = scalar(parameters[0]), scalar(parameters[1])
    sn, co = constant(trig[0]), constant(trig[1])
    c = np.empty((3, 4, 2))
    for i in range(3):
        c[i] = scalar(parameters[i+3])
    dx = jsub(x[0], c[0])
    along = jadd(jmul(co, jsub(x[1], c[1])), jmul(sn, jsub(x[2], c[2])))
    length2 = jadd(jsquare(dx), jsquare(along))
    values, jacobian = np.empty((3, 2)), np.empty((3, 3, 2))
    if length2[0, 0] <= 0:
        return values, jacobian, 0., 1
    length = jroot(length2)
    factor = jmul(radius, jreciprocal(length))
    b = np.empty((3, 4, 2))
    b[0] = jadd(c[0], jmul(factor, dx))
    b[1] = jadd(c[1], jmul(factor, jmul(co, along)))
    b[2] = jadd(c[2], jmul(factor, jmul(sn, along)))
    rx = jsub(jmul(scalar(2.), b[0]), x[0])
    ry = jsub(jmul(scalar(2.), b[1]), x[1])
    lengtha2 = jadd(jsquare(rx), jsquare(ry))
    if lengtha2[0, 0] <= 0:
        return values, jacobian, length[0, 0], 2
    lengtha = jroot(lengtha2)
    inversea = jreciprocal(lengtha)
    a = np.empty((3, 4, 2))
    a[0], a[1], a[2] = jmul(rx, inversea), jmul(ry, inversea), scalar(0.)
    coefficient = jsub(scalar(1.), jmul(scalar(2.), beta))
    for i in range(3):
        result = jadd(jadd(jmul(beta, x[i]), jmul(beta, a[i])), jmul(coefficient, b[i]))
        values[i] = result[0]
        jacobian[i] = result[1:]
    return values, jacobian, min(length[0, 0], lengtha[0, 0]), 0


@njit(cache=True)
def centered_remainder(box, center, image_center, linear_part, parameters, trig, margin):
    """Enclose F(box)+[-margin,margin]^3 by a selected affine part plus a box.

    image_center and linear_part are exact chosen binary64 coefficients, not
    asserted exact function values. The interval calculation encloses their
    errors. The integral mean-value formula gives the stated enclosure.
    """
    _, jacobian, axis, status = evaluate_box(box, parameters, trig)
    point = np.column_stack((center, center))
    values, _, center_axis, center_status = evaluate_box(point, parameters, trig)
    if status or center_status:
        return np.zeros(3), np.zeros(3), min(axis, center_axis), status or center_status
    radii = np.empty(3)
    for i in range(3):
        radii[i] = up(max(up(center[i]-box[i, 0]), up(box[i, 1]-center[i])))
    remainder = np.empty(3)
    for i in range(3):
        mismatch = sub(values[i], np.array([image_center[i], image_center[i]]))
        bound = up(max(abs(mismatch[0]), abs(mismatch[1])) + margin)
        for j in range(3):
            difference = sub(jacobian[i, j], np.array([linear_part[i, j], linear_part[i, j]]))
            error = max(abs(difference[0]), abs(difference[1]))
            bound = up(bound + up(error*radii[j]))
        remainder[i] = bound
    return radii, remainder, min(axis, center_axis), 0
