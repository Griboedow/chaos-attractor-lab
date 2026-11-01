"""Geometric second-order remainder for RAAR's two circle normalizations.

For n(v)=v/|v|, write w=a*u+b with b perpendicular to u=v/|v|.
Then |D^2 n(v)[w,w]|^2=(4*a^2*|b|^2+|b|^4)/|v|^4,
which is at most (4/3)*|w|^4/|v|^4. Taylor's integral remainder
is therefore at most (1/sqrt(3))*|w|^2/s_min^2. We use the larger
exact rational coefficient 3/5. Norms and matrix bounds below are outward.
The smoothness assumptions are checked by positive lower axis distances.
"""
import numpy as np
from numba import njit
from raar_interval64 import add, sub, mul, square, reciprocal, root, up, down, evaluate_box


@njit(cache=True, inline="always")
def point(x):
    return np.array([x, x])


@njit(cache=True)
def norm_upper(v):
    total = point(0.)
    for i in range(len(v)):
        total = add(total, square(v[i]))
    # Added exact zeros can introduce a negative subnormal lower endpoint.
    total[0] = max(0., total[0])
    return root(total)[1]


@njit(cache=True)
def remainder(box, center, image_center, linear_part, p, trig, margin):
    radii = np.empty(3)
    for i in range(3):
        radii[i] = up(max(up(center[i]-box[i, 0]), up(box[i, 1]-center[i])))
    sn, co = trig[0], trig[1]
    dx = sub(box[0], point(p[3]))
    along = add(mul(co, sub(box[1], point(p[4]))), mul(sn, sub(box[2], point(p[5]))))
    q2 = add(square(dx), square(along))
    if q2[0] <= 0:
        return radii, np.zeros(3), 0., 1
    qmin = root(q2)[0]
    # Maximum projected displacement over the symmetric source box.
    projected2 = 0.
    for k in range(8):
        delta = radii*np.array([2*(k&1)-1., 2*((k>>1)&1)-1., 2*((k>>2)&1)-1.])
        along_delta = add(mul(co, point(delta[1])), mul(sn, point(delta[2])))
        projected2 = max(projected2, add(square(point(delta[0])), square(along_delta))[1])
    kappa = mul(point(3.), reciprocal(point(5.)))
    eb = mul(mul(mul(kappa, point(p[1])), point(projected2)), reciprocal(point(q2[0])))[1]
    # The B projection and its derivative at the selected exact center.
    dx0 = sub(point(center[0]), point(p[3]))
    h0 = add(mul(co, sub(point(center[1]), point(p[4]))), mul(sn, sub(point(center[2]), point(p[5]))))
    q0 = root(add(square(dx0), square(h0)))
    iq0 = reciprocal(q0)
    v = np.empty((3, 2))
    v[0], v[1], v[2] = dx0, mul(co, h0), mul(sn, h0)
    u, b = np.empty((3, 2)), np.empty((3, 2))
    for i in range(3):
        u[i] = mul(v[i], iq0)
        b[i] = add(point(p[i+3]), mul(point(p[1]), u[i]))
    normal = np.array([[0., 0.], [-sn[1], -sn[0]], [co[0], co[1]]])
    db = np.empty((3, 3, 2))
    factor = mul(point(p[1]), iq0)
    for i in range(3):
        for j in range(3):
            projection = sub(point((1. if i == j else 0.)), mul(normal[i], normal[j]))
            db[i, j] = mul(factor, sub(projection, mul(u[i], u[j])))
    a_input = np.empty((2, 2))
    for i in range(2):
        a_input[i] = sub(mul(point(2.), b[i]), point(center[i]))
    a02 = add(square(a_input[0]), square(a_input[1]))
    if a02[0] <= 0:
        return radii, np.zeros(3), qmin, 2
    a0 = root(a02)
    # Bound the change in the reflected argument in the A plane.
    linear_dy = 0.
    for k in range(8):
        delta = radii*np.array([2*(k&1)-1., 2*((k>>1)&1)-1., 2*((k>>2)&1)-1.])
        variation = np.empty((2, 2))
        for i in range(2):
            total = point(0.)
            for j in range(3):
                coefficient = sub(mul(point(2.), db[i, j]), point((1. if i == j else 0.)))
                total = add(total, mul(coefficient, point(delta[j])))
            variation[i] = total
        linear_dy = max(linear_dy, norm_upper(variation))
    dy = up(linear_dy+up(2.*eb))
    amin = down(a0[0]-dy)
    if amin <= 0:
        return radii, np.zeros(3), min(qmin, amin), 2
    ea = mul(mul(kappa, square(point(dy))), reciprocal(square(point(amin))))[1]
    # Rows of M=2*beta*D(P_A)(reflected center)+(1-2*beta)*I.
    ua = np.empty((2, 2))
    for i in range(2):
        ua[i] = mul(a_input[i], reciprocal(a0))
    values, jacobian, _, status = evaluate_box(np.column_stack((center, center)), p, trig)
    if status:
        return radii, np.zeros(3), min(qmin, amin), status
    beta = point(p[0])
    twice_beta = mul(point(2.), beta)
    diagonal_coefficient = sub(point(1.), twice_beta)
    result = np.empty(3)
    for i in range(3):
        row = np.empty((3, 2))
        for j in range(3):
            entry = point(0.)
            if i < 2 and j < 2:
                entry = mul(twice_beta, mul(reciprocal(a0), sub(point((1. if i == j else 0.)), mul(ua[i], ua[j]))))
            if i == j:
                entry = add(entry, diagonal_coefficient)
            row[j] = entry
        bound = up(norm_upper(row)*eb)
        if i < 2:
            bound = up(bound+up(p[0]*ea))
        mismatch = sub(values[i], point(image_center[i]))
        bound = up(bound+max(abs(mismatch[0]), abs(mismatch[1])))
        for j in range(3):
            mismatch = sub(jacobian[i, j], point(linear_part[i, j]))
            bound = up(bound+up(max(abs(mismatch[0]), abs(mismatch[1]))*radii[j]))
        result[i] = up(bound+margin)
    return radii, result, min(qmin, amin), 0
