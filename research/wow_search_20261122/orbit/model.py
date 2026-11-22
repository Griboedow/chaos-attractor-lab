"""Native Cartesian orbit, analytic full tangent, and fixed-detuning counterfactual."""

import math
import numpy as np
from numba import njit


@njit(cache=True)
def trace(params, initial, burn, window, windows, keep, frozen=np.nan):
    pump, loss, base, fast, gain, rate = params
    x, y, z = initial
    q = np.eye(3)
    v = np.empty((3, 3))
    j = np.empty((3, 3))
    rates = np.zeros((windows, 3))
    cloud = np.empty((keep, 3))
    count = burn + window * windows
    for i in range(count):
        intensity = x * x + y * y
        ri = 1 / (1 + intensity)
        f = intensity * ri
        phase = base + fast * f + (z if np.isnan(frozen) else frozen)
        co = math.cos(phase)
        sn = math.sin(phase)
        u = loss * (co * x - sn * y)
        w = loss * (sn * x + co * y)
        dx = 2 * x * ri * ri
        dy = 2 * y * ri * ri
        j[0, 0] = loss * co - w * fast * dx
        j[0, 1] = -loss * sn - w * fast * dy
        j[0, 2] = -w if np.isnan(frozen) else 0.0
        j[1, 0] = loss * sn + u * fast * dx
        j[1, 1] = loss * co + u * fast * dy
        j[1, 2] = u if np.isnan(frozen) else 0.0
        j[2, 0] = rate * gain * dx
        j[2, 1] = rate * gain * dy
        j[2, 2] = 1.0 - rate
        zz = (1 - rate) * z + rate * gain * f
        x = pump + u
        y = w
        z = zz
        for a in range(3):
            for b in range(3):
                v[a, b] = j[a, 0] * q[0, b] + j[a, 1] * q[1, b] + j[a, 2] * q[2, b]
        for col in range(3):
            for _ in range(2):
                for prev in range(col):
                    dot = (
                        v[0, col] * q[0, prev]
                        + v[1, col] * q[1, prev]
                        + v[2, col] * q[2, prev]
                    )
                    for a in range(3):
                        v[a, col] -= dot * q[a, prev]
            length = math.sqrt(v[0, col] ** 2 + v[1, col] ** 2 + v[2, col] ** 2)
            for a in range(3):
                q[a, col] = v[a, col] / length
            if i >= burn:
                rates[(i - burn) // window, col] += math.log(length) / window
        if i >= count - keep:
            cloud[i - count + keep] = [x, y, z]
    return cloud, rates


@njit(cache=True)
def independent(params, state, burn, steps, keep):
    pump, loss, base, fast, gain, rate = params
    a = complex(state[0], state[1])
    theta = state[2]
    z = np.empty((keep, 3))
    for i in range(burn + steps):
        I = abs(a) ** 2
        s = 1.0 - 1.0 / (1.0 + I)
        a, theta = (
            pump + loss * a * np.exp(1j * (base + fast * s + theta)),
            theta + rate * (gain * s - theta),
        )
        if i >= burn + steps - keep:
            z[i - burn - steps + keep] = [a.real, a.imag, theta]
    return z
