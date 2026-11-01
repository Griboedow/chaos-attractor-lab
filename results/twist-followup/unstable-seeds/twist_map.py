"""A contractive translation followed by three coordinate-dependent rotations.

Parameters: r, a_xy, a_yz, a_zx, k_xy, k_yz, k_zx, h_xy, h_yz, h_zx.
Angles are a + k*w + h*w**3, where w is the unchanged rotation axis.
There is no clipping, wrapping of coordinates, noise or continuous-time solver.
"""
import numpy as np
from numba import njit, prange


@njit
def step_tangent(z, q, p):
    """Advance the state and any matrix of tangent columns."""
    r = p[0]
    z = r * z
    z[0] += 1 - r
    q = r * q
    for stage in range(3):
        i, j, k = stage, (stage + 1) % 3, (stage + 2) % 3
        w = z[k]
        angle = p[1 + stage] + p[4 + stage] * w + p[7 + stage] * w * w * w
        derivative = p[4 + stage] + 3 * p[7 + stage] * w * w
        c, s = np.cos(angle), np.sin(angle)
        u, v = c * z[i] - s * z[j], s * z[i] + c * z[j]
        for column in range(q.shape[1]):
            du = c * q[i, column] - s * q[j, column] - v * derivative * q[k, column]
            dv = s * q[i, column] + c * q[j, column] + u * derivative * q[k, column]
            q[i, column], q[j, column] = du, dv
        z[i], z[j] = u, v
    return z, q


@njit
def step(z, p):
    r = p[0]
    z = r * z
    z[0] += 1 - r
    for stage in range(3):
        i, j, k = stage, (stage + 1) % 3, (stage + 2) % 3
        w = z[k]
        angle = p[1 + stage] + p[4 + stage] * w + p[7 + stage] * w * w * w
        c, s = np.cos(angle), np.sin(angle)
        z[i], z[j] = c * z[i] - s * z[j], s * z[i] + c * z[j]
    return z


@njit
def orthonormalize(q):
    """Modified Gram-Schmidt, repeated to control near-dependent columns."""
    logs = np.zeros(3)
    for j in range(3):
        for repeat in range(2):
            for k in range(j):
                product = 0.0
                for i in range(3):
                    product += q[i, j] * q[i, k]
                for i in range(3):
                    q[i, j] -= product * q[i, k]
        norm = np.sqrt((q[:, j] ** 2).sum())
        q[:, j] /= norm
        logs[j] = np.log(norm)
    return q, logs


@njit(parallel=True)
def spectra(parameters, starts, warmup, observed, window, keep):
    """Full spectra, disjoint windows and unthinned final trajectories."""
    count = len(parameters)
    windows = np.zeros((count, observed // window, 3))
    tails = np.zeros((count, keep, 3))
    ends = np.zeros((count, 3))
    means = np.zeros((count, observed // window, 3))
    radii = np.zeros(count)
    for lane in prange(count):
        z = starts[lane].copy()
        q = np.eye(3)
        accumulator = np.zeros(3)
        average = np.zeros(3)
        for t in range(warmup + observed):
            z, q = step_tangent(z, q, parameters[lane])
            q, logs = orthonormalize(q)
            if t >= warmup:
                accumulator += logs
                average += z
                radii[lane] = max(radii[lane], np.sqrt((z * z).sum()))
                n = t - warmup + 1
                if n % window == 0:
                    windows[lane, n // window - 1] = accumulator / window
                    means[lane, n // window - 1] = average / window
                    accumulator[:] = 0
                    average[:] = 0
                if t >= warmup + observed - keep:
                    tails[lane, t - (warmup + observed - keep)] = z
        ends[lane] = z
    return windows, tails, ends, means, radii


@njit
def trajectory(p, start, warmup, observed):
    z = start.copy()
    for _ in range(warmup):
        z = step(z, p)
    result = np.empty((observed, 3))
    for i in range(observed):
        z = step(z, p)
        result[i] = z
    return result
