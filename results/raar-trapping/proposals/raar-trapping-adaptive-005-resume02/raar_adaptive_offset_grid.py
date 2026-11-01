"""Curvature-adaptive numerical rasterization; no interval claims."""
import numpy as np
from numba import njit
from propose_raar_trapping_bitset import mark, marked
from search_raar_circles import step
from trace_raar_cycle_seeds import value


@njit(cache=True)
def advance(frontier, mask, shape, origin, h, p, buffer, max_new, tolerance, max_depth, offset):
    ny, nz = shape[1], shape[2]
    new = [np.int64(0) for _ in range(0)]
    minimum = 1e100
    sn, co = np.sin(p[2]), np.cos(p[2])
    splits, leaves, largest_remainder = 0, 0, 0.
    for pos, key in enumerate(frontier):
        ijk = np.array([key//(ny*nz), (key//nz)%ny, key%nz])
        lo, hi = origin+h*(ijk-offset), origin+h*(ijk-offset+1)
        hlo, hhi = co*(lo[1]-p[4])+sn*(lo[2]-p[5]), co*(hi[1]-p[4])+sn*(hi[2]-p[5])
        if lo[0] <= p[3] <= hi[0] and hlo <= 0 <= hhi:
            return np.array(new), 1, pos, key, minimum, splits, leaves, largest_remainder
        mid = lo+.5*h
        stack = [(mid[0], mid[1], mid[2], .5*h, 0)]
        while stack:
            mx, my, mz, half, depth = stack.pop()
            midpoint = np.array([mx, my, mz])
            center, J, axis, flag = step(midpoint, p)
            minimum = min(minimum, axis)
            if flag or not np.isfinite(center).all():
                return np.array(new), 2, pos, key, minimum, splits, leaves, largest_remainder
            remainder = np.zeros(3)
            for k in range(8):
                delta = half*np.array([2*(k&1)-1., 2*((k>>1)&1)-1., 2*((k>>2)&1)-1.])
                y, flag, axis = value(midpoint+delta, p)
                minimum = min(minimum, axis)
                if flag or not np.isfinite(y).all():
                    return np.array(new), 2, pos, key, minimum, splits, leaves, largest_remainder
                remainder = np.maximum(remainder, np.abs(y-center-J@delta))
            remainder = 1.1*remainder+1e-12
            if remainder.max() > tolerance*h and depth < max_depth:
                splits += 1
                for k in range(8):
                    delta = .5*half*np.array([2*(k&1)-1., 2*((k>>1)&1)-1., 2*((k>>2)&1)-1.])
                    v = midpoint+delta
                    stack.append((v[0], v[1], v[2], .5*half, depth+1))
                continue
            leaves += 1
            largest_remainder = max(largest_remainder, remainder.max()/h)
            radius = half*np.sum(np.abs(J), axis=1)+remainder
            a = np.floor((center-radius-origin)/h).astype(np.int64)-buffer+offset
            b = np.floor((center+radius-origin)/h).astype(np.int64)+buffer+offset
            if np.any(a < 0) or np.any(b >= shape):
                return np.array(new), 3, pos, key, minimum, splits, leaves, largest_remainder
            directions = np.zeros((15, 3))
            directions[:3] = np.eye(3)
            k = 3
            for i in range(3):
                for j in range(i+1, 3):
                    directions[k] = np.cross(J[:, i], J[:, j]); k += 1
            for i in range(3):
                for j in range(3):
                    directions[k] = np.cross(J[:, i], np.eye(3)[j]); k += 1
            support = np.zeros(15)
            for k in range(15):
                length = np.linalg.norm(directions[k])
                if length > 1e-15:
                    directions[k] /= length
                support[k] = half*np.sum(np.abs(directions[k]@J))+np.sum(np.abs(directions[k])*(remainder+(.5+buffer)*h))
            for ix in range(a[0], b[0]+1):
                for iy in range(a[1], b[1]+1):
                    for iz in range(a[2], b[2]+1):
                        target = (ix*ny+iy)*nz+iz
                        if marked(mask, target):
                            continue
                        displacement = origin+h*(np.array([ix, iy, iz])-offset+.5)-center
                        intersects = True
                        for k in range(15):
                            if abs(np.dot(directions[k], displacement)) > support[k]:
                                intersects = False; break
                        if intersects:
                            mark(mask, target); new.append(target)
        if len(new) > max_new:
            return np.array(new), 4, pos, key, minimum, splits, leaves, largest_remainder
    return np.array(new), 0, len(frontier), -1, minimum, splits, leaves, largest_remainder
