"""Outward rasterization of RAAR image enclosures onto a packed cubical grid.

Cell coordinates are exact origin + binary64_step * integer. Subdivision uses
integer addresses at dyadic depths. Separating directions are arbitrary chosen
binary64 vectors: their construction need not be exact. Only outward dot
products may reject target cells. The positive image margin is included before
rasterization. A completed forward closure still needs its provenance and
processed-cell accounting audited before an invariant-set claim is made.
"""
import numpy as np
from numba import njit
from propose_raar_trapping_bitset import marked, mark
from search_raar_circles import step
from raar_interval64 import add, sub, mul, up, down
from raar_geometric_enclosure import point, remainder


@njit(cache=True)
def coordinate_interval(origin, step_size, address):
    return add(point(origin), mul(point(step_size), point(address)))


@njit(cache=True)
def target_separated(index, origin, h, center, directions, support):
    displacement = np.empty((3, 2))
    for j in range(3):
        displacement[j] = sub(coordinate_interval(origin[j], h, index[j]+.5), point(center[j]))
    for k in range(len(directions)):
        projection = point(0.)
        for j in range(3):
            projection = add(projection, mul(point(directions[k, j]), displacement[j]))
        if projection[0] > support[k] or projection[1] < -support[k]:
            return True
    return False


@njit(cache=True)
def raster_parameters(center, J, radii, errors, origin, h):
    directions = np.zeros((15, 3))
    directions[:3] = np.eye(3)
    k = 3
    for i in range(3):
        for j in range(i+1, 3):
            directions[k] = np.cross(J[:, i], J[:, j]); k += 1
    for i in range(3):
        for j in range(3):
            directions[k] = np.cross(J[:, i], np.eye(3)[j]); k += 1
    support = np.empty(15)
    for k in range(15):
        length = np.linalg.norm(directions[k])
        if length > 1e-15:
            directions[k] /= length
        bound = 0.
        for j in range(3):
            projected_column = point(0.)
            for i in range(3):
                projected_column = add(projected_column, mul(point(directions[k, i]), point(J[i, j])))
            bound = up(bound+up(max(abs(projected_column[0]), abs(projected_column[1]))*radii[j]))
        for i in range(3):
            bound = up(bound+up(abs(directions[k, i])*up(errors[i]+.5*h)))
        support[k] = bound
    lower_index, upper_index = np.empty(3, np.int64), np.empty(3, np.int64)
    for i in range(3):
        radius = errors[i]
        for j in range(3):
            radius = up(radius+up(abs(J[i, j])*radii[j]))
        lower, upper = down(center[i]-radius), up(center[i]+radius)
        lower_index[i] = int(np.floor(down(down(lower-origin[i])/h)))
        upper_index[i] = int(np.floor(up(up(upper-origin[i])/h)))
    return lower_index, upper_index, directions, support


@njit(cache=True)
def advance(frontier, mask, shape, origin, h, p, trig, margin, max_new, tolerance, max_depth):
    ny, nz = shape[1], shape[2]
    new = [np.int64(0) for _ in range(0)]
    minimum, largest_remainder = 1e100, 0.
    splits, leaves = 0, 0
    for pos, key in enumerate(frontier):
        ix, iy, iz = key//(ny*nz), (key//nz)%ny, key%nz
        stack = [(ix, iy, iz, 0)]
        while stack:
            sx, sy, sz, depth = stack.pop()
            indices = np.array([sx, sy, sz])
            width = h/(2**depth)
            box, mid = np.empty((3, 2)), np.empty(3)
            for j in range(3):
                box[j, 0] = coordinate_interval(origin[j], width, float(indices[j]))[0]
                box[j, 1] = coordinate_interval(origin[j], width, float(indices[j]+1))[1]
                mid[j] = origin[j]+width*(indices[j]+.5)
                assert box[j, 0] <= mid[j] <= box[j, 1]
            center, J, _, flag = step(mid, p)
            if flag or not np.isfinite(center).all() or not np.isfinite(J).all():
                return np.array(new), 2, pos, key, minimum, splits, leaves, largest_remainder
            radii, errors, axis, flag = remainder(box, mid, center, J, p, trig, margin)
            # Split unresolved regularity bounds too; a failed extension does
            # not by itself establish a singularity in the underlying box.
            if (flag or errors.max() > tolerance*h) and depth < max_depth:
                splits += 1
                for k in range(8):
                    stack.append((2*sx+(k&1), 2*sy+((k>>1)&1), 2*sz+((k>>2)&1), depth+1))
                continue
            if flag or not np.isfinite(errors).all():
                return np.array(new), 1, pos, key, minimum, splits, leaves, largest_remainder
            minimum = min(minimum, axis)
            largest_remainder = max(largest_remainder, errors.max()/h)
            leaves += 1
            a, b, directions, support = raster_parameters(center, J, radii, errors, origin, h)
            if np.any(a < 0) or np.any(b >= shape):
                return np.array(new), 3, pos, key, minimum, splits, leaves, largest_remainder
            for tx in range(a[0], b[0]+1):
                for ty in range(a[1], b[1]+1):
                    for tz in range(a[2], b[2]+1):
                        target = (tx*ny+ty)*nz+tz
                        if marked(mask, target):
                            continue
                        if not target_separated(np.array([tx, ty, tz]), origin, h, center, directions, support):
                            mark(mask, target); new.append(target)
        if len(new) > max_new:
            return np.array(new), 4, pos, key, minimum, splits, leaves, largest_remainder
    return np.array(new), 0, len(frontier), -1, minimum, splits, leaves, largest_remainder
