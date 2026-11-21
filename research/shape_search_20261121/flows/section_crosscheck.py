"""Count oriented azimuthal half-plane crossings of the five stored 180844 polygons."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def scaled_integers(points):
    """Use one power-of-two denominator for every native binary coordinate."""
    ratios = [[float(x).as_integer_ratio() for x in point] for point in points]
    denominator = max(d for point in ratios for _, d in point)
    vertices = [tuple(n * (denominator // d) for n, d in point) for point in ratios]
    return vertices, denominator


def intersections(edges, ray, denominator, positions=False):
    u, v = ray
    forward = backward = 0
    hits = []
    for (x, y, z), (xx, yy, zz) in edges:
        start = u * y - v * x
        end = u * yy - v * xx
        if not (start <= 0 < end or end <= 0 < start):
            continue
        delta = end - start
        along = (u * x + v * y) * end - (u * xx + v * yy) * start
        if along * delta <= 0:
            continue
        forward += delta > 0
        backward += delta < 0
        if positions:
            # The primary section y=0, x>0 has native coordinates x and z.
            hits.append((Fraction(x * end - xx * start, delta * denominator),
                         Fraction(z * end - zz * start, delta * denominator)))
    return forward, backward, hits


def check(path, rays):
    with np.load(path) as saved:
        points = saved['points']
    vertices, denominator = scaled_integers(points)
    edges = list(zip(vertices, vertices[1:] + vertices[:1]))
    smallest = min(x * yy - y * xx for (x, y, _), (xx, yy, _) in edges)
    if smallest <= 0:
        raise ValueError('Some edge fails strict angular monotonicity')
    forward, backward, hits = intersections(edges, (1, 0), denominator, True)
    for ray in rays:
        if intersections(edges, ray, denominator)[:2] != (forward, backward):
            raise ValueError('The oriented crossing counts change with section angle')
    if backward or len(set(hits)) != forward:
        raise ValueError('Repeated or backward crossing')
    closest_squared = min((a-c)**2 + (b-d)**2 for i, (a, b) in enumerate(hits)
                          for c, d in hits[i+1:])
    return dict(vertices=len(vertices), dtype=str(points.dtype),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                positive_crossings_per_closed_orbit=forward, negative_crossings=backward,
                winding_about_z=forward-backward, distinct_primary_section_points=len(hits),
                minimum_primary_hit_spacing_native=math.sqrt(float(closest_squared)),
                minimum_primary_hit_spacing_squared_exact=str(closest_squared),
                minimum_native_xy_edge_determinant=float(Fraction(smallest, denominator**2)))


def main(args):
    study = json.loads((args.records/'nested-skeletons.json').read_text())
    files = [(1, 'flows/knot-continuation-a/case-180844.npz')]
    files.extend((r['period_order'], r['cloud_relative']) for r in study['cores'])
    rays = [(u, v) for u in range(-3, 4) for v in range(-3, 4)
            if math.gcd(abs(u), abs(v)) == 1]
    rows = []
    for order, relative in files:
        row = check(args.cache_root/relative, rays)
        rows.append(dict(period_order=order, cache_relative=relative, **row))
    base = rows[0]['positive_crossings_per_closed_orbit']
    for row in rows:
        q, remainder = divmod(row['positive_crossings_per_closed_orbit'], base)
        if remainder or q != row['period_order']:
            raise ValueError('The normalized section counts do not reproduce the stored orders')
        row['section_count_relative_to_base'] = q
    report = dict(case=180844, scope='Exact for the closed piecewise-linear polygons through the stored binary coordinates. ODE orbit closure and interpolation accuracy are separate numerical controls.',
                  section='Vertical half-plane y=0, x>0, counting each crossing on one traversal of the closed polygon. Half-open side tests count a section vertex once.',
                  additional_control='All oriented edges have strictly positive exact planar determinant. Hence their azimuth increases along each straight edge, each ray has the same signed count, and there is exactly one positive intersection per azimuthal revolution.',
                  checked_integer_rays=rays, rows=rows)
    args.out.write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cache-root', required=True, type=Path)
    parser.add_argument('--records', type=Path, default=Path('results/shape_search_20261121/flows'))
    parser.add_argument('--out', required=True, type=Path)
    main(parser.parse_args())
