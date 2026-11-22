"""Verify cellulations by gluing sign charts of the sevenfold torus cover."""

import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path

# Side i joins corner i to i+1. Coordinates are ordered x,y,w.
SIDES = ((0, 1), (2, -1), (1, 1), (0, -1), (2, 1), (1, -1))
DECK = (-1, 1, 1)


def neighbor(chart, side, sheets):
    signs, sheet = chart
    coord, cos_sign = SIDES[side % 6]
    new = list(signs)
    new[coord] *= -1
    shift = DECK[coord] * signs[coord] if cos_sign < 0 else 0
    return tuple(new), (sheet + shift) % sheets


class Partition:
    def __init__(self):
        self.parent = {}

    def find(self, value):
        root = self.parent.setdefault(value, value)
        if root != value:
            self.parent[value] = self.find(root)
        return self.parent[value]

    def join(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def compute(sheets, extras):
    charts = list(product(product([-1, 1], repeat=3), range(sheets)))
    nodes, faces = Partition(), Partition()
    face_degrees = {}

    def corner_face(chart, corner, side):
        sector = side if extras and corner in (0, 3) else "*"
        return chart, corner, sector

    edges = []
    for chart in charts:
        for i in range(6):
            adj = neighbor(chart, i, sheets)
            assert neighbor(adj, i, sheets) == chart
            nodes.join((chart, "M", i), (adj, "M", i))
            for corner in (i, (i + 1) % 6):
                if extras and corner in (0, 3):
                    nodes.join((chart, "V", corner), (adj, "V", corner))
                f = corner_face(chart, corner, i)
                faces.join(f, corner_face(adj, corner, i))
                face_degrees[f] = 2 if extras and corner in (0, 3) else 1
            # One oval arc per chart sector. Its corner index picks the
            # straight continuation at each fold tangency.
            lm, rm = (chart, "M", (i - 1) % 6), (chart, "M", i)
            if extras and i in (0, 3):
                p = (chart, "P", i)
                edges.extend([(lm, p, "oval", i), (p, rm, "oval", i)])
            else:
                edges.append((lm, rm, "oval", i))
        if extras:
            edges.extend(
                [
                    ((chart, "V", 0), (chart, "P", 0), "diameter", -1),
                    ((chart, "P", 0), (chart, "P", 3), "diameter", -1),
                    ((chart, "P", 3), (chart, "V", 3), "diameter", -1),
                ]
            )
        for j in range(2 if extras else 1):
            f = (chart, "central", j)
            faces.find(f)
            face_degrees[f] = 5 if extras else 6
    incidence = defaultdict(list)
    for eid, (a, b, _, _) in enumerate(edges):
        for vertex in (a, b):
            incidence[nodes.find(vertex)].append((eid, vertex))
    assert all(len(ports) == 4 for ports in incidence.values())
    components = Partition()
    for a, b, _, _ in edges:
        components.join(nodes.find(a), nodes.find(b))
    assert len({components.find(v) for v in incidence}) == 1
    paths = Partition()
    for ports in incidence.values():
        for eid, (chart, kind, i) in ports:
            family, corner = edges[eid][2:]

            def continues(port):
                other, (chart2, _, _) = port
                f2, c2 = edges[other][2:]
                if other == eid or f2 != family:
                    return False
                if kind == "P":
                    return True
                if kind == "M":
                    return chart2 != chart and c2 != corner
                expected = neighbor(neighbor(chart, i - 1, sheets), i, sheets)
                return chart2 == expected

            matches = [p for p in ports if continues(p)]
            assert len(matches) == 1
            paths.join(eid, matches[0][0])
    orbit_lengths = Counter(
        (fam, paths.find(i)) for i, (_, _, fam, _) in enumerate(edges)
    )
    contacts = Counter()
    for ports in incidence.values():
        on_site = {paths.find(eid) for eid, _ in ports}
        assert len(on_site) == 2
        contacts[tuple(sorted(on_site))] += 1
    pair_table = defaultdict(Counter)
    for (i, j), count in contacts.items():
        families = []
        for route in (i, j):
            kind = edges[route][2]
            if kind == "diameter" and sheets == 7:
                kind = (
                    "short diameter"
                    if orbit_lengths[kind, route] == 6
                    else "long diameter"
                )
            families.append(kind)
        pair_table[" / ".join(sorted(families))][count] += 1
    assert sum(
        counts * len_pair
        for _, table in pair_table.items()
        for len_pair, counts in table.items()
    ) == len(incidence)
    perimeter = Counter()
    for face, number in face_degrees.items():
        perimeter[faces.find(face)] += number
    degrees = Counter(perimeter.values())
    v = len(incidence)
    e = len(edges)
    f = len(perimeter)
    assert v - e + f == -4 * sheets
    assert sum(k * n for k, n in degrees.items()) == 2 * e
    assert sum(orbit_lengths.values()) == e
    return dict(
        sheets=sheets,
        genus=2 * sheets + 1,
        rays=len(orbit_lengths),
        vertices=v,
        edges=e,
        pair_counts_by_intersection_multiplicity={
            kind: dict(sorted(counts.items()))
            for kind, counts in sorted(pair_table.items())
        },
        cells=f,
        faces_by_edges=dict(sorted(degrees.items())),
        nodes_by_kind=dict(sorted(Counter(vertex[1] for vertex in incidence).items())),
        ray_lengths_in_graph_edges={
            kind: dict(
                sorted(
                    Counter(
                        count for (fk, _), count in orbit_lengths.items() if fk == kind
                    ).items()
                )
            )
            for kind in sorted({k[0] for k in orbit_lengths})
        },
        euler_characteristic=v - e + f,
    )


def native_contacts(b, sheets=7, c=3 / 8):
    """Independently lift analytic crossing coordinates and test every ray in the torus."""
    from math import acos, cos, pi, sin, sqrt, tau

    p = (c * c - 1 - b * b) / 2
    h = sqrt(b * b - p * p)
    U = acos(p / b)
    V = -acos((-1 - p) / c)
    amplitude_sites = []
    for nx, ny in ((1, 0), (0, 1), (-1, -1)):
        x = nx + p * ny
        y = p * nx + b * b * ny
        scale = sqrt(nx * x + ny * y)
        amplitude_sites.extend(((x / scale, y / scale), (-x / scale, -y / scale)))
    if b == 1:
        amplitude_sites.extend(((1, -1), (-1, 1), (h / c, -h / c), (-h / c, h / c)))

    def angles(cosine):
        if abs(abs(cosine) - 1) < 1e-12:
            return (0.0,) if cosine > 0 else (pi,)
        a = acos(cosine)
        return (a, tau - a)

    def match(v):
        return abs(sin(v / 2)) < 1e-9

    pairs = Counter()
    seen = set()
    for A, B in amplitude_sites:
        triples = product(angles(A), angles(B / b), angles((-A - B) / c))
        for x, y, w in triples:
            for sheet in range(sheets):
                z = ((w - x + y + tau * sheet) / sheets) % tau
                native = (x, y, z)
                key = tuple(round(v, 9) for v in native)
                assert key not in seen
                seen.add(key)
                assert abs(cos(x) + b * cos(y) + c * cos(x - y + sheets * z)) < 1e-11
                rays = []
                for s, t, e in product((-1, 1), repeat=3):
                    if match(y - s * (x + e * U)) and match(w - t * (x + e * V)):
                        rays.append(("oval", s, t, e))
                if b == 1:
                    for w0, turn in ((pi / 2, 1), (3 * pi / 2, -1)):
                        if not match(w - w0):
                            continue
                        if match(y + x - pi):
                            rays.append(("long diameter", turn))
                        if match(y - x - pi):
                            lap = round((sheets * z - w0 - pi) / tau) % sheets
                            rays.append(("short diameter", turn, lap))
                assert len(rays) == 2
                pairs[tuple(sorted(rays))] += 1
    table = defaultdict(Counter)
    for (a, b), hits in pairs.items():
        table[" / ".join([a[0], b[0]])][hits] += 1
    return len(seen), {k: dict(v) for k, v in table.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    cases = {
        name: compute(n, extra)
        for n in (1, 7)
        for name, extra in [(f"generic-{n}", False), (f"equal-{n}", True)]
    }
    for name, b in [("generic-7", 9 / 8), ("equal-7", 1)]:
        count, pairs = native_contacts(b)
        assert (
            count == cases[name]["vertices"]
            and pairs == cases[name]["pair_counts_by_intersection_multiplicity"]
        )
        cases[name]["native_torus_coordinate_check"] = True
    for key, values in cases.items():
        print(key, values)
    if args.out:
        args.out.write_text(json.dumps(cases, indent=2) + "\n", encoding="utf-8")
