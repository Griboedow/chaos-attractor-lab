"""Verify first entry of seven launch neighborhoods into the Case 13745 attractor."""

import json
import sys
from pathlib import Path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "research/low_hanging_20261122/annuli"))
import verify as prior

PACKET = ROOT / "results/attractor_followup_2/annuli"


def half_width(box):
    return max(box.real.rad().upper(), box.imag.rad().upper())


def intervals_covering_preimages(target, edge, value, jacobian):
    """Cover every solution F(z) in target with |z| in the given annulus."""
    stack = [acb(arb(0, 1), arb(0, 1))]
    leaves = []
    while stack:
        domain = stack.pop()
        if abs(domain) < edge[0] or abs(domain) > edge[1]:
            continue
        if not value(domain).overlaps(target):
            continue
        center = domain.mid()
        # Any root in domain belongs to this Krawczyk enclosure by the
        # mean-value identity. An invertible preconditioner is optional here.
        J = jacobian(center)
        if abs(J.det()) > 0:
            A = J.inv().mid()
            delta = value(center) - target
            width = domain - center
            K = (arb_mat([[center.real], [center.imag]])
                 - A * arb_mat([[delta.real], [delta.imag]])
                 + (prior.IDENTITY - A*jacobian(domain)) * arb_mat([[width.real], [width.imag]]))
            if not K[0, 0].overlaps(domain.real) or not K[1, 0].overlaps(domain.imag):
                continue
            cropped = acb(K[0, 0].intersection(domain.real), K[1, 0].intersection(domain.imag))
            if half_width(cropped) < half_width(domain)*prior.exact("0.7"):
                domain = cropped
                if half_width(domain) > prior.exact("1e-9"):
                    stack.append(domain)
                    continue
        if half_width(domain) < prior.exact("2e-9"):
            leaves.append(domain)
            continue
        x, y = domain.real, domain.imag
        if x.rad() > y.rad():
            stack.extend((acb(x.lower().union(x.mid()), y), acb(x.mid().union(x.upper()), y)))
        else:
            stack.extend((acb(x, y.lower().union(y.mid())), acb(x, y.mid().union(y.upper()))))
    return leaves


def exclude_previous(box, edges, depth, value, jacobian):
    """Any point of the invariant set must admit a prehistory of every length."""
    r = abs(box)
    if all(r < a or r > b for a, b in edges):
        return
    ring, = [i for i, (a, b) in enumerate(edges) if r > a and r < b]
    boxes = [box]
    for _ in range(depth):
        ring = 1-ring
        boxes = [p for t in boxes for p in intervals_covering_preimages(t, edges[ring], value, jacobian)]
        if not boxes:
            return
    raise ValueError("A possible prehistory survived the exclusion test.")


def verify_capture(build, bits):
    ctx.prec = bits
    base = json.loads((prior.PACKET / "witness.json").read_text())
    trap = json.loads(prior.TRAPS.read_text())["13745"]
    assert base["case"] == 13745 and float.fromhex(trap["binary64_hex"][0]) == 17
    value, jacobian = build([prior.exact(float.fromhex(x)) for x in trap["binary64_hex"][1:]])
    outer = acb(*map(prior.exact, base["centers"][0]))
    outer_width = prior.exact(base["half_width_in_billionths"][0])/10**9
    edges = [tuple(map(prior.exact, edge)) for edge in trap["annuli"]]
    paths = json.loads((PACKET / "capture.json").read_text())
    assert len(paths) == len(trap["entrance_tests"]) == 7
    ratios = []
    for i, (w, launch) in enumerate(zip(paths, trap["entrance_tests"])):
        assert w["seed_idx"] == i
        k = w["start_step"]
        start = acb(*map(prior.exact, launch["start"]))
        neighborhood = prior.exact("1e-12")
        target, derivative = prior.iterate(start, neighborhood, k, value, jacobian)
        radius = prior.exact(w["target_radius"])
        # Every start in this closed neighborhood enters the target's OPEN interior.
        assert prior.norm(derivative)*neighborhood < radius
        if k:
            prev, derivative = prior.iterate(start, neighborhood, k-1, value, jacobian)
            previous_box = prior.square(prev, prior.norm(derivative)*neighborhood)
            exclude_previous(previous_box, edges, {1: 3, 5: 5}.get(i, 0), value, jacobian)

        for steps, coords, extent, tiles in w["path"]:
            assert steps > 0 and tiles > 0
            source = acb(*[prior.exact(float.fromhex(t)) for t in coords])
            width = prior.exact(extent)
            assert width > 0

            def image(c, r):
                return prior.iterate(c, r, steps, value, jacobian)

            ratio, _ = prior.inverse_branch(image, source, width, target, radius, tiles)
            ratios.append(ratio)
            target, radius = source, width

        rotation = (prior.IMAG * 2*arb.pi()*w["rotation"]/17).exp()
        offset = target/rotation - outer
        assert max(abs(offset.real).upper(), abs(offset.imag).upper()) + arb(2).sqrt()*radius < outer_width
    print(f"{build.__name__}, {bits} bits: first entries {[w['start_step'] for w in paths]}, "
          f"{len(ratios)} covers, largest inclusion ratio {float(max(ratios)):.6f}.")


if __name__ == "__main__":
    verify_capture(prior.polynomial, 192)
    verify_capture(prior.expanded_polynomial, 256)
