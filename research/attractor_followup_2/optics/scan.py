"""Certify optical ten-cycle phases at closely spaced slices of the unstable curve."""
import argparse
import json
from pathlib import Path

from flint import arb, ctx

import verify as proof

OUT = Path(__file__).resolve().parents[3] / "results/attractor_followup_2/optics"


def scan(bits, model):
    ctx.prec = bits
    p = proof.cycle.parameters()
    points = proof.fixed_points(p)
    sink, sink_box = points["sink"]
    saddle = points["gap"][1]
    proof.local_sections(p, saddle, model)
    metric, radius = proof.cycle.capture_metric(sink, p)
    basis = proof.axes()
    witnesses = []
    for j in range(41):
        # Signed coordinate along the expanding frame column, in units of 10^-25.
        u = -arb(980+j)/2000
        section = saddle + basis*proof.column([u*arb("1e-25"),
                         proof.centered(arb("1e-44")), proof.centered(arb("1e-44"))])
        n, _ = proof.transport(p, section, sink_box, metric, radius, model)
        witnesses.append({"u_numerator": -980-j, "capture_step": n, "phase": (-n)%10})
    return {"denominator": 2000, "witnesses": witnesses}


def plot(data, path):
    import matplotlib.pyplot as plt
    import numpy as np
    rows = data["witnesses"]
    x = [r["u_numerator"]/data["denominator"] for r in rows]
    y = [r["phase"] for r in rows]
    fig, ax = plt.subplots(figsize=(8,3.5), layout="constrained")
    ax.scatter(x, y, c=y, cmap="tab10", vmin=-.5, vmax=9.5, s=39, zorder=3)
    ax.set(xlabel="Expanding coordinate relative to saddle (units of 10⁻²⁵)",
           ylabel="Ten-cycle phase", yticks=range(10),
           title="A single unstable branch reaches different phases of the same cycle")
    ax.ticklabel_format(useOffset=False, axis="x")
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", alpha=.2)
    fig.savefig(path, dpi=155)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Replace archived results")
    parser.add_argument("--figure", action="store_true", help="Draw the phase samples")
    parser.add_argument("--independent", action="store_true",
                        help="Use the independent complex formula and 512-bit arithmetic")
    args = parser.parse_args()
    receipt = OUT/"phase-slices.json"
    result = scan(512, proof.complex_step) if args.independent else scan(384, proof.cycle.step)
    if args.write:
        receipt.parent.mkdir(parents=True, exist_ok=True)
        receipt.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8",newline="\n")
    else:
        assert result == json.loads(receipt.read_bytes())
    phases = [r["phase"] for r in result["witnesses"]]
    print("Phases:", sorted(set(phases)), "Changes:", sum(x!=y for x,y in zip(phases,phases[1:])))
    if args.figure:
        plot(result, OUT/"phase-slices.png")
