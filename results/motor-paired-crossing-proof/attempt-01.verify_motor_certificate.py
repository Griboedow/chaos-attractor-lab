"""Recompute all interval certificates and link the common periodic branch.

This reruns the same validated algorithm with different order, subdivision
and precision. It is a reproducibility check, not an independent integrator.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx, fmpq
from validate_motor_cycle import ROOT, ball_record, read_ball


def exact_bounds(record):
    """Read stored dyadics without Arb constructor radius inflation."""
    m,e = record["midpoint"]
    r,s = record["radius"]
    midpoint = fmpq(int(m))*fmpq(2)**int(e)
    radius = fmpq(int(r))*fmpq(2)**int(s)
    return midpoint-radius, midpoint+radius


def exact_contains(outer, inner, strict=False):
    lo,hi = exact_bounds(outer)
    a,b = exact_bounds(inner)
    return lo < a and b < hi if strict else lo <= a and b <= hi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"results/motor-proof/verification")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--link-only", action="store_true",
                        help="Recheck exact branch membership of existing recomputations; no new integration")
    args = parser.parse_args()
    if args.output.exists() and not args.link_only:
        raise FileExistsError("Choose an unused output directory to preserve evidence")
    if args.workers < 1:
        parser.error("workers must be positive")
    args.output.mkdir(parents=True, exist_ok=args.link_only)
    ctx.prec = 320
    specs = [
        ("target", "target-recheck.json", []),
        ("family", "crossing-family-attempt-01.json", ["--near-crossing", "--existence-only"]),
        ("low", "crossing-low.json", ["--near-crossing"]),
        ("high", "crossing-high.json", ["--near-crossing"]),
    ]
    def run(spec):
        name, source, extra = spec
        source = ROOT/"results/motor-proof"/source
        output = args.output/(name+".json")
        command = [sys.executable, str(ROOT/"research/validate_motor_cycle.py"),
                   "--verify-input", str(source), "--steps", "160", "--order", "26",
                   "--bits", "320", "--output", str(output)]+extra
        result = subprocess.run(command, capture_output=True, text=True)
        (args.output/(name+".log")).write_text(result.stdout+result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"{name} recomputation failed; read {name}.log")
        record = json.loads(output.read_text())
        print(name, record["status"], flush=True)
        return name, record, dict(input=source.name, input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                                  output=output.name, output_sha256=hashlib.sha256(output.read_bytes()).hexdigest())
    if args.link_only:
        results = []
        for name, source, _ in specs:
            output = args.output/(name+".json")
            results.append((name,json.loads(output.read_text()),
                            dict(output=output.name,output_sha256=hashlib.sha256(output.read_bytes()).hexdigest())))
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(run, specs))
    records = {name:rec for name,rec,_ in results}
    family = records["family"]
    gates = {}
    for name in ("target", "family", "low", "high"):
        rec = records[name]
        gates[name+"_inclusion"] = rec["krawczyk_inclusion"] and all(read_ball(v)<1 for v in rec["krawczyk_ratios"])
        gates[name+"_contraction"] = read_ball(rec["krawczyk_contraction"])<1
    for name in ("low", "high"):
        rec = records[name]
        gates[name+"_same_shooting_branch"] = all(exact_contains(b,x,strict=True)
                                                 for b,x in zip(family["shooting_box"],rec["shooting_box"]))
        gates[name+"_rho_in_family"] = exact_contains(family["parameter_rho"],rec["parameter_rho"])
        gates[name+"_r_in_family"] = exact_contains(family["parameter_r"],rec["parameter_r"])
    expected = {"target":{"unstable":6,"neutral":1,"stable":16},
                "low":{"unstable":6,"neutral":1,"stable":16},
                "high":{"unstable":0,"neutral":1,"stable":22}}
    for name, count in expected.items():
        gates[name+"_floquet_count"] = records[name]["spectrum"]["floquet_count"] == count
    g = {name:read_ball(records[name]["spectrum"]["transverse"]["A"]["det_identity_plus_monodromy"])
         for name in ("low","high")}
    gates["opposite_determinant_signs"] = g["low"] < 0 and g["high"] > 0
    rho_low = arb("3.370218399569")
    rho_high = arb("3.370218399589")
    # The recomputed endpoint records enclose these exact decimal rationals.
    gates["exact_rational_endpoints"] = all(exact_bounds(records[name]["parameter_rho"])[0] <= fmpq(Fraction(value).numerator,Fraction(value).denominator)
                                           <= exact_bounds(records[name]["parameter_rho"])[1]
                                           for name,value in (("low","3.370218399569"),("high","3.370218399589")))
    r = read_ball(family["parameter_r"])
    delta_range = (((r-1)/rho_high)**2-1).union(((r-1)/rho_low)**2-1)
    summary = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   status="crossing_existence_certified" if all(gates.values()) else "not_certified",
                   gates={k:bool(v) for k,v in gates.items()}, files=[item for _,_,item in results],
                   rho_bracket=["3.370218399569", "3.370218399589"], delta_enclosure=ball_record(delta_range),
                   endpoint_det_I_plus_M={name:ball_record(v) for name,v in g.items()},
                   conclusion="A continuous, locally unique half-period-symmetric periodic family has at least one transverse multiplier -1 in the open rho bracket; cluster symmetry repeats it at least six times.",
                   not_proved=["uniqueness or nonzero derivative of the crossing", "least period", "nonlinear flip bifurcation", "hyperchaotic attractor", "scientific priority"],
                   verification="Same algorithm recomputed at 320 bits, 160 steps, order 26; no stored true/false flag is a substitute for recomputation.")
    destination = args.output/("crossing-link-exact.json" if args.link_only else "crossing-link.json")
    if destination.exists():
        raise FileExistsError("Preserve previous branch-link records")
    destination.write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(summary["status"], flush=True)
    if not all(gates.values()):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
