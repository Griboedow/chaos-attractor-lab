"""Certify the least period using validated half-orbit sign enclosures.

Revalidates the half-period shooting root, then records Picard tubes from
the existing interval Taylor integrator. This is an additional theorem
using the same integrator, not an independent numerical implementation.
"""
import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace

from flint import arb, ctx
from validate_motor_cycle import (ROOT, PolynomialODE, ball_record,
                                  certify_target, integrate, read_ball)


class TubeRecorder(PolynomialODE):
    def __init__(self, rho, r, order):
        super().__init__(rho, r, ())
        self.order = order
        self.tubes = []

    def jets(self, x, duration, order):
        # integrate calls order+1 exactly once per step, on its Picard tube,
        # after checking tube inclusion and the contraction inequality.
        if order == self.order+1:
            self.tubes.append(list(x))
        return super().jets(x, duration, order)


def prove(args):
    ctx.prec = args.bits
    anchor = certify_target(SimpleNamespace(
        bits=args.bits, rho="1291/440", rho_radius="0", radius="1e-20",
        verify_input=args.input, near_crossing=False, newton=0,
        steps=args.steps, order=args.order, existence_only=True))
    if anchor["status"] != "periodic_orbit_certified":
        raise ArithmeticError("The input shooting root did not revalidate")
    box = [read_ball(x) for x in anchor["shooting_box"]]
    rho, r = read_ball(anchor["parameter_rho"]), read_ball(anchor["parameter_r"])
    ode = TubeRecorder(rho, r, args.order)
    audit = []
    endpoint = integrate(ode, [arb(0), box[0], box[1]], box[2],
                         args.steps, args.order, audit)
    if len(ode.tubes) != args.steps:
        raise ArithmeticError("The expected validated tube sequence was not recorded")
    fields = [ode.value(tube, arb(1))[0] for tube in ode.tubes]
    prefix = 0
    while prefix < args.steps and fields[prefix] > 0:
        prefix += 1
    suffix = args.steps
    while suffix > prefix and fields[suffix-1] < 0:
        suffix -= 1
    interior = [tube[0] for tube in ode.tubes[prefix:suffix]]
    gates = dict(positive_initial_derivative=prefix > 0,
                 negative_final_derivative=suffix < args.steps,
                 nonempty_interior=bool(interior),
                 positive_interior=all(h > 0 for h in interior),
                 endpoint_encloses_zero=endpoint[0].contains(0),
                 positive_half_period=box[2] > 0,
                 negative_initial_Y=box[0] < 0)
    records = []
    for i, (tube, field) in enumerate(zip(ode.tubes, fields)):
        kind = "positive_derivative_prefix" if i < prefix else (
            "negative_derivative_suffix" if i >= suffix else "positive_h_interior")
        records.append(dict(step=i, normalized_interval=[f"{i}/{args.steps}",f"{i+1}/{args.steps}"],
                            proof_region=kind, h=ball_record(tube[0]), h_prime=ball_record(field)))
    return dict(status="least_period_certified" if all(gates.values()) else "not_certified",
                gates={k:bool(v) for k,v in gates.items()},
                prefix_steps=prefix, interior_steps=suffix-prefix, suffix_steps=args.steps-suffix,
                minimum_interior_h=ball_record(min(h.lower() for h in interior)) if interior else None,
                positive_period=ball_record(2*box[2]), anchor_revalidation=anchor,
                integration=audit, tube_signs=records,
                theorem="For every parameter in the revalidated family, h(t)>0 for 0<t<L. "
                        "Half-period symmetry gives h(t)<0 for L<t<2L. The sole upward "
                        "zero in [0,2L) is t=0, so the least positive period is exactly 2L.",
                limitations=["Uses the same interval integrator as the cycle certificate",
                             "Does not prove a nonlinear bifurcation or a chaotic attractor"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=160)
    parser.add_argument("--order", type=int, default=26)
    parser.add_argument("--bits", type=int, default=320)
    args = parser.parse_args()
    if args.steps < 4 or args.order < 2 or args.bits < 64:
        parser.error("Require steps>=4, order>=2 and bits>=64")
    if args.output.exists():
        raise FileExistsError("Choose a new path to preserve prior evidence")
    start = time.monotonic()
    record = dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  arguments={k:str(v) for k,v in vars(args).items()},
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  source_sha256={p:hashlib.sha256((ROOT/"research"/p).read_bytes()).hexdigest()
                                 for p in ("validate_motor_cycle.py", "validate_motor_period.py")})
    try:
        record.update(prove(args))
    except Exception as error:
        record.update(status="error", error=repr(error))
        raise
    finally:
        record["wall_seconds"] = time.monotonic()-start
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
        args.output.with_suffix(".source.py").write_bytes(Path(__file__).read_bytes())
        print(record.get("status"), flush=True)
    if record["status"] != "least_period_certified":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
