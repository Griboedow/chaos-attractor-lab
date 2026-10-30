"""Validated half-period shooting for the published motor-row polynomial ODE.

Every asserted enclosure uses Arb balls, a Picard existence tube and an
interval Taylor remainder. Ordinary high-precision Newton iterates only
propose centers; the separate Krawczyk inclusion is the existence test.
No NumPy, SciPy, GPU or previous computed trajectory is trusted by the proof.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import platform
import time
from pathlib import Path

import flint
from flint import arb, arb_mat, ctx


ROOT = Path(__file__).resolve().parents[1]
ZERO = arb(0)
ONE = arb(1)


def ball_record(x):
    """Lossless dyadic midpoint/radius, plus a human-readable enclosure."""
    m, e = x.mid().man_exp()
    r, s = x.rad().man_exp()
    return dict(display=x.str(24), midpoint=[str(m), int(e)], radius=[str(r), int(s)])


def matrix_record(a):
    return [[ball_record(a[i, j]) for j in range(a.ncols())]
            for i in range(a.nrows())]


def read_ball(record):
    m,e = record["midpoint"]
    r,s = record["radius"]
    return arb((int(m),int(e)), (int(r),int(s)))


def max_abs(xs):
    return max(x.abs_upper() for x in xs)


class PolynomialODE:
    """Sparse degree-two vector field, optionally with variational blocks.

    A term (c,i,j) denotes c*x_i*x_j; None omits one or both factors.
    The normalized-time duration multiplies the entire field.
    """
    def __init__(self, rho, r, blocks=("core",)):
        self.rho, self.r = rho, r
        self.terms = [
            [(r-rho, 0, None), (-r, 0, 2), (-r, 1, None)],
            [(ONE, 0, None), (-ONE, 0, 2), (-ONE, 1, None)],
            [(ONE, 0, 1), (-ONE, 2, None)],
        ]
        # Each Jacobian entry is an affine polynomial in (h,Y,Z).
        core = [[[(r-rho, None), (-r, 2)], [(-r, None)], [(-r, 0)]],
                [[(ONE, None), (-ONE, 2)], [(-ONE, None)], [(-ONE, 0)]],
                [[(ONE, 1)], [(ONE, 0)], [(-ONE, None)]]]
        self.offsets = {}
        for name in blocks:
            if name == "core":
                a = core
            else:
                # b=1; c_A=1-Z+Y, s_A=1-Z-Y, V_A=-h.
                sign = ONE if name == "A" else -ONE
                c = [(ONE, None), (-ONE, 2), (sign, 1)]
                s = [(ONE, None), (-ONE, 2), (-sign, 1)]
                mul = lambda k, p: [(k*t, i) for t, i in p]
                a = [
                    [[(-ONE, None)], [(sign, 0)]+mul(-r, s), mul(rho, s)],
                    [[(-sign, 0)], [(-ONE, None)]+mul(r, c), mul(-rho, c)],
                    [[], [(r, None)], [(-rho, None)]],
                ]
            off = len(self.terms)
            self.offsets[name] = off
            for i in range(3):
                for j in range(3):
                    self.terms.append([(coef, q, off+3*k+j)
                                       for k in range(3) for coef, q in a[i][k]])

    def initial(self, core):
        return list(core)+[arb(int(i == j)) for _ in self.offsets
                           for i in range(3) for j in range(3)]

    def value(self, x, duration):
        return [duration*sum((c*(ONE if i is None else x[i])
                               *(ONE if j is None else x[j]) for c, i, j in ts), ZERO)
                for ts in self.terms]

    def lipschitz(self, x, duration, weights):
        rows = []
        for ts in self.terms:
            jac = [ZERO for _ in x]
            for c, i, j in ts:
                if i is not None:
                    jac[i] += c*(ONE if j is None else x[j])
                if j is not None:
                    jac[j] += c*(ONE if i is None else x[i])
            row = len(rows)
            rows.append((sum((v.abs_upper()*w for v,w in zip(jac, weights)), ZERO)/weights[row]).upper())
        return duration.abs_upper()*max(rows)

    def jets(self, x, duration, order):
        a = [[v] for v in x]
        for k in range(order):
            for row, ts in enumerate(self.terms):
                total = ZERO
                for c, i, j in ts:
                    if i is None and j is None:
                        v = ONE if k == 0 else ZERO
                    elif i is None or j is None:
                        v = a[j if i is None else i][k]
                    else:
                        v = sum((a[i][l]*a[j][k-l] for l in range(k+1)), ZERO)
                    total += c*v
                a[row].append(duration*total/(k+1))
        return a


def integrate(ode, initial, duration, steps, order, audit=None):
    """Enclose the normalized-time flow on [0,1], including uncertain inputs."""
    dt = ONE/steps
    interval_dt = arb(dt/2, dt/2)
    x = list(initial)
    smallest_margin = None
    largest_lip = ZERO
    largest_radius = ZERO
    largest_remainder = ZERO
    disk_upper = ZERO
    for step in range(steps):
        f = ode.value(x, duration)
        padding = [2*dt*v.abs_upper()+arb("1e-40") for v in f]
        for attempt in range(20):
            tube = [arb(v.mid(), v.rad()+p) for v, p in zip(x, padding)]
            image = [v+interval_dt*w for v, w in zip(x, ode.value(tube, duration))]
            if all(b.contains_interior(v) for b, v in zip(tube, image)):
                break
            padding = [max((arb("5/4")*p).upper(), (2*(v-c).abs_upper()).upper())
                       for p, v, c in zip(padding, image, x)]
        else:
            raise ArithmeticError(f"No Picard tube at step {step}")
        matrix_weight = max(ONE, max_abs(tube[3:])) if len(tube) > 3 else ONE
        weights = [ONE]*3+[matrix_weight]*(len(tube)-3)
        lip = dt*ode.lipschitz(tube, duration, weights)
        if not lip < 1:
            raise ArithmeticError(f"Picard contraction not certified: {lip}, step {step}")
        margins = [min((v.lower()-b.lower()).lower(), (b.upper()-v.upper()).lower())
                   for b, v in zip(tube, image)]
        margin = min(margins)
        smallest_margin = margin if smallest_margin is None else min(smallest_margin, margin)
        largest_lip = max(largest_lip, lip.abs_upper())
        disk_upper = max(disk_upper, ((1-tube[2])**2+tube[1]**2).upper())
        coefficients = ode.jets(x, duration, order)
        remainder_coefficients = ode.jets(tube, duration, order+1)
        out = []
        for a, rem in zip(coefficients, remainder_coefficients):
            value = a[-1]
            for v in reversed(a[:-1]):
                value = value*dt+v
            remainder = rem[-1]*dt**(order+1)
            out.append(value+remainder)
            largest_remainder = max(largest_remainder, remainder.abs_upper())
        x = out
        largest_radius = max(largest_radius, max(v.rad() for v in x))
    if audit is not None:
        audit.append(dict(steps=steps, order=order, dimension=len(x),
                          picard_min_margin=ball_record(smallest_margin),
                          picard_max_contraction=ball_record(largest_lip),
                          maximum_radius=ball_record(largest_radius),
                          maximum_taylor_remainder=ball_record(largest_remainder),
                          core_disk_upper=ball_record(disk_upper)))
    return x


def shooting(guess, rho, r, steps, order, audit=None):
    ode = PolynomialODE(rho, r)
    x = integrate(ode, ode.initial([ZERO, guess[0], guess[1]]), guess[2], steps, order, audit)
    f = PolynomialODE(rho, r, ()).value(x[:3], ONE)
    # Half-period target R(0,Y,Z)=(0,-Y,Z).
    residual = arb_mat([[x[0]], [x[1]+guess[0]], [x[2]-guess[1]]])
    jac = arb_mat([[x[3+3*i+1]+int(i == 1),
                    x[3+3*i+2]-int(i == 2), f[i]] for i in range(3)])
    return residual, jac


def midpoint_matrix(a):
    return arb_mat([[a[i, j].mid() for j in range(a.ncols())] for i in range(a.nrows())])


def certify_spectrum(box, rho, r, steps, order, audit, near_crossing=False):
    ode = PolynomialODE(rho, r, ("core", "A", "B"))
    x = integrate(ode, ode.initial([ZERO, box[0], box[1]]), box[2], steps, order, audit)
    matrices = {name: arb_mat([[x[off+3*i+j] for j in range(3)] for i in range(3)])
                for name, off in ode.offsets.items()}
    symmetry = arb_mat([[-1,0,0], [0,-1,0], [0,0,1]])
    # Equivariance gives the second half from the certified first half.
    core = symmetry*matrices["core"]*symmetry*matrices["core"]
    full = {"A": matrices["B"]*matrices["A"],
            "B": matrices["A"]*matrices["B"]}
    trace_minus_one = sum((core[i,i] for i in range(3)), ZERO)-1
    determinant = core.det()
    schur = [1-determinant, 1-trace_minus_one+determinant, 1+trace_minus_one+determinant]
    core_pass = all(v > 0 for v in schur)
    fixed_brackets = [("-3.71", "-3.70"), ("-0.012", "-0.010"), ("0.420", "0.423")]
    transverse = {}
    transverse_pass = True
    unstable_counts = []
    for name, m in full.items():
        poly = m.charpoly()
        brackets = fixed_brackets
        if near_crossing:
            options = [("-1.001", "-1"), ("-1", "-0.999")]
            first = next(((lo,hi) for lo,hi in options if poly(arb(lo))*poly(arb(hi)) < 0), None)
            brackets = [first or ("-1.001", "-0.999"), ("-0.045", "-0.044"), ("0.587", "0.590")]
        signs = []
        for lo, hi in brackets:
            a, b = poly(arb(lo)), poly(arb(hi))
            passed = (a < 0 and b > 0) or (a > 0 and b < 0)
            signs.append(dict(bracket=[lo,hi], left=ball_record(a), right=ball_record(b), passed=bool(passed)))
            transverse_pass = transverse_pass and passed
        # Brackets touching +/-1 have strict root-interior conclusions.
        outside = sum(arb(hi) <= -1 or arb(lo) >= 1 for lo,hi in brackets)
        inside = sum(arb(lo) >= -1 and arb(hi) <= 1 for lo,hi in brackets)
        transverse_pass = transverse_pass and (outside+inside == 3)
        unstable_counts.append(outside)
        transverse[name] = dict(monodromy=matrix_record(m), root_brackets=signs,
                                unstable_count=outside, stable_count=inside,
                                det_identity_plus_monodromy=ball_record(-poly(-1)),
                                characteristic_coefficients=[ball_record(poly[i]) for i in range(4)])
    unstable = 3*sum(unstable_counts)
    return dict(core_monodromy=matrix_record(core), core_schur_margins=[ball_record(v) for v in schur],
                core_stable=bool(core_pass), transverse=transverse,
                transverse_root_count_certified=bool(transverse_pass),
                period=ball_record(2*box[2]), passive_multiplier=ball_record((-2*box[2]).exp()),
                floquet_count=dict(unstable=unstable, neutral=1, stable=22-unstable) if core_pass and transverse_pass else None)


def certify_target(args):
    ctx.prec = args.bits
    r = 1+arb("19/20")/arb("32001/500000").sqrt()
    rho = arb(arb(args.rho), arb(args.rho_radius))
    history = []
    audit = []
    radius = arb(args.radius).upper()
    if args.verify_input:
        saved = json.loads(args.verify_input.read_text())
        guess = [read_ball(v) for v in saved["shooting_center"]]
        radius = read_ball(saved["shooting_box"][0]).rad()
        rho = read_ball(saved["parameter_rho"])
    else:
        source = json.loads((ROOT/"results/motor-rows/cluster-cycles/analysis.json").read_text())
        case = min(source["cases"], key=lambda c: abs(c["parameters"][0]-float(rho.mid())))
        if args.near_crossing:
            case = json.loads((ROOT/"results/motor-rows/cluster-cycles/crossing-check.json").read_text())
        guess = [arb(str(case["initial_state"][1])).mid(),
                 arb(str(case["initial_state"][2])).mid(), arb(str(case["period"]))/2]
    for iteration in range(0 if args.verify_input else args.newton):
        f, j = shooting(guess, rho.mid(), r, args.steps, args.order, audit)
        correction = midpoint_matrix(j).inv()*midpoint_matrix(f)
        history.append(dict(iteration=iteration, residual=[ball_record(f[i, 0]) for i in range(3)]))
        guess = [(guess[i]-correction[i, 0]).mid() for i in range(3)]
        print(f"Newton {iteration}: residual {max_abs([f[i,0] for i in range(3)])}", flush=True)
    center_f, center_j = shooting(guess, rho, r, args.steps, args.order, audit)
    box = [arb(v, radius) for v in guess]
    _, box_j = shooting(box, rho, r, args.steps, args.order, audit)
    c = midpoint_matrix(midpoint_matrix(center_j).inv())
    defect = arb_mat([[int(i == j) for j in range(3)] for i in range(3)])-c*box_j
    shift = -c*center_f
    delta = arb_mat([[arb(0, radius)] for _ in range(3)])
    image = shift+defect*delta
    ratios = [image[i, 0].abs_upper()/radius for i in range(3)]
    included = all(v < 1 for v in ratios) and not c.det().contains(0)
    contraction = max(sum((defect[i, j].abs_upper() for j in range(3)), ZERO).upper()
                      for i in range(3))
    record = dict(parameter_r=ball_record(r), parameter_rho=ball_record(rho),
                  source_constants="Exact printed decimals; r=1+(19/20)/sqrt(32001/500000)",
                  shooting_center=[ball_record(v) for v in guess],
                  shooting_box=[ball_record(v) for v in box],
                  newton_history=history, krawczyk_inclusion=bool(included),
                  krawczyk_ratios=[ball_record(v) for v in ratios],
                  krawczyk_contraction=ball_record(contraction),
                  preconditioner=matrix_record(c), interval_jacobian=matrix_record(box_j),
                  center_residual=matrix_record(center_f), integrations=audit)
    if not (included and contraction < 1 and box[0] < 0 and box[2] > 0):
        record["status"] = "not_certified"
        return record
    record["status"] = "periodic_orbit_certified"
    print("Periodic orbit Krawczyk inclusion passed", flush=True)
    if args.existence_only:
        return record
    record["spectrum"] = certify_spectrum(box, rho, r, args.steps, args.order, audit, args.near_crossing)
    if record["spectrum"]["floquet_count"] is not None:
        record["status"] = "periodic_orbit_and_floquet_count_certified"
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho", default="1291/440")
    parser.add_argument("--rho-radius", default="0")
    parser.add_argument("--steps", type=int, default=160)
    parser.add_argument("--order", type=int, default=26)
    parser.add_argument("--bits", type=int, default=256)
    parser.add_argument("--newton", type=int, default=4)
    parser.add_argument("--radius", default="1e-20")
    parser.add_argument("--near-crossing", action="store_true")
    parser.add_argument("--existence-only", action="store_true")
    parser.add_argument("--verify-input", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT/"results/motor-proof/target.json")
    args = parser.parse_args()
    if not (args.steps >= 1 and args.order >= 1 and args.bits >= 64 and args.newton >= 0
            and arb(args.radius) > 0 and arb(args.rho_radius) >= 0):
        parser.error("Require steps/order>=1, bits>=64, newton>=0, radius>0 and rho-radius>=0")
    if args.output.exists():
        raise FileExistsError("Preserve previous proof attempts: choose a new output path")
    start = time.monotonic()
    record = dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  arguments={k:str(v) for k,v in vars(args).items()},
                  python=platform.python_version(), python_flint=flint.__version__,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    try:
        record.update(certify_target(args))
    except Exception as error:
        record.update(status="error", error=repr(error))
        raise
    finally:
        record["wall_seconds"] = time.monotonic()-start
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
        args.output.with_suffix(".source.py").write_bytes(Path(__file__).read_bytes())
        print(record.get("status"), args.output, flush=True)
    if not record.get("status", "").endswith("certified") or record.get("status") == "not_certified":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
