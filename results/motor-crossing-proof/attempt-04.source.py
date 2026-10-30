"""Validate a unique transverse multiplier crossing along the motor cycle.

The old shooting center is only an initial guess. Existence, endpoint signs,
spectral separation, and parameter sensitivities are freshly enclosed with
Arb and the interval Taylor/Picard integrator. Finite differences are unused.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import itertools
import json
import platform
import time
from decimal import Decimal, localcontext
from pathlib import Path
from types import SimpleNamespace

import flint
from flint import arb, arb_mat, ctx

import validate_motor_cycle as base
from validate_motor_cycle import (ONE, ZERO, ROOT, PolynomialODE, integrate,
                                 ball_record, matrix_record, read_ball,
                                 shooting, midpoint_matrix, certify_target)
from verify_motor_certificate import exact_contains


def identity(n):
    return arb_mat([[int(i == j) for j in range(n)] for i in range(n)])


class SensitivityODE(PolynomialODE):
    """Append selected total rho derivatives to a polynomial augmented ODE.

    If x_u = L F(x,rho), then z_u = L(F_x z + F_rho + alpha F),
    where alpha = L'/L. Coefficients of F are affine in rho, so their
    derivatives are obtained exactly at rho=0,1 with r=0. State monomials
    are differentiated by the product rule, including repeated factors.
    """
    def __init__(self, rho, r, blocks, selected, alpha=ZERO):
        super().__init__(rho, r, blocks)
        self.base_dimension = len(self.terms)
        self.sensitivity = {row: self.base_dimension+k
                            for k, row in enumerate(selected)}
        derivative_terms = []
        zero = PolynomialODE(ZERO, ZERO, blocks)
        one = PolynomialODE(ONE, ZERO, blocks)
        for row in selected:
            terms = []
            for (c,i,j), (c0,i0,j0), (c1,i1,j1) in zip(
                    self.terms[row], zero.terms[row], one.terms[row]):
                if (i,j) != (i0,j0) or (i,j) != (i1,j1):
                    raise AssertionError("Inconsistent polynomial term ordering")
                dc = c1-c0
                if dc != 0:
                    terms.append((dc,i,j))
                if alpha != 0:
                    terms.append((alpha*c,i,j))
                if i is not None:
                    terms.append((c,self.sensitivity[i],j))
                if j is not None:
                    terms.append((c,i,self.sensitivity[j]))
            derivative_terms.append(terms)
        self.terms.extend(derivative_terms)

    def initial_with_sensitivity(self, core, core_derivative):
        x = super().initial(core)
        return x+[core_derivative[row] if row < 3 else ZERO
                  for row in self.sensitivity]


def block(x, off):
    return arb_mat([[x[off+3*i+j] for j in range(3)] for i in range(3)])


def derivative_block(x, ode, off):
    return arb_mat([[x[ode.sensitivity[off+3*i+j]] for j in range(3)]
                    for i in range(3)])


def determinant_direction(a, da):
    """Derivative of a 3x3 determinant, valid also at singular matrices."""
    result = ZERO
    for permutation in itertools.permutations(range(3)):
        sign = (-1)**sum(permutation[i] > permutation[j]
                         for i in range(3) for j in range(i+1,3))
        for row in range(3):
            term = arb(sign)
            for i,j in enumerate(permutation):
                term *= da[i,j] if i == row else a[i,j]
            result += term
    return result


def uniform_family(center, rho, r, radius, steps, order, audit):
    """Fresh Krawczyk existence/uniqueness uniformly over the rho ball."""
    center_f_mid, center_j = shooting(center,rho.mid(),r,steps,order,audit)
    # A parameter mean-value enclosure retains the tiny rho interval without
    # repeatedly wrapping it into the state box. Both quantities are validated.
    param_ode = SensitivityODE(rho,r,(),range(3))
    param_x = integrate(param_ode,param_ode.initial_with_sensitivity(
        [ZERO,*center[:2]],[ZERO]*3),center[2],steps,order,audit)
    param_derivative = arb_mat([[param_x[param_ode.sensitivity[i]]] for i in range(3)])
    center_f = center_f_mid+(rho-rho.mid())*param_derivative
    box = [arb(v,radius) for v in center]
    _, box_j = shooting(box,rho,r,steps,order,audit)
    c = midpoint_matrix(midpoint_matrix(center_j).inv())
    defect = identity(3)-c*box_j
    shift = -c*center_f
    displacement = arb_mat([[arb(0,radius)] for _ in range(3)])
    image = shift+defect*displacement
    ratios = [image[i,0].abs_upper()/radius for i in range(3)]
    contraction = max(sum((defect[i,j].abs_upper() for j in range(3)),ZERO).upper()
                      for i in range(3))
    passed = (all(v < 1 for v in ratios) and contraction < 1
              and not c.det().contains(0) and box[0] < 0 and box[2] > 0)
    return box,dict(shooting_center=[ball_record(v) for v in center],
                    shooting_box=[ball_record(v) for v in box],
                    center_residual=matrix_record(center_f),
                    parameter_midpoint_residual=matrix_record(center_f_mid),
                    parameter_center_derivative=matrix_record(param_derivative),
                    interval_jacobian=matrix_record(box_j),
                    preconditioner=matrix_record(c),
                    krawczyk_ratios=[ball_record(v) for v in ratios],
                    krawczyk_contraction=ball_record(contraction),
                    passed=bool(passed))


def branch_sensitivity(box,rho,r,steps,order,audit):
    """Use the implicit shooting equations, not differenced orbits."""
    ode = SensitivityODE(rho,r,("core",),range(3))
    x = integrate(ode,ode.initial_with_sensitivity([ZERO,*box[:2]],[ZERO]*3),
                  box[2],steps,order,audit)
    f = PolynomialODE(rho,r,()).value(x[:3],ONE)
    jac = arb_mat([[x[3+3*i+1]+int(i == 1),
                    x[3+3*i+2]-int(i == 2), f[i]] for i in range(3)])
    partial = arb_mat([[x[ode.sensitivity[i]]] for i in range(3)])
    # Default/lu arb_mat.solve has enclosing arithmetic; never use approx.
    derivative = -jac.solve(partial,algorithm="lu")
    if not all(derivative[i,0].is_finite() for i in range(3)):
        raise ArithmeticError("Nonfinite implicit derivative enclosure")
    return [derivative[i,0] for i in range(3)],dict(
        shooting_jacobian=matrix_record(jac),
        fixed_shooting_parameter_derivative=matrix_record(partial),
        branch_derivative=matrix_record(derivative),
        solve_residual=matrix_record(jac*derivative+partial))


def spectrum_and_derivative(box,derivative,rho,r,steps,order,audit):
    selected = list(range(3))+list(range(12,30))
    ode = SensitivityODE(rho,r,("core","A","B"),selected,derivative[2]/box[2])
    x = integrate(ode,ode.initial_with_sensitivity(
        [ZERO,*box[:2]],[ZERO,*derivative[:2]]),box[2],steps,order,audit)
    matrices = {name:block(x,off) for name,off in ode.offsets.items()}
    dmatrices = {name:derivative_block(x,ode,ode.offsets[name]) for name in ("A","B")}
    symmetry = arb_mat([[-1,0,0],[0,-1,0],[0,0,1]])
    core = symmetry*matrices["core"]*symmetry*matrices["core"]
    core_s = sum((core[i,i] for i in range(3)),ZERO)-1
    core_p = core.det()
    margins = [1-core_p,1-core_s+core_p,1+core_s+core_p]
    records = {}
    for name,left,right in (("A","B","A"),("B","A","B")):
        m = matrices[left]*matrices[right]
        dm = dmatrices[left]*matrices[right]+matrices[left]*dmatrices[right]
        polynomial = m.charpoly()
        root_brackets = [("-1.001","-0.999"),("-0.045","-0.044"),("0.587","0.590")]
        tests=[]
        for lo,hi in root_brackets:
            p,q=polynomial(arb(lo)),polynomial(arb(hi))
            passed=(p < 0 and q > 0) or (q < 0 and p > 0)
            tests.append(dict(interval=[lo,hi],left=ball_record(p),right=ball_record(q),
                              passed=bool(passed)))
        determinant_derivative = determinant_direction(identity(3)+m,dm)
        critical_ball = arb(-1,arb("0.001"))
        p_lambda = polynomial.derivative()(critical_ball)
        p_rho = determinant_direction(critical_ball*identity(3)-m,-dm)
        slope = -p_rho/p_lambda if not p_lambda.contains(0) else None
        at_crossing = polynomial.derivative()(-1)
        slope_at_crossing = (determinant_derivative/at_crossing
                             if not at_crossing.contains(0) else None)
        records[name]=dict(monodromy=matrix_record(m),
                           monodromy_parameter_derivative=matrix_record(dm),
                           characteristic_coefficients=[ball_record(polynomial[i]) for i in range(4)],
                           root_brackets=tests,
                           det_I_plus_M=ball_record((identity(3)+m).det()),
                           det_I_plus_M_derivative=ball_record(determinant_derivative),
                           critical_polynomial_lambda_derivative=ball_record(p_lambda),
                           critical_multiplier_derivative=ball_record(slope) if slope is not None else None,
                           multiplier_derivative_at_crossing=(ball_record(slope_at_crossing)
                                                              if slope_at_crossing is not None else None),
                           simple_roots=bool(all(t["passed"] for t in tests)),
                           strictly_increasing_determinant=bool(determinant_derivative > 0),
                           strictly_increasing_multiplier=bool(slope is not None and slope > 0))
    return dict(core_monodromy=matrix_record(core),
                core_schur_margins=[ball_record(v) for v in margins],
                core_stable=bool(all(v > 0 for v in margins)),
                transverse=records,
                half_period_sensitivity=[ball_record(v) for v in derivative])


def refine_center(center,rho,r,args,record):
    """Newton proposals only; all conclusions require later fresh inclusion."""
    history=[]
    for iteration in range(args.refine_iterations):
        for _ in range(3):
            f,j=shooting(center,rho,r,args.steps,args.order)
            correction=midpoint_matrix(j).inv()*midpoint_matrix(f)
            center=[(v-correction[i,0]).mid() for i,v in enumerate(center)]
        q_prime,_=branch_sensitivity(center,rho,r,args.steps,args.order,None)
        spectrum=spectrum_and_derivative(center,q_prime,rho,r,args.steps,args.order,None)
        transverse=spectrum["transverse"]["A"]
        g=read_ball(transverse["det_I_plus_M"])
        dg=read_ball(transverse["det_I_plus_M_derivative"])
        change=(g/dg).mid()
        history.append(dict(iteration=iteration,rho=ball_record(rho),
                            center=[ball_record(v) for v in center],
                            determinant=ball_record(g),derivative=ball_record(dg)))
        rho=(rho-change).mid()
        center=[(v-change*q_prime[i]).mid() for i,v in enumerate(center)]
        print("Proposed rho",rho.str(45,radius=False),flush=True)
    record["center_proposal_history"]=history
    return center,rho


def run(args,record):
    ctx.prec=args.bits
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    source=ROOT/"results/motor-proof/crossing-family-attempt-01.json"
    saved=json.loads(source.read_text())
    center=[read_ball(v).mid() for v in saved["shooting_center"]]
    rho_center=arb(args.rho_center).mid()
    center,rho_center=refine_center(center,rho_center,r,args,record)
    # Convert the proposed dyadic to a declared finite decimal input. Decimal
    # text defines the mathematical center, not a binary floating-point value.
    rho_text=rho_center.str(40,radius=False) if args.refine_iterations else args.rho_center
    rho=arb(arb(rho_text),arb(args.rho_radius))
    with localcontext() as decimal_context:
        decimal_context.prec=120
        endpoints=[str(Decimal(rho_text)+sign*Decimal(args.rho_radius)) for sign in (-1,1)]
    record.update(parameter_r=ball_record(r),parameter_rho=ball_record(rho),
                  original_guess_file=str(source.relative_to(ROOT)),
                  original_guess_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  exact_center_decimal=rho_text,exact_endpoint_decimals=endpoints,
                  integrations=[])
    audit=record["integrations"]
    box,family=uniform_family(center,rho,r,arb(args.radius).upper(),args.steps,args.order,audit)
    record["family"]=family
    if not family["passed"]:
        record["status"]="family_not_certified"
        return
    print("Uniform periodic family certified",flush=True)
    derivative,sensitivity=branch_sensitivity(box,rho,r,args.steps,args.order,audit)
    record["sensitivity"]=sensitivity
    print("Implicit branch derivative",*[v.str(12) for v in derivative],flush=True)
    spectrum=spectrum_and_derivative(box,derivative,rho,r,args.steps,args.order,audit)
    record["spectrum"]=spectrum
    gates={"family":family["passed"],"core_stable":spectrum["core_stable"]}
    for name,rec in spectrum["transverse"].items():
        for field in ("simple_roots","strictly_increasing_determinant","strictly_increasing_multiplier"):
            gates[name+"_"+field]=rec[field]
        print(name,"g'",rec["det_I_plus_M_derivative"]["display"],
              "multiplier'",rec["critical_multiplier_derivative"],flush=True)
    record["endpoints"]={}
    for name,value,expected_sign in (("low",endpoints[0],-1),("high",endpoints[1],1)):
        endpoint_args=SimpleNamespace(bits=args.bits,rho=value,rho_radius="0",radius=args.endpoint_radius,
                                      verify_input=None,near_crossing=True,newton=3,
                                      steps=args.steps,order=args.order,existence_only=False)
        endpoint=certify_target(endpoint_args)
        record["endpoints"][name]=endpoint
        gates[name+"_certified"]=endpoint["status"]=="periodic_orbit_and_floquet_count_certified"
        gates[name+"_inside_same_shooting_box"]=all(exact_contains(b,x,strict=True)
           for b,x in zip(family["shooting_box"],endpoint["shooting_box"]))
        gates[name+"_rho_in_family"]=exact_contains(record["parameter_rho"],endpoint["parameter_rho"])
        gates[name+"_r_in_family"]=exact_contains(record["parameter_r"],endpoint["parameter_r"])
        for block_name,block_record in endpoint["spectrum"]["transverse"].items():
            g=read_ball(block_record["det_identity_plus_monodromy"])
            gates[name+"_"+block_name+"_determinant_sign"]=expected_sign*g > 0
    record["gates"]={k:bool(v) for k,v in gates.items()}
    record["status"]="unique_transverse_crossing_certified" if all(gates.values()) else "not_certified"
    record["statement"]=("Within the declared exact rho bracket, the validated local periodic family "
        "has exactly one parameter with transverse multiplier -1. Its algebraic and geometric "
        "multiplicity in the full 23D system are six. Its derivative with respect to rho is strictly "
        "positive. These conclusions hold only if every recorded gate passes.")
    record["not_proved"]=["least period","nonlinear equivariant flip or doubled branches",
                          "global branch connection","chaotic attractor","scientific novelty"]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits",type=int,default=256)
    parser.add_argument("--steps",type=int,default=160)
    parser.add_argument("--order",type=int,default=26)
    parser.add_argument("--radius",default="1e-10")
    parser.add_argument("--rho-center",default="3.370218399579")
    parser.add_argument("--rho-radius",default="1e-11")
    parser.add_argument("--endpoint-radius",default="1e-20")
    parser.add_argument("--refine-iterations",type=int,default=0)
    parser.add_argument("--output",type=Path,default=ROOT/"results/motor-crossing-proof/attempt-01.json")
    args=parser.parse_args()
    if (args.bits < 64 or args.steps < 1 or args.order < 1 or not arb(args.radius) > 0
            or not arb(args.rho_radius) > 0 or not arb(args.endpoint_radius) > 0
            or args.refine_iterations < 0):
        parser.error("Require bits>=64, steps/order>=1, positive radii, nonnegative refinement count")
    if args.output.exists():
        raise FileExistsError("Choose an unused output path to preserve previous evidence")
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                arguments={k:str(v) for k,v in vars(args).items()},
                python=platform.python_version(),python_flint=flint.__version__,
                source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__),Path(base.__file__),
                                         ROOT/"research/verify_motor_certificate.py")})
    start=time.monotonic()
    try:
        run(args,record)
    except Exception as error:
        record.update(status="error",error=repr(error))
        raise
    finally:
        record["wall_seconds"]=time.monotonic()-start
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
        args.output.with_suffix(".source.py").write_bytes(Path(__file__).read_bytes())
        args.output.with_suffix(".base.py").write_bytes(Path(base.__file__).read_bytes())
        print(record.get("status"),args.output,flush=True)
    if record.get("status") != "unique_transverse_crossing_certified":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
