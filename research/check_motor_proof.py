"""Independent analytic controls for the interval ODE and certificate plumbing."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from flint import arb, ctx
from validate_motor_cycle import PolynomialODE, integrate, ball_record, read_ball, ROOT
from verify_motor_certificate import exact_bounds, exact_contains


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"results/motor-proof/analytic-controls.json")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose a new control output")
    ctx.prec = 256
    one,zero = arb(1),arb(0)
    records = []
    # Uncertain initial values and duration, checked at all extremal corners.
    ode = PolynomialODE(zero,zero,())
    ode.terms = [[(arb(2),0,None)], [(-one,1,None)], [(one,2,None)]]
    initial = [arb(1,arb("0.001")),arb(2,arb("0.002")),arb(-1,arb("0.001"))]
    duration = arb("0.2",arb("0.001"))
    endpoint = integrate(ode,initial,duration,32,18)
    checks=[]
    for i,rate in enumerate((2,-1,1)):
        for x in (initial[i].lower(),initial[i].upper()):
            for t in (duration.lower(),duration.upper()):
                checks.append(endpoint[i].contains(x*(rate*t).exp()))
    records.append(dict(name="linear flow with interval initial data and duration",passed=all(checks),checks=len(checks)))
    # x'=x², m'=2xm: exact flow and its derivative, with genuine nonzero tails.
    ode.terms = [[(one,0,0)],[(arb(2),0,1)],[]]
    endpoint = integrate(ode,[one,one,zero],arb("1/5"),32,16)
    records.append(dict(name="quadratic flow and variational derivative",passed=endpoint[0].contains(arb("5/4")) and endpoint[1].contains(arb("25/16")) and endpoint[2].contains(0)))
    # An oscillator checks signs and coupling; cos/sin are independent Arb routines.
    ode.terms = [[(-one,1,None)],[(one,0,None)],[]]
    endpoint = integrate(ode,[one,zero,zero],arb("3/5"),32,16)
    records.append(dict(name="harmonic oscillator",passed=endpoint[0].contains(arb("3/5").cos()) and endpoint[1].contains(arb("3/5").sin())))
    # Compare assembled augmented rows with direct mathematical matrices.
    rho,r=arb("1291/440"),1+arb("19/20")/arb("32001/500000").sqrt()
    errors=[]
    for numerator in range(1,9):
        h,y,z=arb(numerator)/7,-arb(numerator)/19,arb(numerator)/23
        model=PolynomialODE(rho,r,("core","A","B"))
        values=model.value(model.initial([h,y,z]),one)
        manual={"core":[[r-rho-r*z,-r,-r*h],[1-z,-1,-h],[y,h,-1]]}
        for sign,name in ((1,"A"),(-1,"B")):
            c,s,v=1-z+sign*y,1-z-sign*y,-sign*h
            manual[name]=[[-1,-v-r*s,rho*s],[v,-1+r*c,-rho*c],[0,r,-rho]]
        for name,m in manual.items():
            off=model.offsets[name]
            errors += [(values[off+3*i+j]-m[i][j]).abs_upper() for i in range(3) for j in range(3)]
    error=max(errors)
    records.append(dict(name="direct core and transverse matrix identities",passed=error<arb("1e-65"),checks=len(errors),max_error=ball_record(error)))
    # Radius deserialization may widen: exact stored-box containment must not.
    example={"midpoint":["0",0],"radius":["1063382397",-123]}
    restored=ball_record(read_ball(example))
    records.append(dict(name="exact dyadic containment rejects constructor widening",
                        passed=exact_contains(restored,example) and not exact_contains(example,restored)))
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                status="passed" if all(v["passed"] for v in records) else "failed",checks=records,
                statement="Analytic controls and independent field assembly; these support software review and do not replace the interval existence proof.",
                code_sha256={p:hashlib.sha256((ROOT/"research"/p).read_bytes()).hexdigest()
                             for p in ("validate_motor_cycle.py","verify_motor_certificate.py","check_motor_proof.py")})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    if result["status"] != "passed":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
