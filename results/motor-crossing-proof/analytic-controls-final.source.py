"""Analytic and independent field controls for motor crossing sensitivities.

These controls support implementation review. The crossing proof itself is
the separate interval inclusion and derivative certificate, not these tests.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from flint import arb, arb_mat, ctx
from validate_motor_cycle import ROOT, ONE, ZERO, PolynomialODE, integrate, ball_record
from validate_motor_crossing import SensitivityODE, determinant_direction


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose an unused output path")
    ctx.prec=256
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    errors=[]
    for n in range(1,9):
        rho=2+arb(n)/7
        h,y,z=arb(n)/13,-arb(n)/17,arb(n)/23
        dh,dy,dz=-arb(n)/19,arb(n)/29,arb(n)/31
        alpha=-arb(n)/37
        duration=arb(3)/4
        ode=SensitivityODE(rho,r,("core","A","B"),list(range(3))+list(range(12,30)),alpha)
        x=ode.initial_with_sensitivity([h,y,z],[dh,dy,dz])
        for name,off in ode.offsets.items():
            for i in range(3):
                for j in range(3):
                    x[off+3*i+j]=arb((n+i+2)*(j+1))/41
                    if name != "core":
                        x[ode.sensitivity[off+3*i+j]]=arb((n-j-1)*(i+1))/43
        values=ode.value(x,duration)
        core_f=arb_mat([[(r-rho-r*z)*h-r*y],[h*(1-z)-y],[h*y-z]])
        core_j=arb_mat([[r-rho-r*z,-r,-r*h],[1-z,-1,-h],[y,h,-1]])
        tangent=core_j*arb_mat([[dh],[dy],[dz]])+arb_mat([[-h],[ZERO],[ZERO]])+alpha*core_f
        for i in range(3):
            errors.append((values[i]-duration*core_f[i,0]).abs_upper())
            errors.append((values[ode.sensitivity[i]]-duration*tangent[i,0]).abs_upper())
        for name,off in ode.offsets.items():
            m=arb_mat([[x[off+3*i+j] for j in range(3)] for i in range(3)])
            if name == "core":
                matrix=core_j
            else:
                sign=1 if name == "A" else -1
                c,s,v=1-z+sign*y,1-z-sign*y,-sign*h
                dc,ds,dv=-dz+sign*dy,-dz-sign*dy,-sign*dh
                matrix=arb_mat([[-1,-v-r*s,rho*s],[v,-1+r*c,-rho*c],[0,r,-rho]])
                dmatrix=arb_mat([[0,-dv-r*ds,s+rho*ds],[dv,r*dc,-c-rho*dc],[0,0,-1]])
                dm=arb_mat([[x[ode.sensitivity[off+3*i+j]] for j in range(3)] for i in range(3)])
                expected=duration*(matrix*dm+dmatrix*m+alpha*matrix*m)
                for i in range(3):
                    for j in range(3):
                        errors.append((values[ode.sensitivity[off+3*i+j]]-expected[i,j]).abs_upper())
            expected=duration*matrix*m
            for i in range(3):
                for j in range(3):
                    errors.append((values[off+3*i+j]-expected[i,j]).abs_upper())
    max_error=max(errors)
    checks=[dict(name="independent core, transverse, rho and duration sensitivity fields",
                 cases=8,components=len(errors),maximum_error=ball_record(max_error),
                 passed=bool(max_error < arb("1e-65")))]

    # x'=rho*x, x(0)=1+rho, L=rho/2 at rho=1/3. The total derivative
    # includes both the varying initial state and the varying terminal time.
    # In normalized time, z'=L*(rho*z+2*x). Independent exact solution:
    # x=(4/3)*exp(1/18), z=(13/9)*exp(1/18).
    ode=PolynomialODE(ZERO,ZERO,())
    ode.terms=[[(arb("1/3"),0,None)],[(arb("1/3"),1,None),(arb(2),0,None)],[]]
    result=integrate(ode,[arb("4/3"),ONE,ZERO],arb("1/6"),32,20)
    factor=arb("1/18").exp()
    checks.append(dict(name="linear parameter sensitivity with varying initial state and duration",
                       passed=bool(result[0].contains(arb("4/3")*factor)
                                   and result[1].contains(arb("13/9")*factor)
                                   and result[2].contains(0))))
    a=arb_mat([[0,0,0],[0,2,0],[0,0,3]])
    da=arb_mat([[7,2,3],[5,11,7],[11,13,13]])
    derivative=determinant_direction(a,da)
    checks.append(dict(name="determinant derivative at a singular matrix",
                       value=ball_record(derivative),passed=bool(derivative == 42)))
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                status="passed" if all(c["passed"] for c in checks) else "failed",checks=checks,
                source_sha256={name:hashlib.sha256((ROOT/"research"/name).read_bytes()).hexdigest()
                               for name in ("validate_motor_cycle.py","validate_motor_crossing.py",
                                            "check_motor_crossing.py")},
                statement="Analytic controls and independent field assembly; not a substitute for the interval proof.")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    args.output.with_suffix(".source.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record,indent=2))
    if record["status"] != "passed":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
