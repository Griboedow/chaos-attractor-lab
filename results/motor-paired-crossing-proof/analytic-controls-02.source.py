"""Independent field and exact-rational jet controls for paired crossing proof."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat
from check_motor_paired_proof import manual
from validate_motor_cycle import ROOT, ONE, ZERO, ball_record
from validate_motor_paired_crossing import SensitivityODE
from verify_motor_certificate import exact_bounds


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError("Choose an unused output")
    ctx.prec=384
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    errors=[]
    jet_errors=[]
    for n in range(1,6):
        rho=3+arb(n)/17
        state=[arb((i+2)*n-9)/31 for i in range(11)]
        direction=[arb((i-5)*(n+2))/43 for i in range(11)]
        alpha=arb(0,arb("0.01")) if n==5 else -arb(n)/29
        duration=arb("3/4")
        f,jac,c,s,x,v=manual(state,4,rho,r)
        wx=direction[8:]+[-sum(direction[8:],ZERO)]
        fp=arb_mat([[x[j]*s[j]] for j in range(4)]
                   +[[-x[j]*c[j]] for j in range(4)]+[[-x[j]] for j in range(3)])
        w=jac*arb_mat([[v] for v in direction])+fp+alpha*f
        ode=SensitivityODE(rho,r,alpha,True)
        initial=ode.initial(state,direction)
        values=ode.value(initial,duration)
        errors += [(values[i]-duration*f[i,0]).abs_upper() for i in range(11)]
        # The uncertain-alpha case requires overlap rather than a tiny-width
        # discrepancy; both direct expressions enclose the same exact family.
        if n<5:
            errors += [(values[11+i]-duration*w[i,0]).abs_upper() for i in range(11)]
        off,size=ode.offsets["core"]
        errors += [(values[off+size*i+j]-duration*jac[i,j]).abs_upper() for i in range(11) for j in range(11)]
        for pair in range(4):
            dv=-x[pair]-rho*wx[pair]+r*(direction[4+pair]-sum(direction[4:8],ZERO)/4)
            a=arb_mat([[-1,-v[pair]-r*s[pair],rho*s[pair]],
                       [v[pair],-1+r*c[pair],-rho*c[pair]],[0,r,-rho]])
            da=arb_mat([[0,-dv-r*direction[4+pair],s[pair]+rho*direction[4+pair]],
                        [dv,r*direction[pair],-c[pair]-rho*direction[pair]],[0,0,-1]])
            expected=arb_mat(6,6)
            for i in range(3):
                for j in range(3):
                    expected[i,j]=a[i,j]
                    expected[3+i,3+j]=a[i,j]
                    expected[3+i,j]=da[i,j]+alpha*a[i,j]
            off,size=ode.offsets[str(pair)]
            if n<5:
                errors += [(values[off+size*i+j]-duration*expected[i,j]).abs_upper() for i in range(6) for j in range(6)]
        if n<3:
            fast=ode.jets(initial,duration,6)
            ode.generic=True
            generic=ode.jets(initial,duration,6)
            jet_errors += [(a-b).abs_upper() for row1,row2 in zip(fast,generic) for a,b in zip(row1,row2)]
    checks=[dict(name="independent state, parameter, duration and triangular matrix sensitivity fields",
                 components=len(errors),maximum_error=ball_record(max(errors)),passed=bool(max(errors)<arb("1e-95"))),
            dict(name="nonlinear matrix sensitivity jets agree with full sparse recurrence",
                 components=len(jet_errors),maximum_error=ball_record(max(jet_errors)),
                 passed=bool(max(jet_errors)<arb("1e-90")))]
    # At the exact equilibrium, the triangular variational matrix is constant.
    # All coefficients then have an independent exact rational oracle.
    equilibrium=[ONE]*8+[ZERO]*3
    ode=SensitivityODE(arb(3),arb(4),arb("1/7"),True)
    fast=ode.jets(ode.initial(equilibrium,[ZERO]*11),arb("1/2"),12)
    j=manual(equilibrium,4,arb(3),arb(4))[1]
    refs={"core":fmpq_mat([[j[i,k].fmpq() for k in range(11)] for i in range(11)])}
    a=fmpq_mat([[-1,-4,3],[0,3,-3],[0,4,-3]])
    da=fmpq_mat([[0,0,1],[0,0,-1],[0,0,-1]])
    triangular=fmpq_mat(6,6)
    for i in range(3):
        for k in range(3):
            triangular[i,k]=a[i,k]
            triangular[3+i,3+k]=a[i,k]
            triangular[3+i,k]=da[i,k]+fmpq(1,7)*a[i,k]
    for row in range(4): refs[str(row)]=triangular
    exact_checks=[]
    for name,(off,size) in ode.offsets.items():
        coefficient=fmpq_mat([[int(i==k) for k in range(size)] for i in range(size)])
        for degree in range(13):
            for i in range(size):
                for k in range(size):
                    lo,hi=exact_bounds(ball_record(fast[off+size*i+k][degree]))
                    exact_checks.append(lo<=coefficient[i,k]<=hi)
            coefficient=refs[name]*coefficient*fmpq(1,2*(degree+1))
    checks.append(dict(name="exact rational coefficients of constant triangular sensitivity blocks",
                       components=len(exact_checks),passed=bool(all(exact_checks))))
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                status="passed" if all(x["passed"] for x in checks) else "failed",checks=checks,
                source_sha256={p:hashlib.sha256((ROOT/"research"/p).read_bytes()).hexdigest()
                    for p in ("validate_motor_paired_crossing.py","check_motor_paired_crossing.py",
                              "validate_motor_paired_cycle.py","check_motor_paired_proof.py",
                              "validate_motor_cycle.py","verify_motor_certificate.py")},
                statement="Independent algebra and analytic jet controls; the interval crossing proof is separate.")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    args.output.with_suffix(".source.py").write_bytes(Path(__file__).read_bytes())
    print(record["status"],[(c["name"],c["passed"]) for c in checks])
    if record["status"]!="passed": raise SystemExit(2)


if __name__=="__main__": main()
