"""Independent algebraic controls for the interval four-pair certificate."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat
from validate_motor_cycle import ROOT, ONE, ZERO, ball_record
from validate_motor_paired_cycle import PairedODE, permutation_matrix, schur_certificate, simplify


def manual(state,n,rho,r):
    c,s=state[:n],state[n:2*n]
    x=list(state[2*n:])+[-sum(state[2*n:],ZERO)]
    mean=sum(s,ZERO)/n
    v=[-rho*x[j]+r*(s[j]-mean) for j in range(n)]
    f=[-c[j]-v[j]*s[j]+1 for j in range(n)]
    f += [-s[j]+v[j]*c[j]+1 for j in range(n)]
    f += v[:-1]
    jac=arb_mat(3*n-1,3*n-1)
    for j in range(n):
        jac[j,j]=-1
        jac[n+j,j]=v[j]
        for k in range(n):
            projection=arb(int(j == k))-ONE/n
            jac[j,n+k]=-int(j == k)*v[j]-r*s[j]*projection
            jac[n+j,n+k]=-int(j == k)+r*c[j]*projection
        for k in range(n-1):
            embedding=int(j == k) if j < n-1 else -1
            jac[j,2*n+k]=rho*s[j]*embedding
            jac[n+j,2*n+k]=-rho*c[j]*embedding
    for j in range(n-1):
        for k in range(n):
            jac[2*n+j,n+k]=r*(int(j == k)-ONE/n)
        jac[2*n+j,2*n+j]=-rho
    return arb_mat([[v] for v in f]),jac,c,s,x,v


def embedding():
    b=arb_mat(23,23)
    for j in range(8):
        pair=j//2
        b[j,pair]=1
        b[8+j,4+pair]=1
        if j < 7:
            if pair < 3:
                b[16+j,8+pair]=1
            else:
                for k in range(3):
                    b[16+j,8+k]=-1
    for pair in range(4):
        for row,sign in ((2*pair,1),(2*pair+1,-1)):
            b[row,11+3*pair]=sign
            b[8+row,12+3*pair]=sign
            if row < 7:
                b[16+row,13+3*pair]=sign
    return b


def convolution(a,b):
    c=[ZERO]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]+=x*y
    return c


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose unused output to preserve controls")
    ctx.prec=320
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    b=embedding()
    errors=[]
    checks=[]
    uncertain=arb(0,arb("1e-30"))
    retained=simplify([(uncertain,0,None)])
    checks.append(dict(name="coefficient simplification retains balls containing zero",
                       passed=bool(len(retained)==1 and retained[0][0].contains(uncertain))))
    for numerator in range(1,9):
        rho=3+arb(numerator)/17
        state=[arb((j+2)*(numerator+1)-13)/31 for j in range(11)]
        field,jac,c,s,x,v=manual(state,4,rho,r)
        ode=PairedODE(rho,r,with_transverse=True)
        augmented=ode.initial(state)
        values=ode.value(augmented,ONE)
        errors += [(values[j]-field[j,0]).abs_upper() for j in range(11)]
        off,_=ode.offsets["core"]
        errors += [(values[off+11*i+j]-jac[i,j]).abs_upper()
                   for i in range(11) for j in range(11)]
        diag=arb_mat(23,23)
        for i in range(11):
            for j in range(11):
                diag[i,j]=jac[i,j]
        for pair in range(4):
            a=arb_mat([[-1,-v[pair]-r*s[pair],rho*s[pair]],
                       [v[pair],-1+r*c[pair],-rho*c[pair]],[0,r,-rho]])
            off,_=ode.offsets[str(pair)]
            errors += [(values[off+3*i+j]-a[i,j]).abs_upper()
                       for i in range(3) for j in range(3)]
            for i in range(3):
                for j in range(3):
                    diag[11+3*pair+i,11+3*pair+j]=a[i,j]
        paired=b*arb_mat([[v] for v in state]+[[ZERO] for _ in range(12)])
        full,full_jac,*_=manual([paired[j,0] for j in range(23)],8,rho,r)
        embedded_field=b*arb_mat([[field[j,0]] for j in range(11)]+[[ZERO] for _ in range(12)])
        errors += [(full[i,0]-embedded_field[i,0]).abs_upper() for i in range(23)]
        difference=full_jac*b-b*diag
        errors += [difference[i,j].abs_upper() for i in range(23) for j in range(23)]
        p=permutation_matrix()
        permuted=p*arb_mat([[v] for v in state])
        permuted_field=manual([permuted[j,0] for j in range(11)],4,rho,r)[0]
        difference=permuted_field-p*field
        errors += [difference[i,0].abs_upper() for i in range(11)]
    maximum=max(errors)
    checks.append(dict(name="independent quotient, transverse, permutation and full23D block identities",
                       cases=8,components=len(errors),maximum_error=ball_record(maximum),
                       embedding_determinant=ball_record(b.det()),
                       passed=bool(maximum < arb("1e-85") and not b.det().contains(0))))
    # Exactly specified stable roots: 1/2, -3/4, and 3/10 +/- i/5.
    stable=convolution(convolution([arb("-1/2"),ONE],[arb("3/4"),ONE]),
                       [arb("13/100"),arb("-3/5"),ONE])
    stable_proof=schur_certificate(stable)
    unstable_proof=schur_certificate(convolution(stable,[arb("-101/100"),ONE]))
    boundary_proof=schur_certificate(convolution(stable,[-ONE,ONE]))
    checks.append(dict(name="Schur recursion on known real and complex roots, unstable and boundary controls",
                       passed=bool(stable_proof["passed"] and not unstable_proof["passed"]
                                   and not boundary_proof["passed"]),
                       stable=stable_proof,unstable=unstable_proof,boundary=boundary_proof))
    # At c=s=1, X=0, the core is an exact equilibrium and each matrix jet is
    # (L*J)^k/k!. Integer rho=3,r=4 makes independent reference coefficients
    # exact rational matrices, so this is an enclosure test, not a tolerance.
    equilibrium=[ONE]*8+[ZERO]*3
    ode=PairedODE(arb(3),arb(4),with_transverse=True)
    initial=ode.initial(equilibrium)
    degree=12
    fast=ode.jets(initial,arb("1/2"),degree)
    _,jac,c,s,x,v=manual(equilibrium,4,arb(3),arb(4))
    reference_matrices={"core":jac}
    for pair in range(4):
        reference_matrices[str(pair)]=arb_mat([[-1,-4,3],[0,3,-3],[0,4,-3]])
    jet_checks=[]
    for name,(off,size) in ode.offsets.items():
        matrix=reference_matrices[name]
        exact=fmpq_mat([[matrix[i,j].fmpq() for j in range(size)] for i in range(size)])
        coefficient=fmpq_mat([[int(i==j) for j in range(size)] for i in range(size)])
        for k in range(degree+1):
            jet_checks += [fast[off+size*i+j][k].contains(coefficient[i,j])
                           for i in range(size) for j in range(size)]
            coefficient=exact*coefficient*fmpq(1,2*(k+1))
    ode.generic_jets=True
    generic=ode.jets(initial,arb("1/2"),degree)
    jet_error=max((a-b).abs_upper() for row1,row2 in zip(fast,generic) for a,b in zip(row1,row2))
    checks.append(dict(name="matrix Taylor jets enclose exact rational equilibrium coefficients",
                       coefficient_checks=len(jet_checks),degree=degree,
                       generic_recurrence_difference=ball_record(jet_error),
                       passed=bool(all(jet_checks) and jet_error < arb("1e-85"))))
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                status="passed" if all(x["passed"] for x in checks) else "failed",checks=checks,
                source_sha256={p:hashlib.sha256((ROOT/"research"/p).read_bytes()).hexdigest()
                    for p in ("validate_motor_cycle.py","validate_motor_paired_cycle.py","check_motor_paired_proof.py")},
                statement="Independent field/decomposition and analytic polynomial controls, not a replacement for interval shooting and spectral proof.")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    args.output.with_suffix(".source.py").write_bytes(Path(__file__).read_bytes())
    print(record["status"],"components",len(errors),"maximum error",maximum)
    if record["status"] != "passed":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
