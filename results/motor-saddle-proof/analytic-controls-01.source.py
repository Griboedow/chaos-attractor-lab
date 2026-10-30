"""Independent weighted restriction, full-space and Taylor controls."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
from flint import arb, arb_mat, arb_poly, ctx, fmpq, fmpq_mat
from validate_motor_cycle import ROOT, ONE, ZERO, ball_record
from validate_motor_saddle_cycle import (K,DIM,MULTIPLICITIES,SaddleODE,
    permutation_matrix,saddle_cubic)
from check_motor_paired_proof import manual as full_manual
from verify_motor_certificate import exact_bounds

ROWS=(0,1,2,2,3,4,5,5)

def manual(state,rho,r):
    w=[arb(m)/8 for m in MULTIPLICITIES]
    c,s=state[:K],state[K:2*K]
    x=list(state[2*K:])+[-sum((w[j]*state[2*K+j] for j in range(K-1)),ZERO)/w[-1]]
    mean=sum((w[j]*s[j] for j in range(K)),ZERO)
    v=[-rho*x[j]+r*(s[j]-mean) for j in range(K)]
    field=[1-c[j]-v[j]*s[j] for j in range(K)]+[1-s[j]+v[j]*c[j] for j in range(K)]+v[:-1]
    jac=arb_mat(DIM,DIM)
    for j in range(K):
        jac[j,j]=-1
        jac[K+j,j]=v[j]
        for k in range(K):
            projection=int(j == k)-w[k]
            jac[j,K+k]=-int(j == k)*v[j]-r*s[j]*projection
            jac[K+j,K+k]=-int(j == k)+r*c[j]*projection
        for k in range(K-1):
            embedding=arb(int(j == k)) if j < K-1 else -w[k]/w[-1]
            jac[j,2*K+k]=rho*s[j]*embedding
            jac[K+j,2*K+k]=-rho*c[j]*embedding
    for j in range(K-1):
        for k in range(K): jac[2*K+j,K+k]=r*(int(j == k)-w[k])
        jac[2*K+j,2*K+j]=-rho
    return arb_mat([[z] for z in field]),jac,c,s,x,v

def embedding():
    b=arb_mat(23,23)
    for row,cluster in enumerate(ROWS):
        b[row,cluster]=1
        b[8+row,K+cluster]=1
        if row < 7:
            if cluster < K-1: b[16+row,2*K+cluster]=1
            else:
                for k in range(K-1): b[16+row,2*K+k]=-arb(MULTIPLICITIES[k])/MULTIPLICITIES[-1]
    for index,(a,z) in enumerate(((2,3),(6,7))):
        for row,sign in ((a,1),(z,-1)):
            b[row,DIM+3*index]=sign
            b[8+row,DIM+3*index+1]=sign
            if row < 7: b[16+row,DIM+3*index+2]=sign
    return b

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError("Preserve evidence: choose a new output")
    ctx.prec=384
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    b=embedding()
    checks=[]
    errors=[]
    for numerator in range(1,9):
        rho=3+arb(numerator)/17
        state=[arb((j+2)*(numerator+1)-13)/31 for j in range(DIM)]
        field,jac,c,s,x,v=manual(state,rho,r)
        ode=SaddleODE(rho,r,with_transverse=True)
        values=ode.value(ode.initial(state),ONE)
        errors += [(values[j]-field[j,0]).abs_upper() for j in range(DIM)]
        off,_=ode.offsets["core"]
        errors += [(values[off+DIM*i+j]-jac[i,j]).abs_upper() for i in range(DIM) for j in range(DIM)]
        diag=arb_mat(23,23)
        for i in range(DIM):
            for j in range(DIM): diag[i,j]=jac[i,j]
        for index,cluster in enumerate((2,5)):
            a=arb_mat([[-1,-v[cluster]-r*s[cluster],rho*s[cluster]],
                [v[cluster],-1+r*c[cluster],-rho*c[cluster]],[0,r,-rho]])
            off,_=ode.offsets[str(cluster)]
            errors += [(values[off+3*i+j]-a[i,j]).abs_upper() for i in range(3) for j in range(3)]
            for i in range(3):
                for j in range(3): diag[DIM+3*index+i,DIM+3*index+j]=a[i,j]
        embedded=b*arb_mat([[z] for z in state]+[[ZERO]]*6)
        full,jfull,*_=full_manual([embedded[j,0] for j in range(23)],8,rho,r)
        embedded_field=b*arb_mat([[field[j,0]] for j in range(DIM)]+[[ZERO]]*6)
        errors += [(full[i,0]-embedded_field[i,0]).abs_upper() for i in range(23)]
        difference=jfull*b-b*diag
        errors += [difference[i,j].abs_upper() for i in range(23) for j in range(23)]
        p=permutation_matrix()
        permuted=p*arb_mat([[z] for z in state])
        difference=manual([permuted[j,0] for j in range(DIM)],rho,r)[0]-p*field
        errors += [difference[i,0].abs_upper() for i in range(DIM)]
    maximum=max(errors)
    checks.append(dict(name="weighted field, Jacobian, permutation and full23D decomposition",
        cases=8,components=len(errors),maximum_error=ball_record(maximum),embedding_determinant=ball_record(b.det()),
        passed=bool(maximum < arb("1e-100") and not b.det().contains(0))))
    eq=[ONE]*(2*K)+[ZERO]*(K-1)
    ode=SaddleODE(arb(3),arb(4),with_transverse=True)
    initial=ode.initial(eq)
    degree=12
    fast=ode.jets(initial,arb("1/2"),degree)
    reference={"core":manual(eq,arb(3),arb(4))[1]}
    reference.update({str(j):arb_mat([[-1,-4,3],[0,3,-3],[0,4,-3]]) for j in (2,5)})
    jets=[]
    for name,(off,size) in ode.offsets.items():
        exact=fmpq_mat([[reference[name][i,j].fmpq() for j in range(size)] for i in range(size)])
        coefficient=fmpq_mat([[int(i==j) for j in range(size)] for i in range(size)])
        for k in range(degree+1):
            for i in range(size):
                for j in range(size):
                    lo,hi=exact_bounds(ball_record(fast[off+size*i+j][k]))
                    jets.append(lo <= coefficient[i,j] <= hi)
            coefficient=exact*coefficient*fmpq(1,2*(k+1))
    ode.generic_jets=True
    generic=ode.jets(initial,arb("1/2"),degree)
    err=max((a-z).abs_upper() for u,v in zip(fast,generic) for a,z in zip(u,v))
    checks.append(dict(name="matrix recurrence against exact rational and scalar jets",coefficients=len(jets),
        generic_error=ball_record(err),passed=bool(all(jets) and err < arb("1e-100"))))
    good=arb_poly([arb("-1.125"),1])*arb_poly([arb("-.3"),1])*arb_poly([arb(".2"),1])
    bad=arb_poly([arb("-1.125"),1])*arb_poly([arb("-1.5"),1])*arb_poly([arb(".2"),1])
    good_result,bad_result=saddle_cubic(good),saddle_cubic(bad)
    checks.append(dict(name="one expanding root plus stable quadratic on known polynomials",good=good_result,bad=bad_result,
        passed=bool(good_result["passed"] and not bad_result["passed"])))
    names=("check_motor_saddle_proof.py","validate_motor_saddle_cycle.py","validate_motor_paired_cycle.py",
           "validate_motor_cycle.py","check_motor_paired_proof.py","verify_motor_certificate.py")
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status="passed" if all(c["passed"] for c in checks) else "failed",checks=checks,
        source_sha256={n:hashlib.sha256((ROOT/"research"/n).read_bytes()).hexdigest() for n in names})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    args.output.with_suffix(".source.py").write_bytes(Path(__file__).read_bytes())
    print(record["status"],len(errors),maximum,len(jets),flush=True)
    if record["status"] != "passed": raise SystemExit(2)

if __name__ == "__main__": main()
