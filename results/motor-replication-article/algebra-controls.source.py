"""Independent weighted/full-space algebra and exact rational jet controls."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
from flint import arb,arb_mat,ctx,fmpq,fmpq_mat
import validate_motor_five_cycle as proof
from validate_motor_cycle import ZERO,ONE,ball_record
from verify_motor_certificate import exact_bounds

GROUPS=((0,3),(1,2),(4,),(5,6),(7,))


def full_manual(state,rho,r):
    c,s=state[:8],state[8:16];x=list(state[16:])+[-sum(state[16:],ZERO)]
    v=[-rho*x[i]+r*(s[i]-sum(s,ZERO)/8) for i in range(8)]
    field=arb_mat([[-c[i]-v[i]*s[i]+1] for i in range(8)]+
                  [[-s[i]+v[i]*c[i]+1] for i in range(8)]+[[v[i]] for i in range(7)])
    jac=arb_mat(23,23)
    for i in range(8):
        jac[i,i]=-1
        for j in range(8):
            projection=arb(int(i==j))-arb(1)/8
            jac[i,8+j]=-int(i==j)*v[i]-r*s[i]*projection
            jac[8+i,j]=int(i==j)*v[i]
            jac[8+i,8+j]=-int(i==j)+r*c[i]*projection
            if i<7:jac[16+i,8+j]=r*projection
        for j in range(7):
            embedding=arb(int(i==j)) if i<7 else -ONE
            jac[i,16+j]=rho*s[i]*embedding
            jac[8+i,16+j]=-rho*c[i]*embedding
            if i<7:jac[16+i,16+j]=-rho*int(i==j)
    return field,jac


def maps():
    e=arb_mat(23,14);restrict=arb_mat(14,23);contrasts={}
    for j,group in enumerate(GROUPS):
        for row in group:
            e[row,j]=1;e[8+row,5+j]=1
            if row<7:
                if j<4:e[16+row,10+j]=1
                else:
                    for k in range(4):e[16+row,10+k]=-proof.MULTIPLICITIES[k]
        a=group[0];restrict[j,a]=1;restrict[5+j,8+a]=1
        if j<4:restrict[10+j,16+a]=1
        if len(group)==2:
            a,b=group;contrast=arb_mat(23,3);left=arb_mat(3,23)
            for col,offset in enumerate((0,8,16)):
                contrast[offset+a,col]=1;contrast[offset+b,col]=-1
                left[col,offset+a]=1
            contrasts[str(j)]=(contrast,left)
    return e,restrict,contrasts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose an unused output')
    ctx.prec=384;e,restrict,contrasts=maps();p=proof.permutation_matrix()
    errors=[];jet_errors=[]
    for trial in range(1,9):
        state=[arb((i+3)*trial-17)/37 for i in range(14)]
        rho=arb(3)+arb(trial)/11;r=arb(4)+arb(trial)/13
        full=e*arb_mat([[v] for v in state]);f,jac=full_manual([full[i,0] for i in range(23)],rho,r)
        ode=proof.FiveODE(rho,r,with_transverse=True)
        initial=ode.initial(state);values=ode.value(initial,ONE)
        core=proof.extract(values,ode.offsets['core']);expected=restrict*jac*e
        differences=[(arb_mat([[values[i]] for i in range(14)])-restrict*f),core-expected,jac*e-e*core]
        for name,(embedding,left) in contrasts.items():
            matrix=proof.extract(values,ode.offsets[name])
            differences += [matrix-left*jac*embedding,jac*embedding-embedding*matrix]
        moved=p*arb_mat([[v] for v in state])
        field=proof.FiveODE(rho,r,with_core=False)
        differences.append(arb_mat([[v] for v in field.value([moved[i,0] for i in range(14)],ONE)])-
                           p*arb_mat([[v] for v in field.value(state,ONE)]))
        for matrix in differences:
            errors += [(matrix[i,j]).abs_upper() for i in range(matrix.nrows()) for j in range(matrix.ncols())]
        if trial<3:
            fast=ode.jets(initial,arb('1/2'),5);ode.generic_jets=True
            generic=ode.jets(initial,arb('1/2'),5)
            jet_errors += [(a-b).abs_upper() for x,y in zip(fast,generic) for a,b in zip(x,y)]
    equilibrium=[ONE]*10+[ZERO]*4
    full=e*arb_mat([[v] for v in equilibrium]);_,jf=full_manual([full[i,0] for i in range(23)],arb(3),arb(4))
    references={'core':restrict*jf*e}
    references.update({name:left*jf*embedding for name,(embedding,left) in contrasts.items()})
    ode=proof.FiveODE(arb(3),arb(4),with_transverse=True)
    jets=ode.jets(ode.initial(equilibrium),arb('1/2'),10);exact=[]
    for name,(off,size) in ode.offsets.items():
        matrix=fmpq_mat([[references[name][i,j].fmpq() for j in range(size)] for i in range(size)])
        coefficient=fmpq_mat([[int(i==j) for j in range(size)] for i in range(size)])
        for k in range(11):
            for i in range(size):
                for j in range(size):
                    lo,hi=exact_bounds(ball_record(jets[off+size*i+j][k]));exact.append(lo<=coefficient[i,j]<=hi)
            coefficient=matrix*coefficient*fmpq(1,2*(k+1))
    checks=[dict(name='weighted field, quotient and all three internal contrasts in full23D, and permutation equivariance',
                 components=len(errors),maximum_error=ball_record(max(errors)),passed=bool(max(errors)<arb('1e-95'))),
            dict(name='optimized matrix jets versus full generic sparse recurrence',components=len(jet_errors),
                 maximum_error=ball_record(max(jet_errors)),passed=bool(max(jet_errors)<arb('1e-90'))),
            dict(name='exact rational constant-matrix coefficients',components=len(exact),passed=bool(all(exact))),
            dict(name='exact restriction right inverse and permutation order six',passed=bool(restrict*e==proof.identity(14) and p**6==proof.identity(14)))]
    paths=[Path(__file__),Path(proof.__file__),Path(proof.base.__file__),Path(proof.paired.__file__),Path(proof.spectral.__file__)]
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed' if all(x['passed'] for x in checks) else 'failed',
                checks=checks,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    args.output.with_suffix('.source.py').write_bytes(Path(__file__).read_bytes())
    print(record['status'],[(x['name'],x['passed']) for x in checks])
    if record['status']!='passed':raise SystemExit(2)


if __name__=='__main__':main()
