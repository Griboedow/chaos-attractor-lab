"""Interval proof for a four-pair rotating periodic motor cycle.

The four equally weighted pair states form an exact 11D restriction of the
23D eight-row system. Quarter-period permutation shooting proves a periodic
orbit, and Schur recursion certifies all nonflow Floquet roots inside the
unit disk. Approximate source data is used only to propose a shooting center.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import platform
import time
from pathlib import Path
from types import SimpleNamespace

import flint
from flint import arb, arb_mat, ctx
import validate_motor_cycle as base
from validate_motor_cycle import (ROOT, ONE, ZERO, PolynomialODE, integrate,
                                 ball_record, matrix_record, midpoint_matrix, read_ball)

DIM=11
PERMUTATION=(2,3,1,0)
FREE=(0,1,2,3,4,5,6,7,9,10)


def identity(n):
    return arb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def simplify(terms):
    combined={}
    for c,i,j in terms:
        if i is None:
            i,j=j,None
        if i is not None and j is not None and j < i:
            i,j=j,i
        combined[i,j]=combined.get((i,j),ZERO)+c
    return [(c,i,j) for (i,j),c in combined.items() if not c.is_zero()]


def permutation_matrix():
    p=arb_mat(DIM,DIM)
    for j,k in enumerate(PERMUTATION):
        p[j,k]=1
        p[4+j,4+k]=1
        if j < 3:
            if k < 3:
                p[8+j,8+k]=1
            else:
                for l in range(3):
                    p[8+j,8+l]=-1
    if p**4 != identity(DIM):
        raise AssertionError("The exact row permutation must have order four")
    return p


class PairedODE(PolynomialODE):
    def __init__(self,rho,r,with_core=True,with_transverse=False,order=28,generic_jets=False):
        self.rho,self.r=rho,r
        self.order=order
        self.generic_jets=generic_jets
        velocities=[]
        for j in range(4):
            v=[(r,4+j)]+[(-r/4,4+k) for k in range(4)]
            v += [(-rho,8+j)] if j < 3 else [(rho,8+k) for k in range(3)]
            velocities.append(v)
        terms=[]
        for j in range(4):
            terms.append(simplify([(-ONE,j,None),(ONE,None,None)]
                                 +[(-c,k,4+j) for c,k in velocities[j]]))
        for j in range(4):
            terms.append(simplify([(-ONE,4+j,None),(ONE,None,None)]
                                 +[(c,k,j) for c,k in velocities[j]]))
        for j in range(3):
            terms.append(simplify([(c,k,None) for c,k in velocities[j]]))
        self.field_terms=terms
        self.terms=list(terms)
        self.offsets={}
        self.matrix_jac={}
        if with_core:
            jac=[[[] for _ in range(DIM)] for _ in range(DIM)]
            for row,ts in enumerate(terms):
                for c,i,j in ts:
                    if i is not None:
                        jac[row][i].append((c,j))
                    if j is not None:
                        jac[row][j].append((c,i))
            self.add_matrix("core",jac,DIM)
        if with_transverse:
            for j in range(4):
                v=velocities[j]
                a=[[[(-ONE,None)],[(-c,k) for c,k in v]+[(-r,4+j)],[(rho,4+j)]],
                   [v,[(-ONE,None),(r,j)],[(-rho,j)]],
                   [[],[(r,None)],[(-rho,None)]]]
                self.add_matrix(str(j),a,3)
        self.row_disk_upper=ZERO

    def add_matrix(self,name,jac,size):
        off=len(self.terms)
        self.offsets[name]=(off,size)
        self.matrix_jac[name]=jac
        for i in range(size):
            for j in range(size):
                self.terms.append(simplify([(c,q,off+size*k+j)
                    for k in range(size) for c,q in jac[i][k]]))

    def initial(self,state):
        return list(state)+[arb(int(i == j)) for _,size in self.offsets.values()
                            for i in range(size) for j in range(size)]

    def lipschitz(self,x,duration,unused_weights):
        matrix_weight=max(ONE,max((v.abs_upper() for v in x[DIM:]),default=ONE))
        weights=[ONE]*DIM+[matrix_weight]*(len(x)-DIM)
        rows=[]
        for row,ts in enumerate(self.terms):
            jac={}
            for c,i,j in ts:
                if i is not None:
                    jac[i]=jac.get(i,ZERO)+c*(ONE if j is None else x[j])
                if j is not None:
                    jac[j]=jac.get(j,ZERO)+c*(ONE if i is None else x[i])
            rows.append((sum((v.abs_upper()*weights[i] for i,v in jac.items()),ZERO)/weights[row]).upper())
        return duration.abs_upper()*max(rows)

    def jets(self,x,duration,order):
        if order == self.order+1:
            self.row_disk_upper=max(self.row_disk_upper,
                max((x[j]**2+x[4+j]**2).upper() for j in range(4)))
        if self.generic_jets:
            return super().jets(x,duration,order)
        # The core is independent of all fundamental matrices. Compute its
        # jets once, then evaluate the identical matrix recurrence in Arb's
        # C matrix arithmetic, avoiding repeated Python scalar convolutions.
        core=PolynomialODE.jets(SimpleNamespace(terms=self.field_terms),x[:DIM],duration,order)
        coefficients=list(core)
        for name,(off,size) in self.offsets.items():
            jac=self.matrix_jac[name]
            jac_jets=[arb_mat([[sum((c*(ONE if k==0 else ZERO)
                if q is None else c*core[q][k] for c,q in jac[i][j]),ZERO)
                for j in range(size)] for i in range(size)]) for k in range(order)]
            matrices=[arb_mat([[x[off+size*i+j] for j in range(size)] for i in range(size)])]
            for k in range(order):
                product=arb_mat(size,size)
                for l in range(k+1):
                    product+=jac_jets[l]*matrices[k-l]
                matrices.append((duration/(k+1))*product)
            coefficients += [[m[i,j] for m in matrices] for i in range(size) for j in range(size)]
        return coefficients


def flow(ode,state,duration,steps,order,audit):
    x=integrate(ode,ode.initial(state),duration,steps,order,audit)
    if audit is not None:
        # The base integrator's three-state disk monitor has no meaning here.
        audit[-1].pop("core_disk_upper",None)
        audit[-1]["four_pair_row_disk_tube_upper"]=ball_record(ode.row_disk_upper)
        audit[-1]["jet_implementation"]="generic sparse" if ode.generic_jets else "Arb matrix convolution"
    return x


def extract(x,entry):
    off,size=entry
    return arb_mat([[x[off+size*i+j] for j in range(size)] for i in range(size)])


def initial_state(q):
    state=[ZERO]*DIM
    for i,j in enumerate(FREE):
        state[j]=q[i]
    return state


def shooting(q,rho,r,steps,order,audit=None,generic=False):
    ode=PairedODE(rho,r,order=order,generic_jets=generic)
    state=initial_state(q)
    x=flow(ode,state,q[-1],steps,order,audit)
    p=permutation_matrix()
    target=p*arb_mat([[v] for v in state])
    residual=arb_mat([[x[i]-target[i,0]] for i in range(DIM)])
    m=extract(x,ode.offsets["core"])
    f=PairedODE(rho,r,with_core=False).value(x[:DIM],ONE)
    jac=arb_mat([[m[i,j]-p[i,j] for j in FREE]+[f[i]] for i in range(DIM)])
    return residual,jac


def schur_certificate(coefficients):
    """Real Schur recursion; no approximate root computation or tolerance."""
    a=list(coefficients)
    records=[]
    while len(a)>1:
        if a[-1].contains(0):
            return dict(passed=False,steps=records,reason="leading coefficient overlaps zero")
        a=[v/a[-1] for v in a]
        a[-1]=ONE  # exact mathematical normalization, not a rounded estimate
        margin=1-a[0].abs_upper()
        records.append(dict(degree=len(a)-1,constant=ball_record(a[0]),
                            margin=ball_record(margin),coefficients=[ball_record(v) for v in a]))
        if not margin > 0:
            return dict(passed=False,steps=records,reason="strict Schur margin not certified")
        n=len(a)-1
        a=[a[k+1]-a[0]*a[n-1-k] for k in range(n)]
    return dict(passed=True,steps=records)


def spectrum(box,rho,r,steps,order,audit,generic=False):
    ode=PairedODE(rho,r,with_transverse=True,order=order,generic_jets=generic)
    x=flow(ode,initial_state(box),box[-1],steps,order,audit)
    p=permutation_matrix()
    quarter=extract(x,ode.offsets["core"])
    monodromy=(p**3*quarter)**4
    polynomial=monodromy.charpoly()
    # Exact nonconstant autonomous orbit supplies the factor z-1.
    quotient=[ZERO]*DIM
    quotient[-1]=polynomial[DIM]
    for k in range(DIM-2,-1,-1):
        quotient[k]=polynomial[k+1]+quotient[k+1]
    remainder=polynomial[0]+quotient[0]
    quotient_schur=schur_certificate(quotient)
    matrices={j:extract(x,ode.offsets[str(j)]) for j in range(4)}
    transverse={}
    for j in range(4):
        path=[j]
        for _ in range(3):
            path.append(PERMUTATION[path[-1]])
        m=identity(3)
        for k in path:
            m=matrices[k]*m
        poly=m.charpoly()
        transverse[str(j)]=dict(quarter_path=path,monodromy=matrix_record(m),
                                schur=schur_certificate([poly[i] for i in range(4)]))
    passed=(remainder.contains(0) and quotient_schur["passed"]
            and all(v["schur"]["passed"] for v in transverse.values()))
    return dict(passed=bool(passed),quotient_quarter_matrix=matrix_record(quarter),
                quotient_full_monodromy=matrix_record(monodromy),
                quotient_characteristic_coefficients=[ball_record(polynomial[i]) for i in range(DIM+1)],
                exact_flow_division_remainder=ball_record(remainder),
                nonflow_quotient_schur=quotient_schur,transverse=transverse,
                floquet_count=dict(unstable=0,neutral=1,stable=22) if passed else None)


def prove(args,record):
    ctx.prec=args.bits
    rho=arb(args.rho)
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    source=ROOT/"results/motor-postcritical/paired-cycle-rho-3p3652184-retry/analysis.json"
    seed=json.loads(source.read_text())
    q=[arb(str(seed["initial_state"][j])).mid() for j in FREE]
    q.append((arb(str(seed["period"]))/4).mid())
    if args.verify_input:
        prior=json.loads(args.verify_input.read_text())
        q=[read_ball(v).mid() for v in prior["shooting_center"]]
        record["center_input"]=str(args.verify_input)
        record["center_input_sha256"]=hashlib.sha256(args.verify_input.read_bytes()).hexdigest()
    audit=[]
    record.update(parameter_rho=ball_record(rho),parameter_r=ball_record(r),
                  exact_rho_decimal=args.rho,seed_file=str(source.relative_to(ROOT)),
                  seed_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  quarter_permutation=list(PERMUTATION),integrations=audit,newton_history=[])
    for k in range(args.newton):
        f,j=shooting(q,rho,r,args.steps,args.order,audit,args.generic_jets)
        correction=midpoint_matrix(j).inv()*midpoint_matrix(f)
        record["newton_history"].append(dict(iteration=k,residual=matrix_record(f)))
        q=[(v-correction[i,0]).mid() for i,v in enumerate(q)]
        print("Newton",k,"max residual",max(f[i,0].abs_upper() for i in range(DIM)),flush=True)
    center_f,center_j=shooting(q,rho,r,args.steps,args.order,audit,args.generic_jets)
    radius=arb(args.radius).upper()
    box=[arb(v,radius) for v in q]
    _,box_j=shooting(box,rho,r,args.steps,args.order,audit,args.generic_jets)
    c=midpoint_matrix(midpoint_matrix(center_j).inv())
    defect=identity(DIM)-c*box_j
    image=-c*center_f+defect*arb_mat([[arb(0,radius)] for _ in range(DIM)])
    ratios=[image[i,0].abs_upper()/radius for i in range(DIM)]
    contraction=max(sum((defect[i,j].abs_upper() for j in range(DIM)),ZERO).upper()
                    for i in range(DIM))
    state=initial_state(box)
    section_speed=PairedODE(rho,r,with_core=False).value(state,ONE)[8]
    passed=(all(v < 1 for v in ratios) and contraction < 1
            and not c.det().contains(0) and box[-1] > 0 and section_speed > 0)
    record.update(shooting_center=[ball_record(v) for v in q],
                  shooting_box=[ball_record(v) for v in box],
                  initial_state_box=[ball_record(v) for v in state],
                  center_residual=matrix_record(center_f),interval_jacobian=matrix_record(box_j),
                  preconditioner=matrix_record(c),krawczyk_ratios=[ball_record(v) for v in ratios],
                  krawczyk_contraction=ball_record(contraction),section_speed=ball_record(section_speed),
                  krawczyk_passed=bool(passed),positive_period=ball_record(4*box[-1]))
    if not passed:
        record["status"]="not_certified"
        return
    record["status"]="periodic_four_pair_orbit_certified"
    print("Quarter-period shooting certified",flush=True)
    if args.existence_only:
        return
    record["spectrum"]=spectrum(box,rho,r,args.steps,args.order,audit,args.generic_jets)
    if record["spectrum"]["passed"]:
        record["status"]="orbitally_stable_four_pair_cycle_certified"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho",default="3.3652183995790423")
    parser.add_argument("--bits",type=int,default=320)
    parser.add_argument("--steps",type=int,default=192)
    parser.add_argument("--order",type=int,default=28)
    parser.add_argument("--radius",default="1e-36")
    parser.add_argument("--newton",type=int,default=3)
    parser.add_argument("--existence-only",action="store_true")
    parser.add_argument("--generic-jets",action="store_true",help="Use the slower scalar sparse Taylor recurrence")
    parser.add_argument("--verify-input",type=Path,help="Use only a previous center as a new initial guess")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if (args.bits < 64 or args.steps < 1 or args.order < 2 or args.newton < 0
            or not arb(args.radius) > 0):
        parser.error("Require bits>=64, steps>=1, order>=2, newton>=0, positive radius")
    if args.output.exists():
        raise FileExistsError("Preserve previous evidence: choose a new output path")
    sources={Path(__file__).name:Path(__file__).read_bytes(),Path(base.__file__).name:Path(base.__file__).read_bytes()}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                arguments={k:str(v) for k,v in vars(args).items()},python=platform.python_version(),
                python_flint=flint.__version__,source_sha256={n:hashlib.sha256(b).hexdigest() for n,b in sources.items()},
                limitations=["No least-period proof for the four-pair orbit","No nonlinear branch birth theorem",
                             "No chaotic attractor or scientific priority claim"])
    start=time.monotonic()
    try:
        prove(args,record)
    except Exception as error:
        record.update(status="error",error=repr(error))
        raise
    finally:
        record["wall_seconds"]=time.monotonic()-start
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
        args.output.with_suffix(".source.py").write_bytes(sources[Path(__file__).name])
        args.output.with_suffix(".base.py").write_bytes(sources[Path(base.__file__).name])
        print(record.get("status"),args.output,flush=True)
    if record.get("status") not in ("periodic_four_pair_orbit_certified","orbitally_stable_four_pair_cycle_certified"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
