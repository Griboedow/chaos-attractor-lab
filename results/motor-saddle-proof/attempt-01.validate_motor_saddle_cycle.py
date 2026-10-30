"""Interval half-period proof for the unequal six-cluster motor saddle.

The weights (1,1,2,1,1,2)/8 give a 17D quotient and two 3D pair
contrasts. A half-period row permutation closes the complete orbit.
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
import validate_motor_paired_cycle as paired
from validate_motor_cycle import (ROOT, ONE, ZERO, PolynomialODE, integrate,
    ball_record, matrix_record, midpoint_matrix, read_ball)
from validate_motor_paired_cycle import simplify, identity, extract, schur_certificate

K=6
DIM=17
MULTIPLICITIES=(1,1,2,1,1,2)
PERMUTATION=(1,0,2,4,3,5)
FREE=tuple(j for j in range(DIM) if j != 2*K)

def weights():
    return [arb(m)/8 for m in MULTIPLICITIES]

def positions(state):
    w=weights()
    return list(state[2*K:DIM])+[-sum((w[j]*state[2*K+j] for j in range(K-1)),ZERO)/w[-1]]

def permutation_matrix():
    p=arb_mat(DIM,DIM)
    w=weights()
    for j,k in enumerate(PERMUTATION):
        p[j,k]=1
        p[K+j,K+k]=1
        if j < K-1:
            if k < K-1:
                p[2*K+j,2*K+k]=1
            else:
                for l in range(K-1):
                    p[2*K+j,2*K+l]=-w[l]/w[-1]
    assert p**2 == identity(DIM)
    assert all(MULTIPLICITIES[j] == MULTIPLICITIES[k] for j,k in enumerate(PERMUTATION))
    return p

class SaddleODE(PolynomialODE):
    def __init__(self,rho,r,with_core=True,with_transverse=False,order=48,generic_jets=False):
        self.rho,self.r,self.order,self.generic_jets=rho,r,order,generic_jets
        w=weights()
        velocities=[]
        for j in range(K):
            v=[(r,K+j)]+[(-r*w[k],K+k) for k in range(K)]
            v += [(-rho,2*K+j)] if j < K-1 else [(rho*w[k]/w[-1],2*K+k) for k in range(K-1)]
            velocities.append(v)
        terms=[]
        for j in range(K):
            terms.append(simplify([(-ONE,j,None),(ONE,None,None)]+[(-c,k,K+j) for c,k in velocities[j]]))
        for j in range(K):
            terms.append(simplify([(-ONE,K+j,None),(ONE,None,None)]+[(c,k,j) for c,k in velocities[j]]))
        for j in range(K-1):
            terms.append(simplify([(c,k,None) for c,k in velocities[j]]))
        self.field_terms=terms
        self.terms=list(terms)
        self.offsets={}
        self.matrix_jac={}
        if with_core:
            jac=[[[] for _ in range(DIM)] for _ in range(DIM)]
            for row,ts in enumerate(terms):
                for c,i,j in ts:
                    if i is not None: jac[row][i].append((c,j))
                    if j is not None: jac[row][j].append((c,i))
            self.add_matrix("core",jac,DIM)
        if with_transverse:
            for j,m in enumerate(MULTIPLICITIES):
                if m == 1: continue
                v=velocities[j]
                jac=[[[(-ONE,None)],[(-c,k) for c,k in v]+[(-r,K+j)],[(rho,K+j)]],
                     [v,[(-ONE,None),(r,j)],[(-rho,j)]],
                     [[],[(r,None)],[(-rho,None)]]]
                self.add_matrix(str(j),jac,3)
        self.row_disk_upper=ZERO

    add_matrix=paired.PairedODE.add_matrix
    initial=paired.PairedODE.initial

    def lipschitz(self,x,duration,unused_weights):
        matrix_weight=max(ONE,max((v.abs_upper() for v in x[DIM:]),default=ONE))
        scales=[ONE]*DIM+[matrix_weight]*(len(x)-DIM)
        rows=[]
        for row,ts in enumerate(self.terms):
            jac={}
            for c,i,j in ts:
                if i is not None: jac[i]=jac.get(i,ZERO)+c*(ONE if j is None else x[j])
                if j is not None: jac[j]=jac.get(j,ZERO)+c*(ONE if i is None else x[i])
            rows.append((sum((v.abs_upper()*scales[i] for i,v in jac.items()),ZERO)/scales[row]).upper())
        return duration.abs_upper()*max(rows)

    def jets(self,x,duration,order):
        if order == self.order+1:
            self.row_disk_upper=max(self.row_disk_upper,max((x[j]**2+x[K+j]**2).upper() for j in range(K)))
        if self.generic_jets:
            return super().jets(x,duration,order)
        core=PolynomialODE.jets(SimpleNamespace(terms=self.field_terms),x[:DIM],duration,order)
        coefficients=list(core)
        for name,(off,size) in self.offsets.items():
            jac=self.matrix_jac[name]
            jac_jets=[arb_mat([[sum((c*(ONE if k==0 else ZERO) if q is None else c*core[q][k]
                for c,q in jac[i][j]),ZERO) for j in range(size)] for i in range(size)]) for k in range(order)]
            matrices=[arb_mat([[x[off+size*i+j] for j in range(size)] for i in range(size)])]
            for k in range(order):
                product=arb_mat(size,size)
                for l in range(k+1): product+=jac_jets[l]*matrices[k-l]
                matrices.append((duration/(k+1))*product)
            coefficients += [[m[i,j] for m in matrices] for i in range(size) for j in range(size)]
        return coefficients

def flow(ode,state,duration,steps,order,audit):
    x=integrate(ode,ode.initial(state),duration,steps,order,audit)
    if audit is not None:
        audit[-1].pop("core_disk_upper",None)
        audit[-1]["six_cluster_row_disk_tube_upper"]=ball_record(ode.row_disk_upper)
        audit[-1]["jet_implementation"]="generic sparse" if ode.generic_jets else "Arb matrix convolution"
    return x

def initial_state(q):
    state=[ZERO]*DIM
    for i,j in enumerate(FREE): state[j]=q[i]
    return state

def shooting(q,rho,r,steps,order,audit=None,generic=False):
    ode=SaddleODE(rho,r,order=order,generic_jets=generic)
    state=initial_state(q)
    x=flow(ode,state,q[-1],steps,order,audit)
    p=permutation_matrix()
    target=p*arb_mat([[v] for v in state])
    residual=arb_mat([[x[i]-target[i,0]] for i in range(DIM)])
    m=extract(x,ode.offsets["core"])
    f=SaddleODE(rho,r,with_core=False).value(x[:DIM],ONE)
    return residual,arb_mat([[m[i,j]-p[i,j] for j in FREE]+[f[i]] for i in range(DIM)])

def saddle_cubic(poly):
    # One simple real root lies beyond 1. Divide by this exact, enclosed root;
    # interval synthetic division then proves both other roots lie in |z|<.9.
    lo,hi=arb("1.12"),arb("1.13")
    sign_left,sign_right=poly(lo),poly(hi)
    bracket_pass=bool(sign_left < 0 and sign_right > 0)
    root=arb("1.125","0.0050000001")
    quotient=[ZERO]*3
    quotient[-1]=poly[3]
    for k in (1,0): quotient[k]=poly[k+1]+root*quotient[k+1]
    remainder=poly[0]+root*quotient[0]
    stable=schur_certificate(quotient,arb("0.9"))
    return dict(passed=bool(bracket_pass and stable["passed"] and remainder.contains(0)),
        expanding_root_bracket=["1.12","1.13"],left=ball_record(sign_left),right=ball_record(sign_right),
        exact_root_division_remainder=ball_record(remainder),remaining_quadratic_schur=stable)

def spectrum(box,rho,r,steps,order,audit,generic=False):
    ode=SaddleODE(rho,r,with_transverse=True,order=order,generic_jets=generic)
    x=flow(ode,initial_state(box),box[-1],steps,order,audit)
    p=permutation_matrix()
    half=extract(x,ode.offsets["core"])
    monodromy=(p*half)**2
    polynomial=monodromy.charpoly()
    quotient=[ZERO]*DIM
    quotient[-1]=polynomial[DIM]
    for k in range(DIM-2,-1,-1): quotient[k]=polynomial[k+1]+quotient[k+1]
    remainder=polynomial[0]+quotient[0]
    stable=schur_certificate(quotient,arb("0.8"))
    transverse={}
    for j in (2,5):
        assert PERMUTATION[j] == j
        m=extract(x,ode.offsets[str(j)])**2
        poly=m.charpoly()
        test=schur_certificate([poly[k] for k in range(4)],arb("0.9")) if j == 2 else saddle_cubic(poly)
        transverse[str(j)]=dict(monodromy=matrix_record(m),test=test,
            characteristic_coefficients=[ball_record(poly[k]) for k in range(4)])
    passed=bool(remainder.contains(0) and stable["passed"] and all(v["test"]["passed"] for v in transverse.values()))
    return dict(passed=passed,quotient_half_matrix=matrix_record(half),
        quotient_full_monodromy=matrix_record(monodromy),
        quotient_characteristic_coefficients=[ball_record(polynomial[k]) for k in range(DIM+1)],
        exact_flow_division_remainder=ball_record(remainder),nonflow_quotient_schur=stable,
        transverse=transverse,floquet_count=dict(unstable=1,neutral=1,stable=21) if passed else None)

def prove(args,record):
    ctx.prec=args.bits
    rho=arb(args.rho)
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    source=ROOT/"results/motor-postcritical/six-cluster-cycle-rho-3p268/analysis.json"
    seed=json.loads(source.read_text())
    q=[arb(str(seed["initial_state"][j])).mid() for j in FREE]+[(arb(str(seed["period"]))/2).mid()]
    if args.verify_input:
        prior=json.loads(args.verify_input.read_text())
        q=[read_ball(v).mid() for v in prior["shooting_center"]]
        record.update(center_input=str(args.verify_input),center_input_sha256=hashlib.sha256(args.verify_input.read_bytes()).hexdigest())
    audit=[]
    record.update(parameter_rho=ball_record(rho),parameter_r=ball_record(r),exact_rho_decimal=args.rho,
        seed_file=str(source.relative_to(ROOT)),seed_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        half_permutation=list(PERMUTATION),multiplicities=list(MULTIPLICITIES),integrations=audit,newton_history=[])
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
    contraction=max(sum((defect[i,j].abs_upper() for j in range(DIM)),ZERO).upper() for i in range(DIM))
    state=initial_state(box)
    speed=SaddleODE(rho,r,with_core=False).value(state,ONE)[2*K]
    row_x=positions(state)
    distinctions=[]
    for a in range(K):
        for b in range(a+1,K):
            differences=[state[a]-state[b],state[K+a]-state[K+b],row_x[a]-row_x[b]]
            distinctions.append(dict(rows=[a,b],differences=[ball_record(v) for v in differences],
                distinct=bool(any(not v.contains(0) for v in differences))))
    passed=bool(all(v < 1 for v in ratios) and contraction < 1 and not c.det().contains(0)
        and box[-1] > 0 and speed > 0 and all(v["distinct"] for v in distinctions))
    record.update(shooting_center=[ball_record(v) for v in q],shooting_box=[ball_record(v) for v in box],
        initial_state_box=[ball_record(v) for v in state],center_residual=matrix_record(center_f),
        interval_jacobian=matrix_record(box_j),preconditioner=matrix_record(c),
        krawczyk_ratios=[ball_record(v) for v in ratios],krawczyk_contraction=ball_record(contraction),
        section_speed=ball_record(speed),pairwise_cluster_distinctions=distinctions,
        krawczyk_passed=passed,positive_period=ball_record(2*box[-1]))
    record["status"]="periodic_six_cluster_orbit_certified" if passed else "not_certified"
    if not passed or args.existence_only: return
    print("Half-period shooting certified",flush=True)
    record["spectrum"]=spectrum(box,rho,r,args.steps,args.order,audit,args.generic_jets)
    if record["spectrum"]["passed"]: record["status"]="six_cluster_saddle_certified"

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho",default="3.268")
    parser.add_argument("--bits",type=int,default=512)
    parser.add_argument("--steps",type=int,default=768)
    parser.add_argument("--order",type=int,default=48)
    parser.add_argument("--radius",default="1e-50")
    parser.add_argument("--newton",type=int,default=3)
    parser.add_argument("--existence-only",action="store_true")
    parser.add_argument("--generic-jets",action="store_true")
    parser.add_argument("--verify-input",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.bits < 64 or args.steps < 1 or args.order < 2 or args.newton < 0 or not arb(args.radius)>0:
        parser.error("Require bits>=64, steps>=1, order>=2, newton>=0 and positive radius")
    if args.output.exists(): raise FileExistsError("Preserve evidence: use a new output")
    sources={Path(__file__).name:Path(__file__).read_bytes(),Path(base.__file__).name:Path(base.__file__).read_bytes(),
             Path(paired.__file__).name:Path(paired.__file__).read_bytes()}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        arguments={k:str(v) for k,v in vars(args).items()},python=platform.python_version(),python_flint=flint.__version__,
        source_sha256={n:hashlib.sha256(b).hexdigest() for n,b in sources.items()},
        limitations=["No least-period theorem","No heteroclinic connection theorem","No chaotic attractor or priority claim"])
    start=time.monotonic()
    try: prove(args,record)
    except Exception as error:
        record.update(status="error",error=repr(error))
        raise
    finally:
        record["wall_seconds"]=time.monotonic()-start
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
        for name,data in sources.items(): args.output.with_name(args.output.stem+"."+name).write_bytes(data)
        print(record.get("status"),args.output,flush=True)
    if record.get("status") not in ("periodic_six_cluster_orbit_certified","six_cluster_saddle_certified"): raise SystemExit(2)

if __name__ == "__main__": main()
