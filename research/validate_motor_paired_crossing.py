"""Validate the second, fourfold transverse crossing of the paired motor cycle.

Uniform quarter-period shooting, parameter sensitivities, root separation,
and endpoint determinant signs use interval Taylor/Picard enclosures. Newton
proposals are only guesses. No finite-difference derivative proves a claim.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import platform
import time
from decimal import Decimal, localcontext
from pathlib import Path
from types import SimpleNamespace

import flint
from flint import arb, arb_mat, ctx
import validate_motor_cycle as base
import validate_motor_paired_cycle as paired
from validate_motor_cycle import ROOT, ONE, ZERO, ball_record, matrix_record, read_ball, midpoint_matrix
from validate_motor_crossing import determinant_direction
from verify_motor_certificate import exact_contains

N=11


def affine_jets(jac,core,order):
    size=len(jac)
    return [arb_mat([[sum((c*(ONE if k==0 else ZERO) if q is None else c*core[q][k]
                           for c,q in jac[i][j]),ZERO)
                       for j in range(size)] for i in range(size)]) for k in range(order)]


class SensitivityODE(base.PolynomialODE):
    """Eleven states and their total rho derivatives, plus matrix blocks.

    Each transverse 6x6 triangular block propagates B and its derivative D:
    [A,0; A'+alpha*A,A], alpha=L'/L. The lower-left fundamental block is D.
    """
    def __init__(self,rho,r,alpha=ZERO,matrices=False,order=64,generic=False):
        self.rho,self.r,self.alpha,self.order=rho,r,alpha,order
        self.generic=generic
        model=paired.PairedODE(rho,r,with_transverse=True)
        self.base_terms=model.field_terms
        self.core_jac=model.matrix_jac["core"]
        xs=[[(ONE,8+j)] if j<3 else [(-ONE,8+k) for k in range(3)] for j in range(4)]
        forcing=[]
        for j in range(4): forcing.append([(c,k,4+j) for c,k in xs[j]])
        for j in range(4): forcing.append([(-c,k,j) for c,k in xs[j]])
        for j in range(3): forcing.append([(-c,k,None) for c,k in xs[j]])
        self.parameter_terms=forcing
        tangent=[]
        for i in range(N):
            ts=[(c,q,N+k) for k in range(N) for c,q in self.core_jac[i][k]]
            ts+=forcing[i]
            ts += [(alpha*c,j,k) for c,j,k in self.base_terms[i]]
            tangent.append(paired.simplify(ts))
        self.field_terms=self.base_terms+tangent
        self.terms=list(self.field_terms)
        self.offsets={}
        self.matrix_jac={}
        if matrices:
            paired.PairedODE.add_matrix(self,"core",self.core_jac,N)
            for j in range(4):
                a=model.matrix_jac[str(j)]
                dv=[(-c,k) for c,k in xs[j]]+[( -rho*c,N+k) for c,k in xs[j]]
                dv += [(r,N+4+j)]+[(-r/4,N+4+k) for k in range(4)]
                da=[[[],[(-c,k) for c,k in dv]+[(-r,N+4+j)],[(ONE,4+j),(rho,N+4+j)]],
                    [dv,[(r,N+j)],[(-ONE,j),(-rho,N+j)]],
                    [[],[],[(-ONE,None)]]]
                cmat=[[[] for _ in range(6)] for _ in range(6)]
                for row in range(3):
                    for col in range(3):
                        cmat[row][col]=list(a[row][col])
                        cmat[3+row][3+col]=list(a[row][col])
                        cmat[3+row][col]=da[row][col]+[(alpha*c,k) for c,k in a[row][col]]
                paired.PairedODE.add_matrix(self,str(j),cmat,6)
        self.row_disk_upper=ZERO

    def initial(self,state,direction):
        return list(state)+list(direction)+[arb(int(i==j)) for _,size in self.offsets.values()
                                           for i in range(size) for j in range(size)]

    def lipschitz(self,x,duration,weights):
        # Any positive diagonal scaling gives a valid weighted sup norm.
        # Separate B and D scales avoid charging a large derivative forcing
        # to the same scale as its small fundamental-matrix input.
        tangent_weight=max(ONE,max(v.abs_upper() for v in x[N:2*N]))
        weights=[ONE]*N+[tangent_weight]*N
        for name,(off,size) in self.offsets.items():
            if name=="core":
                scale=max(ONE,max(v.abs_upper() for v in x[off:off+size*size]))
                weights += [scale]*(size*size)
            else:
                diagonal=max(ONE,max(x[off+size*i+j].abs_upper()
                    for i in range(size) for j in range(size) if (i<3)==(j<3)))
                derivative=max(ONE,max(x[off+size*i+j].abs_upper()
                    for i in range(3,6) for j in range(3)))
                zero_block=max(arb("1e-30"),max(x[off+size*i+j].abs_upper()
                    for i in range(3) for j in range(3,6)))
                weights += [diagonal if (i<3)==(j<3) else derivative if i>=3 else zero_block
                            for i in range(size) for j in range(size)]
        rows=[]
        for row,ts in enumerate(self.terms):
            jac={}
            for c,i,j in ts:
                if i is not None: jac[i]=jac.get(i,ZERO)+c*(ONE if j is None else x[j])
                if j is not None: jac[j]=jac.get(j,ZERO)+c*(ONE if i is None else x[i])
            rows.append((sum((v.abs_upper()*weights[i] for i,v in jac.items()),ZERO)/weights[row]).upper())
        return duration.abs_upper()*max(rows)

    def jets(self,x,duration,order):
        if order==self.order+1:
            self.row_disk_upper=max(self.row_disk_upper,max((x[j]**2+x[4+j]**2).upper() for j in range(4)))
        if self.generic:
            return super().jets(x,duration,order)
        core=base.PolynomialODE.jets(SimpleNamespace(terms=self.base_terms),x[:N],duration,order)
        jac=affine_jets(self.core_jac,core,order)
        tangent=[arb_mat([[x[N+i]] for i in range(N)])]
        for k in range(order):
            total=arb_mat(N,1)
            for l in range(k+1): total+=jac[l]*tangent[k-l]
            for row,ts in enumerate(self.parameter_terms):
                for c,i,j in ts:
                    total[row,0]+=c*(core[i][k] if j is None else
                        sum((core[i][l]*core[j][k-l] for l in range(k+1)),ZERO))
            next_coefficient=(duration/(k+1))*total
            for i in range(N): next_coefficient[i,0]+=self.alpha*core[i][k+1]
            tangent.append(next_coefficient)
        all_core=core+[[v[i,0] for v in tangent] for i in range(N)]
        coefficients=list(all_core)
        for name,(off,size) in self.offsets.items():
            field_jets=jac if name=="core" else affine_jets(self.matrix_jac[name],all_core,order)
            matrices=[arb_mat([[x[off+size*i+j] for j in range(size)] for i in range(size)])]
            for k in range(order):
                product=arb_mat(size,size)
                for l in range(k+1): product+=field_jets[l]*matrices[k-l]
                matrices.append((duration/(k+1))*product)
            coefficients += [[v[i,j] for v in matrices] for i in range(size) for j in range(size)]
        return coefficients


def integrate_sensitivity(q,rho,r,args,direction=None,matrices=False,audit=None):
    derivative=[ZERO]*N if direction is None else direction
    alpha=derivative[-1]/q[-1]
    ode=SensitivityODE(rho,r,alpha,matrices,args.order)
    initial=ode.initial(paired.initial_state(q),paired.initial_state(derivative))
    x=base.integrate(ode,initial,q[-1],args.steps,args.order,audit)
    if audit is not None:
        audit[-1].pop("core_disk_upper",None)
        audit[-1]["four_pair_row_disk_tube_upper"]=ball_record(ode.row_disk_upper)
        audit[-1]["jet_implementation"]="Arb matrix and affine sensitivity convolutions"
    return x,ode


def implicit_derivative(q,rho,r,jac,args,audit):
    x,_=integrate_sensitivity(q,rho,r,args,audit=audit)
    partial=arb_mat([[x[N+i]] for i in range(N)])
    derivative=-jac.solve(partial,algorithm="lu")
    if not all(derivative[i,0].is_finite() for i in range(N)):
        raise ArithmeticError("Implicit derivative is not finite")
    return [derivative[i,0] for i in range(N)],dict(
        fixed_shooting_rho_derivative=matrix_record(partial),branch_derivative=matrix_record(derivative),
        shooting_jacobian=matrix_record(jac),solve_residual=matrix_record(jac*derivative+partial))


def spectral_analysis(core,quarters,derivatives=None):
    p=paired.permutation_matrix()
    monodromy=(p**3*core)**4
    polynomial=monodromy.charpoly()
    quotient=[ZERO]*N
    quotient[-1]=polynomial[N]
    for k in range(N-2,-1,-1): quotient[k]=polynomial[k+1]+quotient[k+1]
    remainder=polynomial[0]+quotient[0]
    core_schur=paired.schur_certificate(quotient,arb("0.95"))
    records={}
    for row in range(4):
        path=[row]
        for _ in range(3): path.append(paired.PERMUTATION[path[-1]])
        m=paired.identity(3)
        dm=arb_mat(3,3)
        for j in path:
            if derivatives is not None: dm=derivatives[j]*m+quarters[j]*dm
            m=quarters[j]*m
        poly=m.charpoly()
        tests=[]
        for lo,hi in (("-1.001","-0.999"),("-0.0021","-0.0018"),("0.320","0.323")):
            a,b=poly(arb(lo)),poly(arb(hi))
            tests.append(dict(interval=[lo,hi],left=ball_record(a),right=ball_record(b),
                              passed=bool((a<0 and b>0) or (a>0 and b<0))))
        g=-poly(-1)
        rec=dict(quarter_path=path,monodromy=matrix_record(m),
                 characteristic_coefficients=[ball_record(poly[k]) for k in range(4)],
                 root_brackets=tests,simple_roots=bool(all(t["passed"] for t in tests)),
                 determinant=ball_record(g),determinant_sign=1 if g>0 else -1 if g<0 else 0)
        if derivatives is not None:
            dg=determinant_direction(paired.identity(3)+m,dm)
            plambda=poly.derivative()(-1)
            slope=dg/plambda if not plambda.contains(0) else None
            rec.update(monodromy_derivative=matrix_record(dm),determinant_derivative=ball_record(dg),
                       determinant_strictly_increasing=bool(dg>0),
                       characteristic_lambda_derivative_at_minus_one=ball_record(plambda),
                       multiplier_derivative_at_crossing=ball_record(slope) if slope is not None else None,
                       multiplier_crossing_speed_positive=bool(slope is not None and slope>0))
        records[str(row)]=rec
    return dict(quotient_monodromy=matrix_record(monodromy),nonflow_quotient_schur=core_schur,
                exact_flow_division_remainder=ball_record(remainder),
                quotient_stable=bool(remainder.contains(0) and core_schur["passed"]),transverse=records)


def full_sensitivity_spectrum(q,derivative,rho,r,args,audit=None):
    x,ode=integrate_sensitivity(q,rho,r,args,derivative,True,audit)
    core=paired.extract(x,ode.offsets["core"])
    quarters,derivatives={},{}
    for j in range(4):
        six=paired.extract(x,ode.offsets[str(j)])
        quarters[j]=arb_mat([[six[a,b] for b in range(3)] for a in range(3)])
        derivatives[j]=arb_mat([[six[3+a,b] for b in range(3)] for a in range(3)])
    return spectral_analysis(core,quarters,derivatives)


def plain_spectrum(q,rho,r,args,audit):
    ode=paired.PairedODE(rho,r,with_transverse=True,order=args.order)
    x=paired.flow(ode,paired.initial_state(q),q[-1],args.steps,args.order,audit)
    return spectral_analysis(paired.extract(x,ode.offsets["core"]),
                             {j:paired.extract(x,ode.offsets[str(j)]) for j in range(4)})


def distinctions(q,rho,r):
    state=paired.initial_state(q)
    xs=state[8:]+[-sum(state[8:],ZERO)]
    tests=[]
    for a in range(4):
        for b in range(a+1,4):
            differences=[state[a]-state[b],state[4+a]-state[4+b],xs[a]-xs[b]]
            tests.append(dict(rows=[a,b],differences=[ball_record(v) for v in differences],
                              passed=bool(any(not v.contains(0) for v in differences))))
    speed=paired.PairedODE(rho,r,with_core=False).value(state,ONE)[8]
    return tests,speed


def validate_box(center,rho,r,radius,args,audit,centered_parameter=True):
    center_f,center_j=paired.shooting(center,rho.mid() if centered_parameter else rho,r,args.steps,args.order,audit)
    parameter_record=None
    if centered_parameter:
        x,_=integrate_sensitivity(center,rho,r,args,audit=audit)
        partial=arb_mat([[x[N+i]] for i in range(N)])
        parameter_record=matrix_record(partial)
        center_f+=(rho-rho.mid())*partial
    box=[arb(v,radius) for v in center]
    _,box_j=paired.shooting(box,rho,r,args.steps,args.order,audit)
    c=midpoint_matrix(midpoint_matrix(center_j).inv())
    defect=paired.identity(N)-c*box_j
    image=-c*center_f+defect*arb_mat([[arb(0,radius)] for _ in range(N)])
    ratios=[image[i,0].abs_upper()/radius for i in range(N)]
    contraction=max(sum((defect[i,j].abs_upper() for j in range(N)),ZERO).upper() for i in range(N))
    distinct,speed=distinctions(box,rho,r)
    passed=(all(v<1 for v in ratios) and contraction<1 and not c.det().contains(0)
            and box[-1]>0 and speed>0 and all(t["passed"] for t in distinct))
    rec=dict(shooting_center=[ball_record(v) for v in center],shooting_box=[ball_record(v) for v in box],
             center_residual=matrix_record(center_f),interval_jacobian=matrix_record(box_j),
             parameter_center_derivative=parameter_record,preconditioner=matrix_record(c),
             krawczyk_ratios=[ball_record(v) for v in ratios],krawczyk_contraction=ball_record(contraction),
             cluster_distinctions=distinct,section_speed=ball_record(speed),passed=bool(passed),
             positive_period=ball_record(4*box[-1]))
    return box,box_j,rec


def newton(center,rho,r,args,count,history):
    jac=None
    for k in range(count):
        f,jac=paired.shooting(center,rho,r,args.steps,args.order)
        correction=midpoint_matrix(jac).inv()*midpoint_matrix(f)
        history.append(dict(iteration=k,residual=matrix_record(f)))
        center=[(v-correction[i,0]).mid() for i,v in enumerate(center)]
    return center,jac


def run(args,record):
    ctx.prec=args.bits
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    seed_file=ROOT/"results/motor-postcritical/paired-second-crossing/analysis.json"
    seed=json.loads(seed_file.read_text())
    center=[arb(str(seed["state"][j])).mid() for j in paired.FREE]+[(arb(str(seed["period"]))/4).mid()]
    rho=arb(args.rho_center).mid()
    proposal_settings=SimpleNamespace(steps=args.proposal_steps,order=args.proposal_order)
    history=[]
    record["center_proposals"]=history
    if args.center_input:
        saved=json.loads(args.center_input.read_text())
        if "family" in saved:
            center=[read_ball(v).mid() for v in saved["family"]["shooting_center"]]
            rho=arb(saved["exact_center_decimal"]).mid()
        else:
            prefix="latest_proposed" if "latest_proposed_center" in saved else "in_progress"
            center=[read_ball(v).mid() for v in saved[prefix+"_center"]]
            rho=read_ball(saved[prefix+"_rho"]).mid()
        record["center_input_sha256"]=hashlib.sha256(args.center_input.read_bytes()).hexdigest()
    for iteration in range(args.refine_iterations):
        nh=[]
        center,jac=newton(center,rho,r,proposal_settings,2,nh)
        record["in_progress_center"]=[ball_record(v) for v in center]
        record["in_progress_rho"]=ball_record(rho)
        record["in_progress_newton"]=nh
        derivative,_=implicit_derivative(center,rho,r,jac,proposal_settings,None)
        spectrum=full_sensitivity_spectrum(center,derivative,rho,r,proposal_settings)
        g=read_ball(spectrum["transverse"]["0"]["determinant"])
        dg=read_ball(spectrum["transverse"]["0"]["determinant_derivative"])
        change=(g/dg).mid()
        history.append(dict(iteration=iteration,rho=ball_record(rho),newton=nh,
                            determinant=ball_record(g),proposed_derivative=ball_record(dg),
                            qualification="Center proposal only; not a certified periodic branch derivative"))
        rho=(rho-change).mid()
        center=[(v-change*derivative[i]).mid() for i,v in enumerate(center)]
        record["latest_proposed_center"]=[ball_record(v) for v in center]
        record["latest_proposed_rho"]=ball_record(rho)
        print("Proposed rho",rho.str(65,radius=False),flush=True)
    center_text=rho.str(65,radius=False)
    parameter=arb(arb(center_text),arb(args.rho_radius))
    with localcontext() as dc:
        dc.prec=180
        endpoints=[str(Decimal(center_text)+sgn*Decimal(args.rho_radius)) for sgn in (-1,1)]
    nh=[]
    center,_=newton(center,parameter.mid(),r,args,2,nh)
    audit=[]
    record.update(seed_file=str(seed_file.relative_to(ROOT)),seed_sha256=hashlib.sha256(seed_file.read_bytes()).hexdigest(),
                  center_proposals=history,final_center_newton=nh,exact_center_decimal=center_text,
                  exact_endpoint_decimals=endpoints,parameter_r=ball_record(r),parameter_rho=ball_record(parameter),
                  integrations=audit)
    box,jac,family=validate_box(center,parameter,r,arb(args.radius).upper(),args,audit)
    record["family"]=family
    if not family["passed"]:
        record["status"]="family_not_certified"
        return
    print("Uniform four-pair family certified",flush=True)
    derivative,sensitivity=implicit_derivative(box,parameter,r,jac,args,audit)
    record["sensitivity"]=sensitivity
    spectrum=full_sensitivity_spectrum(box,derivative,parameter,r,args,audit)
    record["spectrum"]=spectrum
    gates=dict(family=family["passed"],quotient_stable=spectrum["quotient_stable"])
    for name,transverse in spectrum["transverse"].items():
        gates[name+"_simple_roots"]=transverse["simple_roots"]
        print("Row",name,"g'",transverse["determinant_derivative"]["display"],flush=True)
    endpoints_records={}
    record["endpoints"]=endpoints_records
    for name,value,expected in (("low",endpoints[0],-1),("high",endpoints[1],1)):
        fixed=arb(value)
        nh=[]
        guess,_=newton(center,fixed,r,args,1,nh)
        small,_,cert=validate_box(guess,fixed,r,arb(args.endpoint_radius).upper(),args,audit,False)
        cert.update(parameter_rho=ball_record(fixed),parameter_r=ball_record(r),newton_proposals=nh)
        cert["spectrum"]=plain_spectrum(small,fixed,r,args,audit)
        endpoints_records[name]=cert
        gates[name+"_periodic"]=cert["passed"]
        gates[name+"_same_shooting_branch"]=all(exact_contains(b,x,strict=True)
             for b,x in zip(family["shooting_box"],cert["shooting_box"]))
        gates[name+"_rho_inside_family"]=exact_contains(record["parameter_rho"],cert["parameter_rho"])
        gates[name+"_quotient_stable"]=cert["spectrum"]["quotient_stable"]
        for j,transverse in cert["spectrum"]["transverse"].items():
            gates[name+"_"+j+"_roots"]=transverse["simple_roots"]
            gates[name+"_"+j+"_sign"]=transverse["determinant_sign"]==expected
        print("Endpoint",name,"periodic",cert["passed"],flush=True)
    record["gates"]={k:bool(v) for k,v in gates.items()}
    exists=all(gates.values())
    monotone=all(v["determinant_strictly_increasing"] and v["multiplier_crossing_speed_positive"]
                 for v in spectrum["transverse"].values())
    record["unique_crossing_gate"]=bool(monotone)
    record["status"]=("unique_fourfold_transverse_crossing_certified" if exists and monotone else
                      "fourfold_transverse_crossing_exists_certified" if exists else "not_certified")
    record["statement"]="A validated continuous paired periodic family has an interior transverse multiplier -1 of exactly fourfold algebraic and geometric multiplicity. Uniqueness and nonzero crossing speed require the separate monotonicity gate."
    record["not_proved"]=["Nonlinear branch birth or period-doubling theorem","Global connection to the first paired-cycle certificate",
                          "Chaotic attractor","Scientific priority"]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho-center",default="3.2748747165675427")
    parser.add_argument("--rho-radius",default="1e-53")
    parser.add_argument("--radius",default="1e-50")
    parser.add_argument("--endpoint-radius",default="1e-75")
    parser.add_argument("--bits",type=int,default=512)
    parser.add_argument("--steps",type=int,default=288)
    parser.add_argument("--order",type=int,default=64)
    parser.add_argument("--proposal-steps",type=int,default=256)
    parser.add_argument("--proposal-order",type=int,default=56)
    parser.add_argument("--refine-iterations",type=int,default=3)
    parser.add_argument("--center-input",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if (args.bits<64 or min(args.steps,args.proposal_steps)<1 or min(args.order,args.proposal_order)<2
            or args.refine_iterations<0 or not all(arb(v)>0 for v in (args.rho_radius,args.radius,args.endpoint_radius))):
        parser.error("Positive radii, steps>=1, orders>=2, bits>=64 and nonnegative refinements required")
    if args.output.exists(): raise FileExistsError("Choose unused output to preserve evidence")
    paths=[Path(__file__),Path(base.__file__),Path(paired.__file__),ROOT/"research/validate_motor_crossing.py",
           ROOT/"research/verify_motor_certificate.py"]
    sources={p.name:p.read_bytes() for p in paths}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                arguments={k:str(v) for k,v in vars(args).items()},python=platform.python_version(),python_flint=flint.__version__,
                source_sha256={n:hashlib.sha256(b).hexdigest() for n,b in sources.items()})
    began=time.monotonic()
    try:
        run(args,record)
    except Exception as error:
        record.update(status="error",error=repr(error))
        raise
    finally:
        record["wall_seconds"]=time.monotonic()-began
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
        for name,data in sources.items(): args.output.with_name(args.output.stem+"."+name).write_bytes(data)
        print(record.get("status"),args.output,flush=True)
    if record.get("status") not in ("unique_fourfold_transverse_crossing_certified","fourfold_transverse_crossing_exists_certified"):
        raise SystemExit(2)


if __name__=="__main__":
    main()
