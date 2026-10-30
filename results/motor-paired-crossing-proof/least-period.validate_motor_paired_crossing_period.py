"""Certify least period uniformly on the second paired crossing family.

This separately revalidates a larger shooting box and links the archived
crossing family by exact dyadic containment. Validated tubes count all X0
zeros, including the periodic boundary ownership, as in the paired theorem.
"""
import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

from flint import arb, ctx
import validate_motor_cycle as base
import validate_motor_paired_cycle as paired
import validate_motor_paired_crossing as crossing
import validate_motor_paired_period as period
from validate_motor_cycle import ROOT, ZERO, ball_record, read_ball
from verify_motor_certificate import exact_contains


def prove(args,record):
    ctx.prec=args.bits
    source=json.loads(args.input.read_text())
    if source.get("status") not in ("unique_fourfold_transverse_crossing_certified",
                                    "fourfold_transverse_crossing_exists_certified"):
        raise ValueError("Input must be a successful crossing-family certificate")
    rho=read_ball(source["parameter_rho"])
    r=1+arb("19/20")/arb("32001/500000").sqrt()
    center=[read_ball(v).mid() for v in source["family"]["shooting_center"]]
    audit=[]
    box,_,anchor=crossing.validate_box(center,rho,r,arb(args.radius).upper(),args,audit)
    linked=all(exact_contains(b,old,strict=True) for b,old in
               zip(anchor["shooting_box"],source["family"]["shooting_box"]))
    linked=linked and exact_contains(ball_record(rho),source["parameter_rho"])
    record.update(anchor_revalidation=anchor,source_family_contained=bool(linked),
                  parameter_rho=ball_record(rho),parameter_r=ball_record(r),integration=audit)
    if not anchor["passed"] or not linked:
        record["status"]="anchor_not_certified"
        return
    initial=paired.initial_state(box)
    ode=period.TubeRecorder(rho,r,args.order)
    endpoint=paired.flow(ode,initial,box[-1],args.steps,args.order,audit)
    if len(ode.tubes)!=args.steps:
        raise ArithmeticError("Expected one validated Picard tube per subdivision")
    initial_x=period.row_positions(initial)
    positions=[period.row_positions(tube) for tube in ode.tubes]
    velocities=[[-rho*x+r*(tube[4+j]-sum(tube[4:8],ZERO)/4)
                 for j,x in enumerate(period.row_positions(tube))] for tube in ode.tubes]
    windows=[]
    for row in (0,2,1,3):
        signs=[period.sign(p[row]) for p in positions]
        i=0
        while i<args.steps:
            if signs[i]!=0:
                i+=1
                continue
            start=i
            while i<args.steps and signs[i]==0:
                i+=1
            end=i
            left=period.sign(initial_x[row]) if start==0 else signs[start-1]
            right=period.sign(initial_x[paired.PERMUTATION[row]]) if end==args.steps else signs[end]
            directions=[period.sign(v[row]) for v in velocities[start:end]]
            direction=1 if all(v==1 for v in directions) else -1 if all(v==-1 for v in directions) else 0
            classification=period.classify_window(direction,left,right,row,start==0,end==args.steps)
            windows.append(dict(row=row,start_step=start,end_step_exclusive=end,
                normalized_interval=[f"{start}/{args.steps}",f"{end}/{args.steps}"],
                left_sign=left,right_sign=right,direction=direction,**classification))
    total_up=sum(w["upward"] for w in windows)
    total_down=sum(w["downward"] for w in windows)
    gates=dict(anchor=anchor["passed"],source_family_contained=linked,
               initial_phase_is_exact_zero=initial_x[0].is_zero(),
               other_initial_positions_nonzero=all(not x.contains(0) for x in initial_x[1:]),
               all_possible_zero_windows_resolved=all(w["passed"] for w in windows),
               exactly_two_upward_zeros=total_up==2,exactly_two_downward_zeros=total_down==2,
               owns_initial_zero_once=sum(w["reason"]=="include the initial upward zero at total time zero" for w in windows)==1,
               excludes_final_zero_once=sum(w["reason"]=="exclude the final upward zero at total time T" for w in windows)==1,
               half_period_X0_strictly_positive=initial_x[1]>0,
               endpoint_encloses_exact_final_zero=period.row_positions(endpoint)[3].contains(0))
    record.update(status="least_second_crossing_family_period_certified" if all(gates.values()) else "not_certified",
                  gates={k:bool(v) for k,v in gates.items()},positive_period=anchor["positive_period"],
                  initial_positions=[ball_record(v) for v in initial_x],
                  half_period_X0=ball_record(initial_x[1]),total_upward_zeros=total_up,
                  total_downward_zeros=total_down,quarter_row_order=[0,2,1,3],windows=windows,
                  tubes=[dict(step=i,X=[ball_record(v) for v in p],
                              V=[ball_record(v) for v in velocities[i]]) for i,p in enumerate(positions)],
                  theorem="Uniformly on the linked family, exactly two upward X0 zeros occur in [0,T). A smaller state period would repeat this count, permitting only T/2. Exact quarter symmetry gives X0(T/2)=X1(0)>0, whereas X0(0)=0, excluding T/2. Thus T is least.",
                  limitations=["No nonlinear flip bifurcation theorem or new branch existence",
                               "Same interval integrator and reviewed zero-window classification as the paired-cycle proof"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--bits",type=int,default=576)
    parser.add_argument("--steps",type=int,default=320)
    parser.add_argument("--order",type=int,default=68)
    parser.add_argument("--radius",default="1e-48")
    args=parser.parse_args()
    if args.bits<64 or args.steps<4 or args.order<2 or not arb(args.radius)>0:
        parser.error("Require bits>=64, steps>=4, order>=2 and positive radius")
    if args.output.exists(): raise FileExistsError("Choose an unused output")
    paths=[Path(__file__),Path(base.__file__),Path(paired.__file__),Path(crossing.__file__),Path(period.__file__),
           ROOT/"research/validate_motor_crossing.py",ROOT/"research/verify_motor_certificate.py"]
    sources={p.name:p.read_bytes() for p in paths}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                arguments={k:str(v) for k,v in vars(args).items()},
                input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                source_sha256={n:hashlib.sha256(b).hexdigest() for n,b in sources.items()})
    began=time.monotonic()
    try:
        prove(args,record)
    except Exception as error:
        record.update(status="error",error=repr(error))
        raise
    finally:
        record["wall_seconds"]=time.monotonic()-began
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
        for name,data in sources.items(): args.output.with_name(args.output.stem+"."+name).write_bytes(data)
        print(record.get("status"),args.output,flush=True)
    if record.get("status")!="least_second_crossing_family_period_certified": raise SystemExit(2)


if __name__=="__main__": main()
