"""Validate the least period of the four-pair rotating motor orbit.

The shooting root is freshly revalidated. Picard tubes count exactly two
upward X0 zeros over four quarters, with explicit endpoint ownership.
The half-period permutation excludes a period T/2, hence T is least.
"""
import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace

from flint import arb, ctx
import validate_motor_cycle as base
import validate_motor_paired_cycle as paired
from validate_motor_cycle import ROOT, ZERO, ball_record, read_ball


def sign(value):
    return 1 if value > 0 else -1 if value < 0 else 0


def classify_window(direction,left,right,row,starts,ends):
    if direction not in (-1,1):
        return dict(passed=False,reason="derivative sign not strict",upward=0,downward=0)
    if left==right and left!=0:
        return dict(passed=True,reason="same strict endpoint signs and monotonicity",upward=0,downward=0)
    if right-left==2*direction:
        return dict(passed=True,reason="unique interior zero by opposite signs and monotonicity",
                    upward=int(direction==1),downward=int(direction==-1))
    if left==0 and right==1 and direction==1 and starts and row==0:
        return dict(passed=True,reason="include the initial upward zero at total time zero",upward=1,downward=0)
    if right==0 and left==-1 and direction==1 and ends and row==3:
        return dict(passed=True,reason="exclude the final upward zero at total time T",upward=0,downward=0)
    return dict(passed=False,reason="endpoint signs do not certify a permitted zero count",upward=0,downward=0)


class TubeRecorder(paired.PairedODE):
    def __init__(self,rho,r,order):
        super().__init__(rho,r,with_core=False,order=order)
        self.tubes=[]

    def jets(self,x,duration,order):
        if order==self.order+1:
            self.tubes.append(list(x))
        return super().jets(x,duration,order)


def row_positions(state):
    return list(state[8:11])+[-sum(state[8:11],ZERO)]


def prove(args,record):
    ctx.prec=args.bits
    source=json.loads(args.input.read_text())
    anchor={}
    paired.prove(SimpleNamespace(bits=args.bits,steps=args.steps,order=args.order,
        radius=args.radius,newton=0,verify_input=args.input,rho=source["exact_rho_decimal"],
        existence_only=True,generic_jets=False),anchor)
    record["anchor_revalidation"]=anchor
    if anchor.get("status")!="periodic_four_pair_orbit_certified":
        record["status"]="anchor_not_certified"
        return
    box=[read_ball(v) for v in anchor["shooting_box"]]
    initial=paired.initial_state(box)
    rho,r=read_ball(anchor["parameter_rho"]),read_ball(anchor["parameter_r"])
    ode=TubeRecorder(rho,r,args.order)
    audit=[]
    endpoint=paired.flow(ode,initial,box[-1],args.steps,args.order,audit)
    if len(ode.tubes)!=args.steps:
        raise ArithmeticError("Expected one validated Picard tube per subdivision")
    initial_x=row_positions(initial)
    positions=[row_positions(tube) for tube in ode.tubes]
    velocities=[[-rho*x+r*(tube[4+j]-sum(tube[4:8],ZERO)/4)
                 for j,x in enumerate(row_positions(tube))] for tube in ode.tubes]
    windows=[]
    for row in (0,2,1,3):
        signs=[sign(p[row]) for p in positions]
        i=0
        while i<args.steps:
            if signs[i]!=0:
                i+=1
                continue
            start=i
            while i<args.steps and signs[i]==0:
                i+=1
            end=i
            left=sign(initial_x[row]) if start==0 else signs[start-1]
            right=sign(initial_x[paired.PERMUTATION[row]]) if end==args.steps else signs[end]
            directions=[sign(v[row]) for v in velocities[start:end]]
            direction=1 if all(v==1 for v in directions) else -1 if all(v==-1 for v in directions) else 0
            classification=classify_window(direction,left,right,row,start==0,end==args.steps)
            windows.append(dict(row=row,start_step=start,end_step_exclusive=end,
                normalized_interval=[f"{start}/{args.steps}",f"{end}/{args.steps}"],
                left_sign=left,right_sign=right,direction=direction,**classification))
    total_up=sum(w["upward"] for w in windows)
    total_down=sum(w["downward"] for w in windows)
    gates=dict(anchor=anchor["krawczyk_passed"],
               initial_phase_is_exact_zero=initial_x[0].is_zero(),
               other_initial_positions_nonzero=all(not x.contains(0) for x in initial_x[1:]),
               all_possible_zero_windows_resolved=all(w["passed"] for w in windows),
               exactly_two_upward_zeros=total_up==2,
               exactly_two_downward_zeros=total_down==2,
               owns_initial_zero_once=sum(w["reason"]=="include the initial upward zero at total time zero" for w in windows)==1,
               excludes_final_zero_once=sum(w["reason"]=="exclude the final upward zero at total time T" for w in windows)==1,
               half_period_X0_strictly_positive=initial_x[1]>0,
               endpoint_encloses_exact_final_zero=row_positions(endpoint)[3].contains(0))
    record.update(status="least_four_pair_period_certified" if all(gates.values()) else "not_certified",
                  gates={k:bool(v) for k,v in gates.items()},positive_period=anchor["positive_period"],
                  initial_positions=[ball_record(v) for v in initial_x],
                  half_period_X0=ball_record(initial_x[1]),total_upward_zeros=total_up,
                  total_downward_zeros=total_down,quarter_row_order=[0,2,1,3],windows=windows,
                  integration=audit,tubes=[dict(step=i,X=[ball_record(v) for v in p],
                                               V=[ball_record(v) for v in velocities[i]])
                                            for i,p in enumerate(positions)],
                  theorem="There are exactly two upward X0 zeros in [0,T). Any smaller state period would repeat that count, so only T/2 is possible. Exact quarter symmetry gives X0(T/2)=X1(0)>0 while X0(0)=0, excluding T/2. Therefore T is the least positive period.",
                  limitations=["No nonlinear bifurcation or chaotic-attractor proof",
                               "The same interval integrator supplies the tube enclosures"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--bits",type=int,default=448)
    parser.add_argument("--steps",type=int,default=288)
    parser.add_argument("--order",type=int,default=42)
    parser.add_argument("--radius",default="1e-38")
    args=parser.parse_args()
    if args.bits<64 or args.steps<4 or args.order<2 or not arb(args.radius)>0:
        parser.error("Require bits>=64, steps>=4, order>=2 and positive radius")
    if args.output.exists():
        raise FileExistsError("Choose an unused output path")
    sources={p.name:p.read_bytes() for p in (Path(__file__),Path(base.__file__),Path(paired.__file__))}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                arguments={k:str(v) for k,v in vars(args).items()},input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                source_sha256={n:hashlib.sha256(b).hexdigest() for n,b in sources.items()})
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
        for name,data in sources.items():
            args.output.with_name(args.output.stem+"."+name).write_bytes(data)
        print(record.get("status"),args.output,flush=True)
    if record.get("status")!="least_four_pair_period_certified":
        raise SystemExit(2)


if __name__=="__main__":
    main()
