"""Validate the saddle's least period and the period of one synchronized row.

Fresh multiple shooting links the root to an existing spectral certificate.
Validated Picard tubes then count strict upward section crossings.
"""
import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace
from flint import ctx
import validate_motor_cycle as base
import validate_motor_paired_cycle as paired
import validate_motor_saddle_cycle as saddle
import verify_motor_certificate as verifier
from validate_motor_cycle import ball_record,read_ball,ZERO
from verify_motor_certificate import exact_bounds

def sign(x):
    return 1 if x>0 else -1 if x<0 else 0

class Recorder(saddle.SaddleODE):
    def __init__(self,rho,r,order):
        super().__init__(rho,r,with_core=False,order=order)
        self.tubes=[]

    def jets(self,x,duration,order):
        if order==self.order+1: self.tubes.append(list(x))
        return super().jets(x,duration,order)

def count_windows(row,positions,velocities,start_x,end_x,own_initial=False,exclude_final=False):
    signs=[sign(x[row]) for x in positions]
    windows=[]
    i=0
    while i<len(signs):
        if signs[i]!=0:
            i+=1
            continue
        start=i
        while i<len(signs) and signs[i]==0: i+=1
        end=i
        left=sign(start_x) if start==0 else signs[start-1]
        right=sign(end_x) if end==len(signs) else signs[end]
        ds=[sign(x[row]) for x in velocities[start:end]]
        direction=1 if all(v==1 for v in ds) else -1 if all(v==-1 for v in ds) else 0
        up=down=0
        passed=False
        reason="unresolved zero window"
        if direction:
            if left==right and left!=0:
                passed=True;reason="same strict endpoint signs and monotonicity"
            elif right-left==2*direction:
                passed=True;reason="one interior zero by opposite signs and monotonicity"
                up=int(direction==1);down=int(direction==-1)
            elif own_initial and start==0 and left==0 and right==1 and direction==1:
                passed=True;reason="owned initial upward zero";up=1
            elif exclude_final and end==len(signs) and left==-1 and right==0 and direction==1:
                passed=True;reason="excluded final upward zero"
        windows.append(dict(row=row,start_step=start,end_step_exclusive=end,left_sign=left,
            right_sign=right,derivative_sign=direction,passed=passed,reason=reason,upward=up,downward=down))
    return dict(row=row,windows=windows,passed=all(x["passed"] for x in windows),
        upward=sum(x["upward"] for x in windows),downward=sum(x["downward"] for x in windows))

def prove(args,record):
    source=json.loads(args.input.read_text())
    if source.get("status")!="six_cluster_saddle_certified":
        raise ValueError("Input must have a successful complete saddle certificate")
    segments=source["shooting_segments"]
    if segments!=4: raise ValueError("This period certificate uses four segments per half period")
    a=source["arguments"]
    ctx.prec=int(a["bits"])
    anchor={}
    saddle.prove(SimpleNamespace(bits=ctx.prec,rho=source["exact_rho_decimal"],
        steps=int(a["steps"]),order=int(a["order"]),segments=segments,radius=a["radius"],
        newton=0,verify_input=args.input,existence_only=True,generic_jets=False),anchor)
    record["anchor_revalidation"]=anchor
    if not anchor.get("krawczyk_passed"):
        record["status"]="anchor_not_certified"
        return
    contained=all(exact_bounds(x)[0]<=exact_bounds(y)[0] and exact_bounds(y)[1]<=exact_bounds(x)[1]
        for x,y in zip(source["shooting_box"],anchor["shooting_box"]))
    if not contained:
        record.update(status="anchor_link_not_certified",anchor_box_contained=False)
        return
    box=[read_ball(x) for x in anchor["shooting_box"]]
    rho,r=read_ball(anchor["parameter_rho"]),read_ball(anchor["parameter_r"])
    initial_states=saddle.segment_states(box,segments)
    positions=[];velocities=[];audit=[]
    for state in initial_states:
        ode=Recorder(rho,r,int(a["order"]))
        saddle.flow(ode,state,box[-1],int(a["steps"]),int(a["order"]),audit)
        if len(ode.tubes)!=int(a["steps"]): raise ArithmeticError("Missing validated tube")
        for tube in ode.tubes:
            x=saddle.positions(tube)
            mean=sum((w*tube[saddle.K+j] for j,w in enumerate(saddle.weights())),ZERO)
            v=[-rho*x[j]+r*(tube[saddle.K+j]-mean) for j in range(saddle.K)]
            positions.append(x);velocities.append(v)
    x0=saddle.positions(initial_states[0])
    xquarter=saddle.positions(initial_states[2])
    full_parts=[count_windows(0,positions,velocities,x0[0],x0[1],own_initial=True),
        count_windows(1,positions,velocities,x0[1],x0[0],exclude_final=True)]
    row2=count_windows(2,positions,velocities,x0[2],x0[2])
    up=sum(x["upward"] for x in full_parts)
    down=sum(x["downward"] for x in full_parts)
    gates=dict(anchor_box_contained=contained,initial_section_exact_zero=x0[0].is_zero(),
        full_zero_windows_resolved=all(x["passed"] for x in full_parts),
        full_four_upward=up==4,full_four_downward=down==4,
        owns_initial_once=sum(w["reason"]=="owned initial upward zero" for part in full_parts for w in part["windows"])==1,
        excludes_final_once=sum(w["reason"]=="excluded final upward zero" for part in full_parts for w in part["windows"])==1,
        half_period_not_state_period=not x0[1].contains(0),
        row2_initial_position_nonzero=not x0[2].contains(0),row2_windows_resolved=row2["passed"],
        row2_two_upward_in_half=row2["upward"]==2,row2_two_downward_in_half=row2["downward"]==2,
        row2_quarter_difference_nonzero=not (xquarter[2]-x0[2]).contains(0))
    record.update(status="least_saddle_and_row_periods_certified" if all(gates.values()) else "not_certified",
        gates={k:bool(v) for k,v in gates.items()},positive_period=anchor["positive_period"],
        row2_least_period=ball_record(4*box[-1]),initial_positions=[ball_record(x) for x in x0],
        quarter_period_row2_difference=ball_record(xquarter[2]-x0[2]),full_crossings=full_parts,
        row2_crossings=row2,integration=audit,
        tubes=[dict(step=i,X=[ball_record(x) for x in xs],V=[ball_record(v) for v in velocities[i]]) for i,xs in enumerate(positions)],
        theorem=["X0 has exactly four upward zeros on [0,T). Any smaller state period divides the count; T/2 is not a period, also excluding T/4. Thus T is least.",
            "Row2 is fixed by the half-period permutation, so its X component is T/2 periodic. It has exactly two upward zeros there; its nonzero change at T/4 excludes a smaller period. Thus X2 and the complete row2 state have least period T/2."],
        limitations=["No heteroclinic connection or attraction theorem","Uses the same interval integrator and source spectral certificate"])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError("Preserve evidence: choose a new output")
    sources={p.name:p.read_bytes() for p in (Path(__file__),Path(base.__file__),Path(paired.__file__),Path(saddle.__file__),Path(verifier.__file__),Path(saddle.spectral.__file__))}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        input_path=str(args.input),input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        source_sha256={n:hashlib.sha256(b).hexdigest() for n,b in sources.items()})
    start=time.monotonic()
    try: prove(args,record)
    except Exception as error:
        record.update(status="error",error=repr(error))
        raise
    finally:
        record["wall_seconds"]=time.monotonic()-start
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
        for n,b in sources.items(): args.output.with_name(args.output.stem+"."+n).write_bytes(b)
        print(record.get("status"),args.output,flush=True)
    if record.get("status")!="least_saddle_and_row_periods_certified": raise SystemExit(2)

if __name__=="__main__": main()
