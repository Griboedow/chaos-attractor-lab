"""Assemble and recheck the nine rigorously computed transfers and dimension bounds."""
import argparse
import base64
import hashlib
import json
import zlib
from pathlib import Path
from fractions import Fraction
from entropy import INPUT,FACES,SIDES,BRIDGES,regions
from affine_arb import region
from affine_directed import replay
from dimension_lower import bound as global_bound
from phase_dimension import bound as phase_bound
from upper import verify as verify_upper
from upper_directed import replay as replay_upper
from flint import arb,ctx
from verify_boxes import column,frames,rotate,values

def get(f):
    return json.loads(f.read_text())

def write(f,x):
    f.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")

def digest(o):
    return hashlib.sha256(json.dumps(o,sort_keys=True).encode()).hexdigest()

def controls(box):
    ctx.prec=256;o,M,V=frames(box);o,M,V=o[0],M[0],V[0];p=[arb(float(v)) for v in box["params"]]
    out=[]
    for i,t in enumerate(FACES):
        x=o+M*column([0.,0.,t])
        for _ in range(20):x,_=rotate(values(x),p)
        u=(V*(x-o))[2,0]
        if not (u<FACES[0] if SIDES[i]<0 else u>FACES[-1]):raise ValueError("A face-center control failed")
        out.append(dict(face=t,exact_input_normalized=[0.,0.,t],output_u_approx=float(u.mid()),
            opposite_crossing_ruled_out=True))
    return out

def assemble(cache,out,bits):
    out.mkdir(parents=True,exist_ok=True);box=get(INPUT)["box"];tests=regions();proof=[]
    for t in tests:
        key=t["name"]+".json"
        a=get(cache/f"affine-{bits}"/key);d=get(cache/f"affine-directed-{bits}"/key)
        if a["region"]!=t or d["name"]!=t["name"] or a["bits"]!=bits or a["box_sha256"]!=digest(box):
            raise ValueError("Mismatched region")
        proof.append(dict(arb=a,directed_binary64=d))
    entropy=dict(case=box["case"],base_order=box["case_order"],symbolic_order=2*box["case_order"],
        radius_each_of_ten_binary64_coefficients=f"2^-{bits}",block_radius_bits=bits,source_trap_box=box,
        source_box_sha256=digest(box),axis_zero_based=2,faces=FACES,orientations=SIDES,intermediate_boxes=BRIDGES,
        complete_region_checks=proof,strict_trapping_compositions=2,contained_full_slabs=2,
        entropy_lower="log(2)/20",face_center_negative_controls=controls(box))
    write(out/"entropy-certificate.json",entropy)
    low=get(cache/f"dimension-lower-{bits}.json");phase=get(cache/f"phase-dimension-{bits}.json")
    candidates=[cache/f"recenter-upper-{bits}-17.json",*cache.glob("upper-*-*.json")]
    suitable=[]
    for f in candidates:
        if not f.exists():continue
        a=get(f);d=f.with_name(f.stem+"-directed.json")
        if a.get("block_radius_bits")==bits and "partition_zlib_base64" in a and d.exists():
            replayed=get(d)
            if replayed.get("dimension_upper")==a.get("dimension_upper") and replayed.get("maximum_power_upper",1)<1:
                suitable.append((Fraction(a["dimension_upper"]),a,replayed))
    _,upper,second=min(suitable,key=lambda x:x[0])
    dimension=dict(case=box["case"],block_radius_bits=bits,global_lower=low,phase_lower=phase,
        upper_Arb=upper,upper_directed=second,conservative_lower=phase["conservative_dimension_lower"],
        upper=upper["dimension_upper"])
    write(out/"dimension-certificate.json",dimension)
    report=dict(case=box["case"],block_radius_bits=bits,components=10,topological_dimension_per_component=1,
        entropy_lower="log(2)/20",strict_Hausdorff_dimension_lower=dimension["conservative_lower"],
        Hausdorff_dimension_upper=dimension["upper"],files={f.name:hashlib.sha256(f.read_bytes()).hexdigest()
        for f in [out/"entropy-certificate.json",out/"dimension-certificate.json"]})
    write(out/"manifest.json",report)
    return report

def check(out):
    e=get(out/"entropy-certificate.json");d=get(out/"dimension-certificate.json");summary=get(out/"manifest.json")
    box=e["source_trap_box"];bits=e["block_radius_bits"]
    if box!=get(INPUT)["box"] or e["source_box_sha256"]!=digest(box):raise ValueError("Different original model")
    if (e["faces"],e["orientations"],e["intermediate_boxes"],e["axis_zero_based"],e["base_order"],e["symbolic_order"],
        e["strict_trapping_compositions"],e["contained_full_slabs"])!=(FACES,SIDES,BRIDGES,2,10,20,2,2):raise ValueError("Geometry mismatch")
    if len(e["complete_region_checks"])!=9 or bits!=d["block_radius_bits"] or bits!=summary["block_radius_bits"]:raise ValueError("Wrong evidence scope")
    for name,checksum in summary["files"].items():
        if hashlib.sha256((out/name).read_bytes()).hexdigest()!=checksum:raise ValueError("Packet checksum mismatch")
    for t,checkset in zip(regions(),e["complete_region_checks"]):
        a=checkset["arb"];b=checkset["directed_binary64"]
        other=region(box,t,bits,80000)
        if (other["region"]!=a["region"] or other["bits"]!=bits or a["box_sha256"]!=digest(box) or
            not 0<a["normalized_strict_margin_above"]<=other["normalized_strict_margin_above"]):raise ValueError("Arb transfer mismatch")
        result=replay(box,t,a,bits)
        if result["name"]!=b["name"] or result["normalized_margin_lower_directed"]<=0:raise ValueError("Directed transfer failed")
    if not all(item["opposite_crossing_ruled_out"] for item in controls(box)):raise ValueError("Control failure")
    if global_bound(box,bits)["dimension_lower"]<d["global_lower"]["dimension_lower"]:raise ValueError("Global dimension lower mismatch")
    tree=e["complete_region_checks"][-1]["arb"]
    phase=phase_bound(box,tree,bits,4)
    if phase["conservative_dimension_lower"]<summary["strict_Hausdorff_dimension_lower"]:raise ValueError("Phase dimension lower mismatch")
    a=d["upper_Arb"];other=verify_upper(box,bits,a["dimension_upper"],100000)
    if other["dimension_upper"]!=summary["Hausdorff_dimension_upper"]:raise ValueError("Upper dimension mismatch")
    if replay_upper(box,a)["maximum_power_upper"]>=1:raise ValueError("Directed dimension replay failed")
    return summary

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest="action",required=True)
    build=sub.add_parser("assemble");build.add_argument("cache",type=Path);build.add_argument("out",type=Path);build.add_argument("--bits",type=int,default=23)
    v=sub.add_parser("verify");v.add_argument("out",type=Path)
    args=ap.parse_args()
    print(json.dumps(assemble(args.cache,args.out,args.bits) if args.action=="assemble" else check(args.out),indent=2),flush=True)
