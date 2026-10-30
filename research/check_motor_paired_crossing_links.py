"""Audit exact containment and immutable inputs of the paired crossing proofs."""
import argparse
import datetime
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import verify_motor_certificate as verifier


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first",type=Path,required=True)
    parser.add_argument("--recheck",type=Path,required=True)
    parser.add_argument("--period",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError("Choose an unused output")
    paths={key:getattr(args,key) for key in ("first","recheck","period")}
    records={key:json.loads(path.read_text()) for key,path in paths.items()}
    first,recheck,period=(records[key] for key in ("first","recheck","period"))
    gates={}
    for key in ("first","recheck"):
        rec=records[key]
        gates[key+"_unique_crossing"]=rec.get("status")=="unique_fourfold_transverse_crossing_certified"
        gates[key+"_all_recorded_gates"]=bool(rec["gates"]) and all(rec["gates"].values()) and rec["unique_crossing_gate"]
        gates[key+"_shooting_box_inside_period_anchor"]=(len(period["anchor_revalidation"]["shooting_box"])==len(rec["family"]["shooting_box"])==11) and all(verifier.exact_contains(b,x,strict=True)
            for b,x in zip(period["anchor_revalidation"]["shooting_box"],rec["family"]["shooting_box"]))
        gates[key+"_endpoint_parameter_boxes_inside_period_interval"]=set(rec["endpoints"])=={"low","high"} and all(
            verifier.exact_contains(period["parameter_rho"],endpoint["parameter_rho"])
            for endpoint in rec["endpoints"].values())
    gates["period_success"]=period.get("status")=="least_second_crossing_family_period_certified"
    gates["period_all_recorded_gates"]=bool(period["gates"]) and all(period["gates"].values())
    gates["identical_exact_crossing_interval"]=(len(first["exact_endpoint_decimals"])==len(recheck["exact_endpoint_decimals"])==2) and all(Fraction(a)==Fraction(b) for a,b in
        zip(first["exact_endpoint_decimals"],recheck["exact_endpoint_decimals"]))
    gates["recheck_references_first_bytes"]=recheck["center_input_sha256"]==digest(paths["first"])
    gates["period_references_first_bytes"]=period["input_sha256"]==digest(paths["first"])
    shared=set(first["source_sha256"]) & set(recheck["source_sha256"]) & set(period["source_sha256"])
    required={"validate_motor_cycle.py","validate_motor_paired_cycle.py","validate_motor_paired_crossing.py",
              "validate_motor_crossing.py","verify_motor_certificate.py"}
    gates["shared_model_and_proof_sources_identical"]=required.issubset(shared) and all(
        first["source_sha256"][name]==recheck["source_sha256"][name]==period["source_sha256"][name]
        for name in shared)
    checked=0
    for key,path in paths.items():
        valid=True
        for name,expected in records[key]["source_sha256"].items():
            snapshot=path.with_name(path.stem+"."+name)
            valid=valid and snapshot.is_file() and digest(snapshot)==expected
            checked+=1
        gates[key+"_archived_source_hashes"]=valid
    sources={p.name:p.read_bytes() for p in (Path(__file__),Path(verifier.__file__))}
    planned=[args.output.with_name(args.output.stem+"."+name) for name in sources]
    if any(path.exists() for path in planned): raise FileExistsError("A source snapshot already exists")
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status="paired_crossing_proof_links_verified" if all(gates.values()) else "failed",
        inputs={key:dict(path=str(path),sha256=digest(path)) for key,path in paths.items()},
        gates={key:bool(value) for key,value in gates.items()},archived_source_files_checked=checked,
        shared_source_names=sorted(shared),source_sha256={name:hashlib.sha256(data).hexdigest() for name,data in sources.items()},
        statement="Both crossing certificates use the same exact parameter interval and have shooting boxes inside the separately unique least-period anchor. Their endpoint parameter enclosures lie in its uniform parameter interval. Hence both certify the same family and inherit its least period.",
        limitation="This checks exact rational containment and provenance; it does not repeat interval integration.")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for name,data in sources.items(): args.output.with_name(args.output.stem+"."+name).write_bytes(data)
    args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(record["status"],{key:value for key,value in gates.items() if not value})
    if not all(gates.values()): raise SystemExit(2)


if __name__=="__main__": main()
