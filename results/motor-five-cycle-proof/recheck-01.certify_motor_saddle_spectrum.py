"""Certify an archived saddle monodromy after exact dyadic similarity.

An ordinary real-Schur calculation proposes a basis only. Arb verifies its
invertibility and encloses the complete similarity transform. No approximate
eigenvalue or orthogonality assertion is used as evidence.
"""
import argparse
import copy
import datetime
import hashlib
import json
from pathlib import Path
from flint import arb,arb_mat,ctx
import validate_motor_cycle as base
import validate_motor_paired_cycle as paired
from validate_motor_cycle import ZERO,ONE,ball_record,matrix_record,read_ball

def similarity_schur(matrix,basis_path,radius):
    data=json.loads(Path(basis_path).read_text())
    values=[]
    for row in data["basis"]:
        out=[]
        for numerator,denominator in row:
            n,d=int(numerator),int(denominator)
            if d<=0 or d&(d-1): raise ValueError("Basis denominators must be positive powers of two")
            x=arb(n)/arb(d)
            if not x.rad().is_zero(): raise ArithmeticError("Basis must be exact at the working precision")
            out.append(x)
        values.append(out)
    q=arb_mat(values)
    if q.nrows()!=matrix.nrows() or q.ncols()!=matrix.ncols(): raise ValueError("Wrong basis dimension")
    determinant=q.det()
    if determinant.contains(0): raise ArithmeticError("Basis invertibility not certified")
    inverse=q.inv()
    if not all(inverse[i,j].is_finite() for i in range(q.nrows()) for j in range(q.ncols())):
        raise ArithmeticError("Basis inverse not enclosed")
    transformed=inverse*matrix*q
    polynomial=transformed.charpoly()
    n=q.nrows()
    quotient=[ZERO]*n
    quotient[-1]=polynomial[n]
    for k in range(n-2,-1,-1): quotient[k]=polynomial[k+1]+quotient[k+1]
    remainder=polynomial[0]+quotient[0]
    schur=paired.schur_certificate(quotient,radius)
    return dict(passed=bool(remainder.contains(0) and schur["passed"]),
        basis_path=str(basis_path),basis_sha256=hashlib.sha256(Path(basis_path).read_bytes()).hexdigest(),
        exact_dyadic_basis=matrix_record(q),basis_determinant=ball_record(determinant),
        inverse_enclosure=matrix_record(inverse),similar_monodromy=matrix_record(transformed),
        characteristic_coefficients=[ball_record(polynomial[k]) for k in range(n+1)],
        exact_flow_division_remainder=ball_record(remainder),nonflow_schur=schur,
        mathematical_basis="Exact similarity preserves characteristic polynomial. The nonconstant autonomous orbit supplies an exact z-1 factor. All input radii remain in the transformed matrix.")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--basis",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError("Preserve evidence: choose a new output")
    old=json.loads(args.input.read_text())
    if old.get("status") not in ("periodic_six_cluster_orbit_certified","six_cluster_saddle_certified"):
        raise ValueError("Input must certify periodic existence")
    ctx.prec=int(old["arguments"]["bits"])
    if not(old["krawczyk_passed"] and read_ball(old["section_speed"])>0
        and all(x["distinct"] for x in old["pairwise_cluster_distinctions"])):
        raise ValueError("Missing periodic orbit gates")
    # Reuse the validated monodromy enclosure, not its floating-point midpoint.
    matrix=arb_mat([[read_ball(v) for v in row] for row in old["spectrum"]["quotient_full_monodromy"]])
    result=similarity_schur(matrix,args.basis,arb("0.8"))
    record=copy.deepcopy(old)
    record["parent_certificate"]=dict(path=str(args.input),sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        status=old["status"],integration_repeated=False)
    record["postprocessed_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
    record["spectrum"]["similarity_certificate"]=result
    transverse=all(x["test"]["passed"] for x in record["spectrum"]["transverse"].values())
    passed=bool(result["passed"] and transverse)
    record["spectrum"]["raw_characteristic_test_passed"]=old["spectrum"]["passed"]
    record["spectrum"]["passed"]=passed
    record["spectrum"]["floquet_count"]=dict(unstable=1,neutral=1,stable=21) if passed else None
    record["status"]="six_cluster_saddle_certified" if passed else "periodic_six_cluster_orbit_certified"
    record["postprocessor_source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    record["postprocessor_arguments"]={k:str(v) for k,v in vars(args).items()}
    current_snapshot=args.output.with_name(args.output.stem+".executed-postprocessor.py")
    record["postprocessor_source_snapshot"]=current_snapshot.name
    args.output.parent.mkdir(parents=True,exist_ok=True)
    planned=[args.output.with_name(args.output.stem+"."+name) for name in old["source_sha256"]]+[current_snapshot]
    if any(p.exists() for p in planned): raise FileExistsError("A planned source snapshot already exists")
    for name,expected in old["source_sha256"].items():
        snapshot=args.input.with_name(args.input.stem+"."+name)
        data=snapshot.read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected: raise ValueError("Source snapshot hash mismatch: "+name)
        args.output.with_name(args.output.stem+"."+name).write_bytes(data)
    # The inherited integrator may have used an older version of this module.
    # Keep that snapshot under its original source name and store the executing
    # postprocessor separately, so neither provenance chain is overwritten.
    current_snapshot.write_bytes(Path(__file__).read_bytes())
    args.output.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(record["status"],result["nonflow_schur"]["steps"][-1]["margin"]["display"],flush=True)
    if not passed: raise SystemExit(2)

if __name__=="__main__": main()
