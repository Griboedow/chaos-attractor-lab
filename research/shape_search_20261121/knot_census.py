"""Compare core and cycle invariants with the pinned KnotInfo prime-knot census."""
import argparse
import ast
import csv
import hashlib
import gzip
import json
from collections import Counter,defaultdict
from pathlib import Path

from index_campaign import read,write,RESULTS

COMMIT="219c2d982b89fb92b6ab0c556f28cdd766a18734"
SHA256="ed4447ff81dc50e2561764760d631356186c34e2dd35f89cc84e0a4d8fedd3d0"


def node(n):
    if isinstance(n,ast.Constant) and type(n.value) is int:
        return {0:n.value}
    if isinstance(n,ast.Name) and n.id=="t":
        return {1:1}
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):
        sign=-1 if isinstance(n.op,ast.USub) else 1
        return {k:sign*v for k,v in node(n.operand).items()}
    assert isinstance(n,ast.BinOp),ast.dump(n)
    a,b=node(n.left),node(n.right)
    if isinstance(n.op,(ast.Add,ast.Sub)):
        sign=-1 if isinstance(n.op,ast.Sub) else 1
        return {k:a.get(k,0)+sign*b.get(k,0) for k in a.keys()|b.keys()}
    if isinstance(n.op,ast.Pow):
        assert len(a)==1 and 0 in b and len(b)==1
        (power,coeff),=a.items()
        assert coeff==1
        return {power*b[0]:1}
    if isinstance(n.op,ast.Div):
        assert len(b)==1
        (den,coeff),=b.items()
        assert coeff in (-1,1)
        return {k-den:v*coeff for k,v in a.items()}
    assert isinstance(n.op,ast.Mult)
    out=Counter()
    for k,v in a.items():
        for l,w in b.items():
            out[k+l]+=v*w
    return dict(out)


def get(row,stem):
    vector=row[stem+"_polynomial_vector"]
    if vector:
        if vector=="1":
            return {0:1}
        lo,hi,*values=ast.literal_eval(vector)
        assert len(values)==hi-lo+1
        return {lo+i:v for i,v in enumerate(values) if v}
    expr=row[stem+"_polynomial"].strip().replace("^","**")
    return {power:v for power,v in node(ast.parse(expr,mode="eval").body).items() if v}


def composite_signatures(primes,fields,limit=13):
    joint=set()
    alex=set()
    factor_lists=0
    def extend(start,total,n,coeff,values):
        nonlocal factor_lists
        for i in range(start,len(primes)):
            crossing,poly,pair=primes[i]
            if crossing+total>limit:
                break
            product=[0]*(len(coeff)+len(poly)-1)
            for k,v in enumerate(coeff):
                for l,w in enumerate(poly):
                    product[k+l]+=v*w
            product=tuple(product)
            options=set()
            # Mirrors can be chosen independently before the connected sum.
            for old,new in values:
                for left,right in (pair,pair[::-1]):
                    f=tuple(v*w%m for v,w,m in zip(old,left,fields))
                    g=tuple(v*w%m for v,w,m in zip(new,right,fields))
                    options.add(tuple(sorted([f,g])))
            if n:
                factor_lists+=1
                alex.add(product)
                joint.update((product,s) for s in options)
            extend(i,total+crossing,n+1,product,options)
    one=(1,)*len(fields)
    extend(0,0,0,(1,),{(one,one)})
    return joint,alex,factor_lists


def mirror_signature(j,evals):
    values=[]
    for direction in (-1,1):
        xyz=[]
        for e in evals:
            p=e["prime"]
            t=pow(e["bracket_variable_A"],4*direction,p)
            xyz.append(sum(v*pow(t,k,p) for k,v in j.items())%p)
        values.append(tuple(xyz))
    return tuple(sorted(values))


def compare(path):
    with path.open("rb") as f:
        assert hashlib.file_digest(f,"sha256").hexdigest()==SHA256
    jones=read(RESULTS/"braids/jones.json")
    evals=jones["evaluations"]
    index=defaultdict(list)
    coarse=defaultdict(list)
    counts=Counter()
    primes=[]
    with path.open(encoding="utf-8",newline="") as f:
        for item in csv.DictReader(f,delimiter="|"):
            if not item["crossing_number"].isdigit():
                assert item["name"]=="Name"
                continue
            crossing=int(item["crossing_number"])
            counts[crossing]+=1
            a,j=get(item,"alexander"),get(item,"jones")
            lo,hi=min(a),max(a)
            sign=1 if a[lo]>0 else -1
            poly=tuple(sign*a.get(i,0) for i in range(lo,hi+1))
            # The snapshot stores zero in the unknot determinant field.
            # Both published polynomials equal one at that row.
            det=1 if crossing==0 else int(item["determinant"])
            assert abs(sum(v*(-1 if k%2 else 1) for k,v in enumerate(poly)))==det
            assert poly==poly[::-1] and abs(sum(poly))==1
            assert sum(j.values())==1 and abs(sum(v*(-1 if k%2 else 1) for k,v in j.items()))==det
            signature=mirror_signature(j,evals)
            if 0<crossing<=10:
                primes.append((crossing,poly,signature))
            coarse[poly].append(item["name"])
            index[(poly,signature)].append(item["name"])
    assert counts[13]==9988 and counts[12]==2176 and counts[0]==1
    composite,composite_poly,n_lists=composite_signatures(sorted(primes),[q["prime"] for q in evals])
    pmap={d["case"]:tuple(d["alexander_low_to_high"]) for d in read(RESULTS/"braids/core-knots.json")["rows"]}
    braid=[]
    for row in jones["rows"]:
        case=row["case"]
        poly=pmap[case]
        sig=tuple(tuple(v) for v in row["jones_pair_mod_primes"])
        hits=index.get((poly,sig),[])
        braid.append({"case":case,"alexander_matches":len(coarse.get(poly,[])),
                      "alexander_and_jones_matches":hits,"compatible_low_crossing_composite":(poly,sig) in composite})
    assert "8_19" in next(r["alexander_and_jones_matches"] for r in braid if r["case"]==91625)
    flow=[]
    fj=json.loads(gzip.decompress((RESULTS/"flows/jones-polynomials.json.gz").read_bytes()))
    flow_jones={r["case"]:r for r in fj["resolved"]}
    with (RESULTS/"flows/periodic-catalog.csv").open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            poly=tuple(json.loads(r["alexander_coefficients"]))
            case=int(r["case"])
            proof=flow_jones[case]
            assert tuple(proof["alexander"])==poly
            sig=mirror_signature(dict(proof["views"][0]["handed_polynomial"]),evals)
            hits=index.get((poly,sig),[])
            flow.append({"case":case,"alexander_matches":coarse.get(poly,[]),"alexander_and_jones_matches":hits,
                         "compatible_low_crossing_composite":(poly,sig) in composite})
    summary={"repository":"https://github.com/soehms/database_knotinfo","snapshot_commit":COMMIT,
             "csv_sha256":SHA256,"rows_by_crossing_number":dict(sorted(counts.items())),
             "prime_rows":sum(counts.values())-1,"unknot_rows":counts[0],
             "braid_cases":len(braid),"braids_without_any_census_match":sum(not r["alexander_and_jones_matches"] for r in braid),
             "cycle_cases":len(flow),"cycles_without_any_census_match":sum(not r["alexander_matches"] for r in flow),
             "cycles_without_any_joint_census_match":sum(not r["alexander_and_jones_matches"] for r in flow),
             "connected_sum_factor_multisets_checked":n_lists,
             "braids_without_prime_or_tested_composite_match":sum(not r["alexander_and_jones_matches"] and not r["compatible_low_crossing_composite"] for r in braid),
             "cycles_without_prime_or_tested_composite_match":sum(not r["alexander_and_jones_matches"] and not r["compatible_low_crossing_composite"] for r in flow),
             "scope":"Reference matches are necessary invariant tests. The catalog contains prime knots through 13 crossings and the unknot. The additional composite test covers connected-sum diagrams with factor crossing counts summing to at most 13 and arbitrary mirrors. It assumes no crossing-number additivity theorem."}
    write(RESULTS/"knot-census.json",summary)
    write(RESULTS/"knot-census-matches.json.gz",{"braids":braid,"cycles":flow})
    print(summary,flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv",type=Path)
    compare(ap.parse_args().csv)
