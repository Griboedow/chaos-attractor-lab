"""Cross-check the first-core polynomials using arc diagrams of the closed braids."""
import argparse
import json
from pathlib import Path

from flint import fmpz_mat

from search import write_json


def arcs(strands,word):
    histories=[[] for _ in range(strands)]
    positions=list(range(strands))
    for k,g in enumerate(word):
        i=abs(g)-1
        left,right=positions[i:i+2]
        over,under=(left,right) if g>0 else (right,left)
        histories[over].append((k,True))
        histories[under].append((k,False))
        positions[i:i+2]=[right,left]
    visit,seen=[],set()
    strand=0
    while strand not in seen:
        seen.add(strand)
        visit.extend(histories[strand])
        strand=positions.index(strand)
    assert strand==0 and len(seen)==strands and len(visit)==2*len(word)
    active=next(k for k,over in reversed(visit) if not over)
    edges=[[None,None,k] for k in range(len(word))]
    for k,over in visit:
        if over:
            edges[k][0]=active
        else:
            edges[k][1]=active
            active=k
    return edges


def reduced_arc_minor(word,edges,t):
    n=len(word)
    matrix=[[0]*n for _ in range(n)]
    for k,(over,inc,out) in enumerate(edges):
        weights=(1-t,t,-1) if word[k]>0 else (t-1,1,-t)
        for j,value in zip((over,inc,out),weights):
            matrix[k][j]+=value
    if n==1:
        return 1
    return int(fmpz_mat([r[:-1] for r in matrix[:-1]]).det())


def remove_units(value,t):
    value=abs(value)
    while value and abs(t)>1 and value % abs(t)==0:
        value//=abs(t)
    return value


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates",type=Path)
    ap.add_argument("polynomials",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    polys={r["case"]:r["alexander_low_to_high"] for r in json.loads(args.polynomials.read_text())["rows"]}
    rows=json.loads(args.certificates.read_text())["rows"]
    evaluations=(-1,2,3,5)
    for row in rows:
        word=row["exact_word"]
        edges=arcs(row["degree"],word)
        coeff=polys[row["case"]]
        for t in evaluations:
            actual=reduced_arc_minor(word,edges,t)
            expect=sum(c*t**k for k,c in enumerate(coeff))
            assert remove_units(actual,t)==remove_units(expect,t),(row["case"],t,actual,expect)
    write_json(args.out,{"cases":len(rows),"values":list(evaluations),
                         "evaluations_checked":len(rows)*len(evaluations),
                         "method":"Independent arc labels on the braid closure; exact integer Fox minors, up to Laurent units."})
    print(len(rows)*len(evaluations),"closed-diagram checks passed")
