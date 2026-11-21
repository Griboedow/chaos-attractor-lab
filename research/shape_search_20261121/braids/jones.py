"""Distinguish the certified braid closures by normalized Kauffman-bracket evaluations."""
import argparse
import json
from functools import lru_cache
from pathlib import Path


def identity(n):
    return tuple(range(n,2*n))+tuple(range(n))


@lru_cache(maxsize=None)
def cap(diagram,i):
    n=len(diagram)//2
    a,b=i+n,i+n+1
    c,d=diagram[a],diagram[b]
    if c==b:
        return diagram,1
    out=list(diagram)
    out[a],out[b]=b,a
    out[c],out[d]=d,c
    return tuple(out),0


@lru_cache(maxsize=None)
def closed_count(diagram):
    n=len(diagram)//2
    seen=set()
    loops=0
    for i in range(2*n):
        if i in seen:
            continue
        loops+=1
        j=i
        while j not in seen:
            seen.add(j);seen.add(diagram[j])
            j=(diagram[j]+n)%(2*n)
    return loops


def bracket(n,word,a,p):
    inv=pow(a,-1,p)
    delta=(-a*a-inv*inv)%p
    states={identity(n):1}
    for gen in word:
        f,g=(a,inv) if gen>0 else (inv,a)
        result={}
        for d,coeff in states.items():
            result[d]=(result.get(d,0)+f*coeff)%p
            e,loop=cap(d,abs(gen)-1)
            result[e]=(result.get(e,0)+g*coeff*(delta if loop else 1))%p
        states={d:v for d,v in result.items() if v}
    writhe=sum(1 if gen>0 else -1 for gen in word)
    return sum(v*pow(delta,closed_count(d)-1,p) for d,v in states.items())*pow(-a**3,-writhe,p)%p


def torus_signature(u,v,evals):
    pair=[]
    for direction in (-1,1):
        readings=[]
        for prime,a in evals:
            t=pow(a,4*direction,prime)
            val=pow(t,(u-1)*(v-1)//2,prime)
            val*=1-pow(t,u+1,prime)-pow(t,v+1,prime)+pow(t,u+v,prime)
            readings.append(val*pow(1-t*t,-1,prime)%prime)
        pair.append(tuple(readings))
    return tuple(sorted(pair))


def run(certs,polys,out):
    braids=json.loads(certs.read_text())["rows"]
    cmap={r["case"]:r for r in json.loads(polys.read_text())["rows"]}
    evals=[(998244353,2),(1000000007,3),(1000000009,11)]
    eighth=pow(3,(998244353-1)//8,998244353)
    assert pow(eighth,4,998244353)==998244352
    records=[]
    combined=set()
    by_poly={}
    excluding_torus=set()
    for row in braids:
        n,word=row["degree"],row["exact_word"]
        poly=tuple(cmap[row["case"]]["alexander_low_to_high"])
        d=row["verified_core_determinant"]
        value=bracket(n,word,eighth,998244353)
        assert value in (d%998244353,(-d)%998244353),(row["case"],value,d)
        direct=[bracket(n,word,a,p) for p,a in evals]
        mirror=[bracket(n,word,pow(a,-1,p),p) for p,a in evals]
        # Moving the cut in a closed braid and either stabilization preserve it.
        p,a=evals[0]
        assert direct[0]==bracket(n,word[1:]+word[:1],a,p)
        assert direct[0]==bracket(n+1,word+[n],a,p)==bracket(n+1,word+[-n],a,p)
        signature=tuple(sorted([tuple(direct),tuple(mirror)]))
        matches=[pair for pair in cmap[row["case"]]["torus_polynomial_matches"]
                 if torus_signature(*pair,evals)==signature]
        if not matches:
            excluding_torus.add((poly,signature))
        combined.add((poly,signature))
        by_poly.setdefault(poly,set()).add(signature)
        records.append({"case":row["case"],"jones_pair_mod_primes":signature,"determinant_check":d,
                        "torus_parameters_compatible_with_both_invariants":matches})
    out.write_text(json.dumps({"evaluations":[{"prime":p,"bracket_variable_A":a} for p,a in evals],
                   "mirror_convention":"Unordered pair of V at A and its inverse",
                   "different_core_types_lower_bound":len(combined),
                   "different_non_torus_core_types_lower_bound":len(excluding_torus),
                   "alexander_classes_split_by_jones":sum(len(x)>1 for x in by_poly.values()),
                   "independent_per_case_checks":["A^4=-1 determinant","cyclic cut","positive and negative Markov stabilization"],
                   "rows":records},indent=2)+"\n",encoding="utf-8")
    print(len(records),"closures,",len(combined),"distinct mirror-insensitive signatures",flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates",type=Path)
    ap.add_argument("polynomials",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    run(args.certificates,args.polynomials,args.out)
