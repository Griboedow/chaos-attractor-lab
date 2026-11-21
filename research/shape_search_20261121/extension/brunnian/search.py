"""Search balanced pure three-strand braids and verify their integral triple invariants."""
import argparse
from collections import Counter
from math import gcd, isqrt
import json
from pathlib import Path

# Truncate the noncommutative Magnus expansion after degree two, setting X3=0.
def multiply(left, right):
    x,y,a,b,c,d=left
    u,v,e,f,g,h=right
    return x+u,y+v,a+e+x*u,b+f+x*v,c+g+y*u,d+h+y*v


def inverse(word):
    x,y,a,b,c,d=word
    return -x,-y,x*x-a,x*y-b,x*y-c,y*y-d


ZERO=(0,0,0,0,0,0)
MERIDIANS=[(1,0,0,0,0,0),(0,1,0,0,0,0),ZERO]


def magnus(word):
    strands=[(i,ZERO) for i in range(3)]
    for letter in word:
        slot=abs(letter)-1
        for _ in range(2):
            (i,p),(j,q)=strands[slot:slot+2]
            if letter>0:
                mer=multiply(multiply(inverse(p),MERIDIANS[i]),p)
                strands[slot:slot+2]=[(j,multiply(q,inverse(mer))),(i,p)]
            else:
                mer=multiply(multiply(inverse(q),MERIDIANS[j]),q)
                strands[slot:slot+2]=[(j,q),(i,multiply(p,mer))]
    assert [i for i,_ in strands]==[0,1,2]
    return strands[2][1]


def rank(word):
    # Prefer short diagrams; at fixed length distribute the exchanges around the ring.
    changes=sum(abs(a)!=abs(b) for a,b in zip(word,word[1:]+word[:1]))
    return len(word),-changes,word


def boundary(area):
    # A near-square integer rectangle with one extra partial row.
    width=isqrt(area)
    height,extra=divmod(area,width)
    steps=[1]*width+[2]*height
    if extra:
        steps+= [-1]*(width-extra)+[2]+[-1]*extra+[-2]*(height+1)
    else:
        steps+= [-1]*width+[-2]*height
    return tuple(steps)


def prime_to_thirty(value):
    for prime in (2,3,5):
        while value%prime==0:
            value//=prime
    return value


def search(max_blocks, max_area):
    counts=Counter()
    best={}
    def visit(word,x,y,area):
        n=len(word)
        if n and not x and not y and word[0]!=-word[-1]:
            value=abs(area)
            counts[n,value]+=1
            if value and (value not in best or rank(word)<rank(best[value])):
                best[value]=word
        if n+abs(x)+abs(y)>max_blocks or n==max_blocks:
            return
        for letter,dx,dy in ((1,1,0),(-1,-1,0),(2,0,1),(-2,0,-1)):
            if word[-1]==-letter:
                continue
            visit(word+(letter,),x+dx,y+dy,area+x*dy)
    # Fix a positive first A block. The reported counts use this convention.
    visit((1,),1,0,0)
    enumerated=set(best)
    for area in range(1,max_area+1):
        if gcd(area,30)==1 and area not in best:
            best[area]=boundary(area)
    rows=[]
    for value,word in sorted(best.items()):
        w=magnus(word)
        assert w[:3]==(0,0,0) and w[5]==0 and w[3]==-w[4] and abs(w[3])==value,(word,w,value)
        rows.append(dict(case=30300000+value,word=list(word),mu_123=w[3],mu_213=w[4],
                         origin="word enumeration" if value in enumerated else "integer boundary",
                         prime_to_30=prime_to_thirty(value),width=0.27,radius=3.0,degrees=[2,3,5],fiber_contraction=.105,
                         fiber_offset=.55,base_ripple=.12,base_shift=.137))
    return dict(max_pure_blocks=max_blocks,alphabet={"1":"sigma_1 squared","2":"sigma_2 squared"},
                first_block=1,maximum_added_area=max_area,number_of_distinct_prime_to_30_parts=len({r["prime_to_30"] for r in rows}),
                cyclically_reduced_balanced_words=sum(counts.values()),
                counts_by_length_and_abs_mu={f"{n}/{v}":z for (n,v),z in sorted(counts.items())},
                cases=rows)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-blocks",type=int,default=16)
    ap.add_argument("--max-area",type=int,default=256)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    result=search(args.max_blocks,args.max_area)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(result["cyclically_reduced_balanced_words"],[(r["case"],r["mu_123"]) for r in result["cases"]],flush=True)
