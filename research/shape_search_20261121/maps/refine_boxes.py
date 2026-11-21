"""Refine return-box proposals with dense source grids and multiple phase frames."""
import argparse
import itertools
import json
from pathlib import Path
from search import CACHE,START,sample
from box_scout import frames,points_step
from return_scout import attempt
from cover_scout import charts
from validate import write_json
import numpy as np


def scoring(box):
    o,m,inv=charts(box);p=np.array(box["params"]);period=box["case_order"]
    u=np.array(list(itertools.product(np.linspace(-1,1,19),repeat=3)))
    x=o+u@m.T
    gap=1e10
    for step in range(1,period+1):
        x=points_step(x,p)
        z=(x-o)@inv.T
        if step<period:
            gap=min(gap,float(np.min(np.max(abs(z),axis=1)-1)))
    returned=1-float(np.max(np.abs(z)))
    if min(gap,returned)<=.001:
        return None
    shape=z.reshape(19,19,19,3)
    derivative=sum(np.max(abs(np.diff(shape,axis=j)),axis=(0,1,2))*9 for j in range(3))
    return dict(dense_return_margin=returned,dense_avoidance_margin=gap,
                selection_score=float(returned/max(1.,derivative.max())))


def refine(cache,preview,cases):
    lookup={r["case"]:r for r in json.loads(preview.read_text())}
    grid=np.array(list(itertools.product(np.linspace(-1,1,9),repeat=3)))
    for case in cases:
        row=lookup[case];p=np.array(row["params"]);order=row["order"]
        _,tails=sample(p[None,:],START[None,:],100_000,50_000,5,50_000)
        centers,bases,width=frames(tails[0],order);hits=[]
        for phase in range(order):
            for floor in (.005,.03,.08,.18,.35,.60):
                r=attempt(p,order,centers[phase],bases[phase],width[phase],floor,grid)
                if not r:
                    continue
                r.update(case=case,params=p.tolist(),source_phase=phase)
                score=scoring(r)
                if score:
                    r.update(score);hits.append(r)
        hits.sort(key=lambda r:-r["selection_score"])
        for i,hit in enumerate(hits[:6]):
            write_json(cache/"return-refined"/f"{case}-{i}.json",hit)
        print(case,len(hits),[(r["source_phase"],round(r["dense_return_margin"],4)) for r in hits[:3]],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--preview",type=Path,required=True)
    p.add_argument("--cases",type=int,nargs="+",required=True)
    a=p.parse_args()
    refine(a.cache,a.preview,a.cases)
