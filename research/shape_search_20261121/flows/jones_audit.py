"""Compute mirror-normalized Jones polynomials by exact bracket contraction."""
import argparse
import gzip
import json
import os
import random
from collections import defaultdict
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from exact_knot import exact_diagram,prepare


def crossings(events):
    n=len(events);slots=[{} for _ in range(n//2)];signs=[None]*(n//2)
    for k,(idx,over,sign) in enumerate(events):
        slots[idx][over]=((k-1)%n,k);signs[idx]=sign
    res=[]
    for slot,s in zip(slots,signs):
        oi,oo=slot[True];ui,uo=slot[False]
        res.append((oo,uo,oi,ui) if s>0 else (oo,ui,oi,uo))
    return res,signs


def ordering(xs):
    parity=[]
    for edges in xs:
        s=set()
        for e in edges:s.symmetric_difference_update([e])
        parity.append(s)
    rng=random.Random(5779);best=None;best_value=None
    for _ in range(40):
        remaining=set(range(len(xs)));front=set();order=[];sizes=[]
        while remaining:
            choice=min(remaining,key=lambda i:(len(front^parity[i]),rng.random()))
            front^=parity[choice];sizes.append(len(front))
            remaining.remove(choice);order.append(choice)
        score=(max(sizes),sum(2**(x//2) for x in sizes))
        if best is None or score<best_value:best=order;best_value=score
    return best,best_value[0]


def connection(pairing,smoothing,frontier):
    connected={}
    def root(x):
        while connected.get(x,x)!=x:
            x=connected[x]
        return x
    loops=0
    for u,v in (*pairing,*smoothing):
        a,b=root(u),root(v)
        if a==b:loops+=1
        else:connected[a]=b
    groups=defaultdict(list)
    for e in frontier:groups[root(e)].append(e)
    if any(len(v)!=2 for v in groups.values()):raise ValueError('Broken frontier arc')
    return tuple(sorted(tuple(sorted(v)) for v in groups.values())),loops


def polynomial(events,normalize_sign):
    if not events:return {0:1},0,0
    cs,signs=crossings(events);order,width=ordering(cs)
    if width>22:raise ValueError('Frontier too wide')
    frontier=set();states={():{0:1}};peak=0
    for t,i in enumerate(order):
        a,b,c,d=cs[i]
        for e in cs[i]:frontier.symmetric_difference_update([e])
        new={}
        for state,poly in states.items():
            for shift,sm in [(1,((a,b),(c,d))),(-1,((b,c),(d,a)))]:
                pairing,loops=connection(state,sm,frontier)
                if t==len(order)-1:loops-=1
                if loops<0:raise ValueError('Missing final closed component')
                target=new.setdefault(pairing,defaultdict(int))
                factors=[(0,1)]
                for _ in range(loops):factors=[(val+delta,-factor) for val,factor in factors for delta in (-2,2)]
                for k,v in poly.items():
                    for offset,mult in factors:target[k+shift+offset]+=v*mult
        states=new;peak=max(peak,len(states))
        if peak>60000:raise ValueError('State budget')
    w=sum(signs)
    result={int(k+normalize_sign*3*w):int(v*(-1 if w%2 else 1)) for k,v in states[()].items() if v}
    if any(k%4 for k in result):raise ValueError('Non-integral knot exponents')
    return {k//4:v for k,v in result.items()},width,peak


def canon(poly):
    direct=sorted(poly.items())
    mirror=sorted((-k,v) for k,v in poly.items())
    return [list(row) for row in min(direct,mirror)]


def controls():
    t=np.linspace(0,2*np.pi,280,endpoint=False)+.0031
    def torus(p,q):
        return np.column_stack(((2+np.cos(q*t))*np.cos(p*t),(2+np.cos(q*t))*np.sin(p*t),np.sin(q*t)))
    eight=np.column_stack(((2+np.cos(2*t))*np.cos(3*t),(2+np.cos(2*t))*np.sin(3*t),np.sin(4*t)))
    tests=[(torus(2,3),canon({1:1,3:1,4:-1})),(eight,canon({-2:1,-1:-1,0:1,1:-1,2:1})),
           (torus(2,5),canon({2:1,4:1,5:-1,6:1,7:-1}))]
    for norm in (-1,1):
        passed=True
        for points,target in tests:
            for seed in (17,35,81):
                v,m,_,_=prepare(points,seed)
                try:got=canon(polynomial(exact_diagram(v,m,gauss=True),norm)[0])
                except ValueError:passed=False;break
                if got!=target:passed=False;break
            if not passed:break
        if passed:return norm
    raise ValueError('Independent known Jones polynomials failed')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--catalog',required=True,type=Path)
    args=ap.parse_args();args.out.mkdir(exist_ok=False,parents=True)
    norm=controls();print('Knot controls passed',norm,flush=True)
    cases={str(json.loads(t)['case']) for t in args.catalog.read_text().splitlines()}
    distinctions=set();groups=set();handed=set();resolved=[];unresolved=[]
    for file in sorted(args.source.glob('case-*.json')):
        name=file.stem[5:]
        if name not in cases:continue
        data=json.loads(file.read_text());vertices=data['final_vertices']
        projections=sorted(data['views'],key=lambda v:v['crossings'])[:2];outputs=[]
        try:
            for view in projections:
                events=exact_diagram(vertices,view['camera'],gauss=True)
                p,width,peak=polynomial(events,norm)
                expected=abs(sum((-1)**(k%2)*v for k,v in p.items()))
                alex=view['coefficients']
                if expected!=abs(sum((-1)**k*v for k,v in enumerate(alex))):raise ValueError('Determinant disagreement')
                x,y,z=view['camera']
                orientation=x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0])
                physical=sorted([[k if orientation>0 else -k,v] for k,v in p.items()])
                outputs.append(dict(polynomial=canon(p),handed_polynomial=physical,width=width,states=peak,
                                   camera=view['camera'],crossings=len(events)//2))
            if outputs[0]['polynomial']!=outputs[1]['polynomial']:raise ValueError('View disagreement')
            if outputs[0]['handed_polynomial']!=outputs[1]['handed_polynomial']:raise ValueError('Camera handedness disagreement')
            result=dict(case=int(name),views=outputs,alexander=alex)
            resolved.append(result);distinctions.add((tuple(alex),str(outputs[0]['polynomial'])))
            groups.add(str(outputs[0]['polynomial']))
            handed.add(str(outputs[0]['handed_polynomial']))
        except ValueError as exc:unresolved.append(dict(case=int(name),reason=str(exc)))
        payload=dict(resolved=resolved,unresolved=unresolved,normalization_sign=norm,
                     different_mirror_quotiented_Jones=len(groups),different_handed_Jones=len(handed),
                     different_Alexander_Jones_pairs=len(distinctions))
        (args.out/'summary.json.gz').write_bytes(gzip.compress(json.dumps(payload,separators=(',',':')).encode(),mtime=0))
        if len(resolved)%25==0:print(len(resolved),len(groups),len(distinctions),flush=True)


if __name__=='__main__':
    main()
