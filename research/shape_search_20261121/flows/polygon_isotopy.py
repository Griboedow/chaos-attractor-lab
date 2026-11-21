"""Reduce a frozen dense integer polygon by certified elementary triangle moves."""
import argparse
import gzip
import hashlib
import itertools
import json
import os
from collections import defaultdict
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from exact_knot import exact_diagram,prepare


def sub(a,b):
    return (a[0]-b[0],a[1]-b[1],a[2]-b[2])


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def bounds(vs):
    return tuple(x for coord in zip(*vs) for x in (min(coord),max(coord)))


def buckets(bound,width):
    return itertools.product(*(range(bound[2*i]//width,bound[2*i+1]//width+1) for i in range(3)))


def hit(a,b,c,normal,p,q):
    v=dot(normal,sub(p,a));w=dot(normal,sub(q,a))
    if v*w>0:return False
    if v==w==0:
        # Treat coplanar overlap as an obstacle. A false positive only declines this move.
        return True
    if v==0 and p in (a,c):return False
    if w==0 and q in (a,c):return False
    # Exact projective coordinates of the segment/plane intersection.
    den=w-v
    xyz=tuple(w*x-v*y for x,y in zip(p,q))
    if den<0:den=-den;xyz=tuple(-x for x in xyz)
    for f,g in ((a,b),(b,c),(c,a)):
        relative=tuple(x-den*y for x,y in zip(xyz,f))
        if dot(normal,cross(sub(g,f),relative))<0:return False
    return True


def pass_once(vertices,offset):
    v=vertices[offset:]+vertices[:offset]
    n=len(v);planned=list(range(1,n-1,2))
    # Include every segment before the pass and every proposed chord. Consequently a
    # triangle disjoint from these edges is safe under any ordering of the accepted moves.
    edges=[(i,(i+1)%n) for i in range(n)]+[(i-1,i+1) for i in planned]
    box=[bounds((v[i],v[j])) for i,j in edges]
    span=max(q-p for p,q in zip(map(min,zip(*v)),map(max,zip(*v))))
    width=max(1,span//max(6,min(550,int(n**.58))))
    grid=defaultdict(list)
    for j,eb in enumerate(box):
        for key in buckets(eb,width):grid[key].append(j)
    removed=set()
    for i in planned:
        a,b,c=v[i-1:i+2];normal=cross(sub(b,a),sub(c,a))
        if normal==(0,0,0):continue
        bb=bounds((a,b,c))
        candidates={j for cell in buckets(bb,width) for j in grid[cell]}
        for j in candidates:
            first,last=edges[j]
            if (first,last) in ((i-1,i),(i,i+1),(i-1,i+1)):continue
            eb=box[j]
            if any(bb[t+1]<eb[t] or eb[t+1]<bb[t] for t in (0,2,4)):continue
            if hit(a,b,c,normal,v[first],v[last]):break
        else:
            removed.add(i)
    return [pt for i,pt in enumerate(v) if i not in removed],len(removed)


def reduce(vertices):
    rows=[];no_change=0
    while len(vertices)>5 and no_change<3:
        n=len(vertices)
        vertices,count=pass_once(vertices,1 if len(rows)%2 else 0)
        rows.append([n,count,len(vertices)])
        no_change=no_change+1 if count==0 else 0
    return vertices,rows


def calculate(source,expected,dyadic=False,offset=0):
    with np.load(source) as f:points=f['points']
    if dyadic:
        # All binary64 numbers are dyadic rationals. Multiplication by the largest
        # denominator converts the complete source polygon exactly to integer coordinates.
        ratios=[[float(x).as_integer_ratio() for x in row] for row in points]
        denominator=max(d for row in ratios for _,d in row)
        vertices=[tuple(n*(denominator//d) for n,d in row) for row in ratios]
        encoding=dict(coordinate_mode='exact binary64 dyadics',common_denominator=str(denominator))
    else:
        shift=points.mean(0);scale=np.linalg.norm(np.std(points,axis=0))
        native=np.rint((points-shift)/scale*1e9).astype(np.int64)
        vertices=[tuple(map(int,row)) for row in native.tolist()]
        encoding=dict(coordinate_mode='rounded lattice',integer_sha256=hashlib.sha256(
            native.astype('<i8').tobytes()).hexdigest(),center=shift.tolist(),scale=scale,quantization=10**9)
    if len(set(vertices))!=len(vertices):raise RuntimeError('Coincident input vertices')
    offset=offset%len(vertices)
    joined,rows=reduce(vertices[offset:]+vertices[:offset])
    for i,p in enumerate(joined):
        incoming=sub(p,joined[i-1]);outgoing=sub(joined[(i+1)%len(joined)],p)
        if cross(incoming,outgoing)==(0,0,0) and dot(incoming,outgoing)<=0:
            raise RuntimeError('Overlapping adjacent segments in reduced polygon')
    views=[]
    for rep in (17,35,81,192):
        _,mat,_,_=prepare(points[::128],rep)
        try:d=exact_diagram(joined,mat)
        except ValueError:continue
        if d['coefficients']!=expected:raise RuntimeError('Dense polygon disagrees with orbit witness')
        views.append(dict(camera=mat.tolist(),crossings=d['crossings'],coefficients=d['coefficients']))
    if len(views)<2:raise RuntimeError('Insufficient generic projections')
    # A generic embedding of the final polygon follows from a verified diagram with
    # no spatial edge intersections. The elementary moves then certify the original.
    return dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),**encoding,
        old_vertices=len(vertices),reduction_offset=offset,passes=rows,final_vertices=joined,views=views)


def selfcheck():
    a,b,c=(0,0,0),(10,0,0),(0,10,0)
    normal=cross(sub(b,a),sub(c,a))
    assert hit(a,b,c,normal,(1,1,-1),(1,1,1))
    assert not hit(a,b,c,normal,(8,8,-1),(8,8,1))
    assert not hit(a,b,c,normal,a,(1,2,1))
    assert hit(a,b,c,normal,(5,0,-1),(5,0,1))
    assert not hit(a,b,c,normal,(2,2,1),(1,1,2))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path)
    ap.add_argument('--cache',required=True,type=Path)
    ap.add_argument('--out',type=Path)
    ap.add_argument('--verify',type=Path,help='Rerun every integer triangle test and compare a frozen audit')
    ap.add_argument('--ids',nargs='+',type=int)
    ap.add_argument('--catalog',type=Path,help='Include all accepted stable cycles in the JSONL index')
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--dyadic',action='store_true',help='Use exact binary64 coordinates without quantization')
    ap.add_argument('--offset',type=int,default=0,help='Cyclic vertex offset before applying the first reduction')
    ap.add_argument('--extra',type=Path,help='Period-doubling run folder carrying a refined doubled cycle')
    args=ap.parse_args();selfcheck()
    if args.verify:
        with gzip.open(args.verify,'rt') as f:records=json.load(f)['success']
        for key,row in records.items():
            if args.ids and (not key.isdigit() or int(key) not in args.ids):continue
            got=calculate(args.cache/row['source_relative'],row['views'][0]['coefficients'],
                          row['coordinate_mode']=='exact binary64 dyadics',row.get('reduction_offset',0))
            if got['source_sha256']!=row['source_sha256'] or got['passes']!=row['passes'] or [list(pt) for pt in got['final_vertices']]!=row['final_vertices']:
                raise ValueError('Dense reduction verification failed for '+key)
        print('Verified elementary reductions and final polynomials.');return
    if args.out is None or args.manifest is None:ap.error('Provide --out and --manifest when building')
    args.out.mkdir(exist_ok=args.resume,parents=True)
    manifest=json.loads(args.manifest.read_text());tasks=[]
    for kind in ('stable_cycles','knotted_periodic_skeletons'):
        for r in manifest[kind]:
            if args.ids is None or r['case'] in args.ids:
                tasks.append((str(r['case']),args.cache/r['cloud_relative'],r['alexander_coefficients']))
    if args.catalog:
        seen={t[0] for t in tasks}
        for line in args.catalog.read_text().splitlines():
            r=json.loads(line);name=str(r['case'])
            if name not in seen and (not args.ids or r['case'] in args.ids):
                tasks.append((name,args.cache/r['cloud_relative'],r['alexander_coefficients']))
    if args.extra:
        child=next(r for r in json.loads((args.extra/'report.json').read_text())['trials'] if r.get('cloud_file'))
        tasks.append(('doubled',args.extra/child['cloud_file'],child['coefficients']))
    for name,source,expected in tasks:
        target=args.out/f'case-{name}.json'
        if args.resume and target.exists():continue
        try:report=calculate(source,expected,args.dyadic,args.offset)
        except (ArithmeticError,RuntimeError,ValueError) as exc:
            if not args.catalog:raise
            (args.out/f'failed-{name}.json').write_text(json.dumps(dict(case=name,reason=str(exc)))+'\n')
            print('Unresolved',name,str(exc),flush=True);continue
        report['source_relative']=source.resolve().relative_to(args.cache.resolve()).as_posix()
        target.write_text(json.dumps(report,separators=(',',':'))+'\n')
        print(name,report['old_vertices'],len(report['final_vertices']),[v['crossings'] for v in report['views']],flush=True)


if __name__=='__main__':
    main()
