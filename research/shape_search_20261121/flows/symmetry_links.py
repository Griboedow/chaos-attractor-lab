"""Compute exact pairwise linking numbers of symmetry-related polygonal cycles."""
import argparse
import json
import os
from collections import defaultdict
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.spatial import cKDTree
from cycle_scan import resample


def link_integer(a,b,matrix):
    def transform(v):
        return [[sum(int(p[k])*int(matrix[k][j]) for k in range(3)) for j in range(3)] for p in v]
    def edges(v):
        out=[]
        for i,p in enumerate(v):
            q=v[(i+1)%len(v)]
            out.append((p,tuple(y-x for x,y in zip(p,q)),min(p[0],q[0]),max(p[0],q[0]),min(p[1],q[1]),max(p[1],q[1])))
        return out
    first=edges(transform(a));second=edges(transform(b))
    # Integer planar bounding boxes preserve all possible intersections at dense meshes.
    # An edge can occupy several cells. A set removes repeated candidate pairs.
    span=max(max(e[j+1] for e in second)-min(e[j] for e in second) for j in (2,4))
    cell=max(1,span//max(5,min(600,int(len(second)**.56))))
    def keys(edge):
        x0,x1,y0,y1=edge[-4:]
        return ((x,y) for x in range(x0//cell,x1//cell+1) for y in range(y0//cell,y1//cell+1))
    grid=defaultdict(list)
    for i,edge in enumerate(second):
        for pos in keys(edge):grid[pos].append(i)
    signed=0;crossings=0
    for edge in first:
        p,u,x0,x1,y0,y1=edge
        indices={j for pos in keys(edge) for j in grid[pos]}
        for i in indices:
            q,v,s0,s1,t0,t1=second[i]
            if x1<s0 or s1<x0 or y1<t0 or t1<y0:continue
            dx=q[0]-p[0];dy=q[1]-p[1]
            den=u[0]*v[1]-u[1]*v[0]
            sn=dx*v[1]-dy*v[0];tn=dx*u[1]-dy*u[0]
            if den==0:
                if sn==tn==0:raise ValueError('Collinear projected edges')
                continue
            sign=1 if den>0 else -1
            if den<0:den=-den;sn=-sn;tn=-tn
            if sn<0 or sn>den or tn<0 or tn>den:continue
            if sn in (0,den) or tn in (0,den):raise ValueError('Endpoint crossing')
            depth=(p[2]-q[2])*den+sn*u[2]-tn*v[2]
            if depth==0:raise ValueError('Intersecting components')
            crossings+=1;signed+=sign if depth>0 else -sign
    if signed%2:raise ValueError('Odd intercomponent crossing sum')
    return signed//2,crossings


def check_witness(record):
    for view in record['views']:
        for i,first in enumerate(record['integer_polygons']):
            for j,second in enumerate(record['integer_polygons'][:i]):
                lk,cross=link_integer(first,second,view['camera'])
                assert (lk,cross)==(view['linking'][i][j],view['crossings'][i][j])


def check_controls():
    t=np.arange(120)*2*np.pi/120+.00172
    first=np.rint(1e8*np.array([np.cos(t),np.sin(t),0*t]).T).astype(np.int64)
    matrix=[[3811,1113,1822],[-3325,7591,2299],[8135,4192,5341]]
    for shift,expected in [(1,1),(3,0)]:
        other=np.rint(1e8*np.array([shift+.5*np.cos(t),0*t,.5*np.sin(t)]).T).astype(np.int64)
        assert abs(link_integer(first,other,matrix)[0])==expected


def matrices(case, p):
    if p.shape==(3,9) and not np.any(p!=p[0]):
        cyclic=np.eye(3)[:,[1,2,0]]
        positive=[np.linalg.matrix_power(cyclic,k) for k in range(3)]
        return positive+([-m for m in positive] if np.all(p[:,4:7]==0) else [])
    if p.shape==(11,):
        return [np.array([[np.cos(a),np.sin(a),0],[-np.sin(a),np.cos(a),0],[0,0,1.]]) for a in 2*np.pi*np.arange(int(p[0]))/p[0]]
    return [np.eye(3)]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--in',dest='sources',nargs='+',type=Path)
    ap.add_argument('--out',type=Path)
    ap.add_argument('--verify',type=Path)
    ap.add_argument('--vertices',type=int,default=640)
    ap.add_argument('--ids',type=int,nargs='+')
    args=ap.parse_args()
    check_controls()
    if args.verify:
        for f in args.verify.glob('case-*.json'):
            check_witness(json.loads(f.read_text()))
        print('Exact linking witnesses and Hopf/unlink controls passed.')
        return
    args.out.mkdir(parents=True,exist_ok=False)
    for source in args.sources:
        for path in sorted(source.glob('case-*.npz')):
            case=int(path.stem[5:])
            if args.ids and case not in args.ids:continue
            record=json.loads(path.with_suffix('.json').read_text());p=np.array(record['parameters'])
            z=np.load(path)['points'];scale=np.linalg.norm(np.std(z,axis=0))
            transforms=matrices(case,p);clouds=[];trees=[];using=[]
            for r in transforms:
                v=z@r
                if any(np.quantile(tree.query(v[::32],workers=1)[0],.95)/scale<.006 for tree in trees):continue
                clouds.append(v);trees.append(cKDTree(v[::4]));using.append(r)
            if len(clouds)<2:continue
            polygons=[np.rint(resample(v,args.vertices)/scale*1e9).astype(np.int64) for v in clouds]
            views=[];consistent=True
            for seed in [313,718,2179]:
                camera=np.rint(np.linalg.qr(np.random.default_rng(seed).normal(size=(3,3)))[0]*10000).astype(int)
                lk=np.zeros((len(clouds),len(clouds)),int)
                cross=np.zeros_like(lk)
                try:
                    for i in range(len(clouds)):
                        for j in range(i):
                            val,count=link_integer(polygons[i],polygons[j],camera)
                            lk[i,j]=lk[j,i]=val;cross[i,j]=cross[j,i]=count
                except ValueError:continue
                if views and np.any(lk!=np.array(views[0]['linking'])):consistent=False
                views.append(dict(camera=camera.tolist(),linking=lk.tolist(),crossings=cross.tolist()))
            if consistent and len(views)>=2 and np.max(np.abs(views[0]['linking']))>0:
                out=dict(case=case,family=record['family'],period=record['period'],component_count=len(clouds),
                         transforms=[r.tolist() for r in using],common_rms=scale,views=views,
                         integer_polygons=[p.tolist() for p in polygons])
                (args.out/f'case-{case}.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
                if args.vertices==640:
                    np.savez_compressed(args.out/f'case-{case}.npz',points=np.concatenate(clouds),
                                        component=np.repeat(np.arange(len(clouds)),len(z)),
                                        offsets=np.arange(len(clouds)+1)*len(z))
                print(case,len(clouds),views[0]['linking'],flush=True)


if __name__=='__main__':
    main()
