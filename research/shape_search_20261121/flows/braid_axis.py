"""Recognize a closed two-strand braid of a fixed integer polygon, when possible."""
import argparse
import gzip
import json
import os
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.optimize import linprog


def certify(vertices,projection,center):
    u,v=projection
    points=[(sum(x*y for x,y in zip(p,u))-center[0],sum(x*y for x,y in zip(p,v))-center[1]) for p in vertices]
    ends=zip(points,points[1:]+points[:1]);areas=[];w=0
    for a,b in ends:
        value=a[0]*b[1]-a[1]*b[0]
        if value<=0:return None
        areas.append(value)
        if a[1]<=0<b[1]:w+=1
    return dict(winding=w,minimum_oriented_area=str(min(areas)))


def search(vertices,count,seed):
    v=np.asarray(vertices,dtype=float);scale=max(np.ptp(v,axis=0));w=v/scale
    rng=np.random.default_rng(seed)
    dirs=[(1,0,0),(0,1,0),(0,0,1),*rng.normal(size=(count,3))]
    for axis in dirs:
        axis=np.asarray(axis,dtype=float);axis/=np.linalg.norm(axis)
        u=np.cross(axis,np.eye(3)[np.argmin(abs(axis))]);u/=np.linalg.norm(u)
        orient=np.cross(axis,u)
        for sign in (1,-1):
            matrix=np.rint(100000*np.array([u,orient*sign])).astype(int)
            p=w@matrix.T/100000
            q=np.roll(p,-1,axis=0);edge=q-p;dist=np.linalg.norm(edge,axis=1)
            if min(dist)==0:continue
            a=np.column_stack((edge[:,1]/dist,-edge[:,0]/dist,np.ones(len(p))))
            b=(p[:,0]*q[:,1]-p[:,1]*q[:,0])/dist
            out=linprog([0,0,-1],A_ub=a,b_ub=b,bounds=[(-3,3),(-3,3),(.00000001,None)],method='highs')
            if not out.success:continue
            center=[int(round(x*scale*100000)) for x in out.x[:2]]
            projection=matrix.tolist()
            record=certify(vertices,projection,center)
            if record and record['winding']==2:
                return dict(projection=projection,center=center,**record)
    return None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path)
    ap.add_argument('--out',type=Path)
    ap.add_argument('--verify',type=Path,help='Verify the merged axis certificates in the packet')
    ap.add_argument('--reductions',type=Path,help='Compressed dense-polygon reductions')
    ap.add_argument('--directions',type=int,default=700)
    ap.add_argument('--names',nargs='+')
    a=ap.parse_args()
    if a.verify:
        with gzip.open(a.reductions,'rt') as stream:reduced=json.load(stream)['success']
        for name,row in json.loads(a.verify.read_text()).items():
            cert=row['certificate']
            if cert is None:continue
            vertices=reduced[name]['final_vertices']
            got=certify(vertices,cert['projection'],cert['center'])
            if got is None or got['winding']!=2 or got['minimum_oriented_area']!=cert['minimum_oriented_area']:
                raise ValueError('Angular witness failed for '+name)
            coeff=reduced[name]['views'][0]['coefficients'];det=abs(sum((-1)**i*x for i,x in enumerate(coeff)))
            if cert['torus_parameters_up_to_mirror']!=[2,det]:raise ValueError('Torus classification disagrees')
        print('All stored angular and winding witnesses verified.');return
    if a.source is None or a.out is None:ap.error('Provide --source and --out when searching')
    a.out.mkdir(exist_ok=False,parents=True)
    inputs=sorted(a.source.glob('case-*.json'))
    for f in inputs:
        name=f.stem[5:]
        if a.names and name not in a.names:continue
        data=json.loads(f.read_text());v=data['final_vertices']
        result=search(v,a.directions,774815)
        if result:
            coeff=data['views'][0]['coefficients']
            det=abs(sum(x*(-1)**k for k,x in enumerate(coeff)))
            result.update(torus_parameters_up_to_mirror=[2,det],coefficients=coeff)
        out=dict(name=name,directions=a.directions,certificate=result)
        (a.out/f.name).write_text(json.dumps(out,indent=2)+'\n')
        print(name,(result or {}).get('torus_parameters_up_to_mirror'),flush=True)


if __name__=='__main__':
    main()
