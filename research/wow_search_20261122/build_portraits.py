"""Bin native planar points or an orthographic view of a 3D cloud for the portrait atlas."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def bin_cloud(source, output, meta, key, view, virtual):
    with np.load(source) as data:
        points=np.asarray(data[key],dtype=float)
    if points.ndim!=2 or points.shape[1] not in (2,3) or not np.all(np.isfinite(points)):
        raise ValueError('Expected finite native plane or space coordinates')
    if view is None:
        if points.shape[1]!=2:
            raise ValueError('Specify an orthographic camera for a spatial cloud')
        projection=np.eye(2)
    else:
        e,a=np.deg2rad(view)
        projection=np.array([[-np.sin(a),np.cos(a),0],[-np.cos(a)*np.sin(e),-np.sin(a)*np.sin(e),np.cos(e)]])
    center=(points.max(0)+points.min(0))/2
    xy=(points-center)@projection.T
    half=abs(xy).max()*1.04
    h=np.histogram2d(xy[:,0],xy[:,1],bins=1000,range=[[-half,half]]*2)[0].T[::-1].astype(np.uint32)
    assert h.sum()==len(points)
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,counts=h)
    row=dict(source=virtual,key=key,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
             center=center.tolist(),basis=np.eye(points.shape[1]).tolist(),projection=projection.tolist(),
             points=len(points),grid_size=1000,half_width=float(half))
    meta.parent.mkdir(parents=True,exist_ok=True)
    with meta.open('w',encoding='utf-8',newline='\n') as f:
        json.dump(row,f,indent=2);f.write('\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('output',type=Path);p.add_argument('meta',type=Path)
    p.add_argument('--key',default='points')
    p.add_argument('--view',type=float,nargs=2,metavar=('ELEVATION','AZIMUTH'))
    p.add_argument('--virtual',required=True)
    a=p.parse_args()
    bin_cloud(a.source,a.output,a.meta,a.key,a.view,a.virtual)
