"""Curvature-adaptive numerical forward-closure proposals, never interval proofs.

Vertex/center sampling is only a proposal heuristic. All failures are retained.
A box intersecting the B projection axis is rejected exactly at the geometric
formula level (binary64 evaluation here still needs interval confirmation).
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit
from trace_raar_cycle_seeds import value
from search_raar_circles import step

@njit(cache=True,inline="always")
def marked(mask,key):
    return (mask[key//8] & (1 << (key%8))) != 0

@njit(cache=True,inline="always")
def mark(mask,key):
    mask[key//8] |= np.uint8(1 << (key%8))

@njit(cache=True)
def occupied_keys(mask):
    count=0
    for byte in mask:
        v=int(byte)
        while v:
            count+=1;v &= v-1
    result=np.empty(count,np.int64);at=0
    for i in range(len(mask)):
        if not mask[i]:continue
        for bit in range(8):
            if mask[i] & (1 << bit):result[at]=8*i+bit;at+=1
    return result

@njit(cache=True)
def initial_mask(keys,shape,pad):
    mask=np.zeros((shape[0]*shape[1]*shape[2]+7)//8,np.uint8);ny,nz=shape[1],shape[2]
    for key in keys:
        ix=key//(ny*nz);iy=(key//nz)%ny;iz=key%nz
        for a in range(max(0,ix-pad),min(shape[0],ix+pad+1)):
            for b in range(max(0,iy-pad),min(ny,iy+pad+1)):
                for c in range(max(0,iz-pad),min(nz,iz+pad+1)):mark(mask,(a*ny+b)*nz+c)
    return mask

from raar_adaptive_grid import advance

def main(args):
    repo=Path(__file__).resolve().parents[1];args.out.mkdir(parents=True,exist_ok=False);begin=time.time();reference=repo/'results/raar-circles/render/case-1051-animation.npz';(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());data=np.load(reference);z=data['states'];p=data['params'];h=args.step
    origin=np.floor((z.min(0)-args.domain_pad)/h)*h;shape=np.ceil((z.max(0)+args.domain_pad-origin)/h).astype(np.int64);cells=int(np.prod(shape));assert cells<args.max_grid_cells,(shape,cells)
    ijk=np.floor((z-origin)/h).astype(np.int64);keys=np.unique((ijk[:,0]*shape[1]+ijk[:,1])*shape[2]+ijk[:,2]);mask=initial_mask(keys,shape,args.initial_pad)
    windows=json.loads((repo/'results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json').read_text())
    for c,B,r in zip(windows['centers'],windows['bases'],windows['radii']):
        radius=np.sum(abs(np.array(B)*np.array(r)[None,:]),axis=1);a=np.floor((np.array(c)-radius-origin)/h).astype(int);b=np.floor((np.array(c)+radius-origin)/h).astype(int)
        assert (a>=0).all() and (b<shape).all()
        for i in range(a[0],b[0]+1):
            for j in range(a[1],b[1]+1):
                for k in range(a[2],b[2]+1):mark(mask,(i*shape[1]+j)*shape[2]+k)
    frontier=occupied_keys(mask);initial=len(frontier);history=[];flag=5;key=-1;position=0
    for generation in range(args.generations):
        new,flag,position,key,axis,splits,leaves,max_remainder=advance(frontier,mask,shape,origin,h,p,args.buffer,args.max_new,args.curvature_tolerance,args.max_depth)
        row=dict(subdivision_splits=int(splits),subdivision_leaves=int(leaves),largest_leaf_remainder_over_grid_step=float(max_remainder),generation=generation,sources=len(frontier),processed=int(position),added=len(new),status=int(flag),minimum_sampled_axis=float(axis),seconds=time.time()-begin);history.append(row);print(json.dumps(row),flush=True);(args.out/'progress.json').write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8')
        if flag:break
        frontier=new
        if not len(frontier):break
    occupied=occupied_keys(mask);np.savez_compressed(args.out/'proposal.npz',params=p,origin=origin,shape=shape,step=h,occupied=occupied,pending=frontier,initial_reference_keys=keys)
    result=dict(curvature_tolerance=args.curvature_tolerance,max_depth=args.max_depth,adaptive_source_sha256=hashlib.sha256(Path(__file__).with_name('raar_adaptive_grid.py').read_bytes()).hexdigest(),status='numerically_closed' if flag==0 and not len(frontier) else 'not_closed',reason={0:'Iteration limit or sampled closure',1:'Occupied cube intersects B projection axis',2:'Undefined or nonfinite sampled update',3:'Proposed image bounds leave declared bounding domain; actual escape is not established',4:'New-cell budget exceeded',5:'No generations attempted'}[int(flag)],step=h,initial_pad=args.initial_pad,image_buffer=args.buffer,image_shape=args.image_shape,domain_pad=args.domain_pad,grid_cells=cells,initial_occupied=initial,final_occupied=len(occupied),history=history,failing_key=int(key),failing_position=int(position),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),scope='Numerical forward-closure proposal from box vertices and centers. Zonotope remainder is sampled, not an interval bound. No invariant-region theorem, attraction proof or historical novelty claim. All attempted resolutions and failures are retained.')
    (args.out/'raar_adaptive_grid.py').write_bytes(Path(__file__).with_name('raar_adaptive_grid.py').read_bytes());(args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:v for k,v in result.items() if k!='history'}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--curvature-tolerance',type=float,default=.2);p.add_argument('--max-depth',type=int,default=3);p.add_argument('--out',type=Path,required=True);p.add_argument('--step',type=float,required=True);p.add_argument('--initial-pad',type=int,default=2);p.add_argument('--buffer',type=int,default=1);p.add_argument('--image-shape',choices=['box','zonotope'],default='box');p.add_argument('--domain-pad',type=float,default=.35);p.add_argument('--generations',type=int,default=80);p.add_argument('--max-new',type=int,default=2000000);p.add_argument('--max-grid-cells',type=int,default=600000000);main(p.parse_args())
