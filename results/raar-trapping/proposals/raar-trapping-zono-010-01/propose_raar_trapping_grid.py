"""Numerical cubical forward-closure proposals, never interval trapping proofs.

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

@njit(cache=True)
def initial_mask(keys,shape,pad):
    mask=np.zeros(shape[0]*shape[1]*shape[2],np.bool_);ny,nz=shape[1],shape[2]
    for key in keys:
        ix=key//(ny*nz);iy=(key//nz)%ny;iz=key%nz
        for a in range(max(0,ix-pad),min(shape[0],ix+pad+1)):
            for b in range(max(0,iy-pad),min(ny,iy+pad+1)):
                for c in range(max(0,iz-pad),min(nz,iz+pad+1)):mask[(a*ny+b)*nz+c]=True
    return mask

@njit(cache=True)
def advance(frontier,mask,shape,origin,h,p,buffer,max_new,zonotope):
    ny,nz=shape[1],shape[2];new=[np.int64(0) for _ in range(0)];minimum=1e100;sn,co=np.sin(p[2]),np.cos(p[2]);corner=np.empty(3)
    for pos,key in enumerate(frontier):
        ijk=np.array([key//(ny*nz),(key//nz)%ny,key%nz]);lo=origin+h*ijk;hi=lo+h
        # dx and h_B are independent coordinate blocks, so this detects an axis crossing.
        hlo=co*(lo[1]-p[4])+sn*(lo[2]-p[5]);hhi=co*(hi[1]-p[4])+sn*(hi[2]-p[5])
        if lo[0]<=p[3]<=hi[0] and hlo<=0<=hhi:return np.array(new),1,pos,key,minimum
        lower=np.full(3,np.inf);upper=-lower;midpoint=lo+.5*h;center,J,_,_=step(midpoint,p);remainder=np.zeros(3)
        for k in range(9):
            if k<8:
                corner[0]=lo[0]+h*(k&1);corner[1]=lo[1]+h*((k>>1)&1);corner[2]=lo[2]+h*((k>>2)&1)
            else:corner=lo+.5*h
            y,flag,axis=value(corner,p);minimum=min(minimum,axis)
            if flag or not np.isfinite(y).all():return np.array(new),2,pos,key,minimum
            lower=np.minimum(lower,y);upper=np.maximum(upper,y)
            remainder=np.maximum(remainder,abs(y-center-J@(corner-midpoint)))
        directions=np.zeros((15,3));support=np.zeros(15)
        if zonotope:
            # Sampled curvature only: this is still a proposal, not an interval bound.
            remainder=1.1*remainder+1e-12;radius=.5*h*np.sum(abs(J),axis=1)+remainder;lower=center-radius;upper=center+radius
            directions[:3]=np.eye(3);k=3
            for i in range(3):
                for j in range(i+1,3):directions[k]=np.cross(J[:,i],J[:,j]);k+=1
            for i in range(3):
                for j in range(3):directions[k]=np.cross(J[:,i],np.eye(3)[j]);k+=1
            for k in range(15):
                length=np.linalg.norm(directions[k])
                if length>1e-15:directions[k]/=length
                support[k]=.5*h*np.sum(abs(directions[k]@J))+np.sum(abs(directions[k])*(remainder+(.5+buffer)*h))
        a=np.floor((lower-origin)/h).astype(np.int64)-buffer;b=np.floor((upper-origin)/h).astype(np.int64)+buffer
        if np.any(a<0) or np.any(b>=shape):return np.array(new),3,pos,key,minimum
        for ix in range(a[0],b[0]+1):
            for iy in range(a[1],b[1]+1):
                for iz in range(a[2],b[2]+1):
                    target=(ix*ny+iy)*nz+iz
                    if zonotope:
                        displacement=origin+h*(np.array([ix,iy,iz])+.5)-center;intersects=True
                        for k in range(15):
                            if abs(np.dot(directions[k],displacement))>support[k]:intersects=False;break
                        if not intersects:continue
                    if not mask[target]:mask[target]=True;new.append(target)
        if len(new)>max_new:return np.array(new),4,pos,key,minimum
    return np.array(new),0,len(frontier),-1,minimum

def main(args):
    repo=Path(__file__).resolve().parents[1];args.out.mkdir(parents=True,exist_ok=False);begin=time.time();reference=repo/'results/raar-circles/render/case-1051-animation.npz';data=np.load(reference);z=data['states'];p=data['params'];h=args.step
    origin=np.floor((z.min(0)-args.domain_pad)/h)*h;shape=np.ceil((z.max(0)+args.domain_pad-origin)/h).astype(np.int64);cells=int(np.prod(shape));assert cells<args.max_grid_cells,(shape,cells)
    ijk=np.floor((z-origin)/h).astype(np.int64);keys=np.unique((ijk[:,0]*shape[1]+ijk[:,1])*shape[2]+ijk[:,2]);mask=initial_mask(keys,shape,args.initial_pad)
    windows=json.loads((repo/'results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json').read_text())
    for c,B,r in zip(windows['centers'],windows['bases'],windows['radii']):
        radius=np.sum(abs(np.array(B)*np.array(r)[None,:]),axis=1);a=np.floor((np.array(c)-radius-origin)/h).astype(int);b=np.floor((np.array(c)+radius-origin)/h).astype(int)
        assert (a>=0).all() and (b<shape).all()
        for i in range(a[0],b[0]+1):
            for j in range(a[1],b[1]+1):mask[(i*shape[1]+j)*shape[2]+a[2]:(i*shape[1]+j)*shape[2]+b[2]+1]=True
    frontier=np.flatnonzero(mask);initial=len(frontier);history=[];flag=5;key=-1;position=0
    for generation in range(args.generations):
        new,flag,position,key,axis=advance(frontier,mask,shape,origin,h,p,args.buffer,args.max_new,args.image_shape=='zonotope')
        row=dict(generation=generation,sources=len(frontier),processed=int(position),added=len(new),status=int(flag),minimum_sampled_axis=float(axis),seconds=time.time()-begin);history.append(row);print(json.dumps(row),flush=True)
        if flag:break
        frontier=new
        if not len(frontier):break
    occupied=np.flatnonzero(mask);np.savez_compressed(args.out/'proposal.npz',params=p,origin=origin,shape=shape,step=h,occupied=occupied,pending=frontier,initial_reference_keys=keys)
    result=dict(status='numerically_closed' if flag==0 and not len(frontier) else 'not_closed',reason={0:'Iteration limit or sampled closure',1:'Occupied cube intersects B projection axis',2:'Undefined or nonfinite sampled update',3:'Sampled image leaves declared bounding domain',4:'New-cell budget exceeded',5:'No generations attempted'}[int(flag)],step=h,initial_pad=args.initial_pad,image_buffer=args.buffer,image_shape=args.image_shape,domain_pad=args.domain_pad,grid_cells=cells,initial_occupied=initial,final_occupied=len(occupied),history=history,failing_key=int(key),failing_position=int(position),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),scope='Numerical forward-closure proposal from box vertices and centers. Zonotope remainder is sampled, not an interval bound. No invariant-region theorem, attraction proof or historical novelty claim. All attempted resolutions and failures are retained.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:v for k,v in result.items() if k!='history'}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--step',type=float,required=True);p.add_argument('--initial-pad',type=int,default=2);p.add_argument('--buffer',type=int,default=1);p.add_argument('--image-shape',choices=['box','zonotope'],default='box');p.add_argument('--domain-pad',type=float,default=.35);p.add_argument('--generations',type=int,default=80);p.add_argument('--max-new',type=int,default=2000000);p.add_argument('--max-grid-cells',type=int,default=600000000);main(p.parse_args())
