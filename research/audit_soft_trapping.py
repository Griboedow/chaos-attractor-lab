"""Check a fixed cubical proposal with outward soft-map image bounds.

Every cube, including its boundary, must map with a positive coordinate
margin into the recorded union. Subdivision tightens enclosures; it never
removes source points. An unresolved leaf is a failed check, not an escape.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit
from search_soft_constellations import step
from soft_interval import remainder
from propose_raar_trapping_bitset import marked,mark
from raar_interval_scalar_grid import coordinate_interval,raster_parameters,target_separated


@njit(cache=False)
def audit(occupied,mask,shape,origin,h,beta,tau,points,margin,max_depth,max_leaves):
    ny,nz=shape[1],shape[2];leaves=0;splits=0;axis_min=1e100;deepest=0;max_error=0.
    for pos,key in enumerate(occupied):
        ix,iy,iz=key//(ny*nz),(key//nz)%ny,key%nz
        stack=[(ix,iy,iz,0)]
        while stack:
            sx,sy,sz,depth=stack.pop();indices=np.array([sx,sy,sz]);width=h/(2**depth)
            box=np.empty((3,2));mid=np.empty(3)
            for j in range(3):
                box[j,0]=coordinate_interval(origin[j],width,float(indices[j]))[0]
                box[j,1]=coordinate_interval(origin[j],width,float(indices[j]+1))[1]
                mid[j]=origin[j]+width*(indices[j]+.5)
                assert box[j,0]<=mid[j]<=box[j,1]
            center,J,r,peak,entropy,flag=step(mid,beta,points,tau)
            if flag:return 1,pos,key,-1,leaves,splits,axis_min,deepest,max_error
            radii,errors,axis,flag=remainder(box,mid,center,J,beta,tau,points,margin)
            missing=-1;contained=not flag
            if contained:
                a,b,directions,support=raster_parameters(center,J,radii,errors,origin,h)
                if np.any(a<0) or np.any(b>=shape):contained=False
                else:
                    for tx in range(a[0],b[0]+1):
                        if not contained:break
                        for ty in range(a[1],b[1]+1):
                            if not contained:break
                            for tz in range(a[2],b[2]+1):
                                target=(tx*ny+ty)*nz+tz
                                if marked(mask,target):continue
                                if not target_separated(np.array([tx,ty,tz]),origin,h,center,directions,support):
                                    contained=False;missing=target;break
            if not contained:
                if depth>=max_depth:return 2,pos,key,missing,leaves,splits,axis_min,deepest,max_error
                splits+=1
                for k in range(8):stack.append((2*sx+(k&1),2*sy+((k>>1)&1),2*sz+((k>>2)&1),depth+1))
                continue
            leaves+=1;axis_min=min(axis_min,axis);deepest=max(deepest,depth);max_error=max(max_error,errors.max())
            if leaves>=max_leaves:return 3,pos,key,-1,leaves,splits,axis_min,deepest,max_error
    return 0,len(occupied),-1,-1,leaves,splits,axis_min,deepest,max_error


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();sources={}
    for name in [Path(__file__).name,'soft_interval.py','raar_interval_scalar.py','raar_interval_scalar_grid.py','raar_geometric_scalar.py','search_soft_constellations.py','search_raar_constellations.py','propose_raar_trapping_bitset.py']:
        data=Path(__file__).with_name(name).read_bytes();(args.out/name).write_bytes(data);sources[name]=hashlib.sha256(data).hexdigest()
    with np.load(args.source) as f:
        beta,tau=f['params'];points=f['points'];occupied=f['occupied'];origin=f['origin'];shape=f['shape'];h=float(f['step'])
    assert np.all(np.diff(occupied)>0);mask=np.zeros((int(np.prod(shape))+7)//8,np.uint8)
    for k in occupied:mark(mask,k)
    history=[];minimum=1e100;leaves=0;splits=0;processed=0;deepest=0;flag=0
    for offset in range(0,len(occupied),args.batch):
        block=occupied[offset:offset+args.batch]
        flag,pos,key,missing,nleaf,nsplit,axis,depth,error=audit(block,mask,shape,origin,h,beta,tau,points,args.margin,args.max_depth,args.max_leaves)
        minimum=min(minimum,axis);leaves+=nleaf;splits+=nsplit;processed+=pos;deepest=max(deepest,depth)
        row=dict(offset=offset,count=len(block),status=int(flag),processed=int(pos),failing_key=int(key),first_unresolved_target=int(missing),leaves=int(nleaf),splits=int(nsplit),minimum_axis_lower_bound=float(axis),maximum_depth=int(depth),largest_error=float(error),seconds=time.time()-begin)
        history.append(row);print(json.dumps(row),flush=True)
        (args.out/'progress.json').write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8',newline='\n')
        if flag:break
    result=dict(status='all_cubes_passed' if flag==0 and processed==len(occupied) else 'not_certified',source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),sources=sources,beta=float(beta),tau=float(tau),step=h,margin=args.margin,cubes=len(occupied),processed=processed,leaves=leaves,splits=splits,minimum_axis_lower_bound=float(minimum),maximum_depth=deepest,history=history,seconds=time.time()-begin,scope='Outward interval check of the fixed union and a positive coordinate margin. Requires source/provenance and independent replay audit before publication as a trapping certificate. Does not establish chaos or historical uniqueness.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({k:v for k,v in result.items() if k!='history'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--margin',type=float,default=2.**-35);p.add_argument('--max-depth',type=int,default=8);p.add_argument('--max-leaves',type=int,default=2000000);p.add_argument('--batch',type=int,default=1000);main(p.parse_args())
