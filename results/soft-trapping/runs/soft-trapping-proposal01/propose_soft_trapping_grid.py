"""Sampled grid-closure proposals for the soft circle map, never proofs.

Corner residuals are inflated heuristically. A successful closure must be
rechecked with outward interval bounds before any invariant-set assertion.
"""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import njit
from search_soft_constellations import step
from propose_raar_trapping_bitset import marked, mark, occupied_keys, initial_mask


@njit(cache=True)
def value(x, beta, points, tau):
    logits=(2*points@x-np.sum(points**2,axis=1))/tau
    weights=np.exp(logits-logits.max()); weights/=weights.sum()
    m=weights@points; v=2*m-x; r=np.hypot(v[0],v[1])
    if r<1e-12:return x,r,1
    return beta*x+beta*np.array([v[0]/r,v[1]/r,0.])+(1-2*beta)*m,r,0


@njit(cache=True)
def advance(frontier,mask,shape,origin,h,beta,tau,points,max_new,slack,curvature,split,cut):
    ny,nz=shape[1],shape[2];new=[np.int64(0) for _ in range(0)]
    minimum=1e100; flip_margin=1e100; largest_remainder=0.
    for pos,key in enumerate(frontier):
        ijk=np.array([key//(ny*nz),(key//nz)%ny,key%nz]);mid=origin+h*(ijk+.5)
        center,J,r,peak,entropy,flag=step(mid,beta,points,tau)
        if flag:return np.array(new),1,pos,key,minimum,flip_margin,largest_remainder
        minimum=min(minimum,r);remainder=np.zeros(3)
        source_center=np.dot(split,mid)-cut;source_radius=.5*h*np.sum(abs(split))
        if abs(source_center)<=source_radius:return np.array(new),4,pos,key,minimum,flip_margin,largest_remainder
        for corner in range(8):
            delta=.5*h*np.array([2*(corner&1)-1,2*((corner>>1)&1)-1,2*((corner>>2)&1)-1])
            y,r,flag=value(mid+delta,beta,points,tau);minimum=min(minimum,r)
            if flag:return np.array(new),1,pos,key,minimum,flip_margin,largest_remainder
            remainder=np.maximum(remainder,abs(y-center-J@delta))
        remainder=curvature*remainder+slack
        largest_remainder=max(largest_remainder,np.max(remainder))
        image_center=np.dot(split,center)-cut
        image_radius=.5*h*np.sum(abs(split@J))+np.sum(abs(split)*remainder)
        flip_margin=min(flip_margin,-np.sign(source_center)*image_center-image_radius)
        if -np.sign(source_center)*image_center<=image_radius:
            return np.array(new),5,pos,key,minimum,flip_margin,largest_remainder
        radius=.5*h*np.sum(abs(J),axis=1)+remainder
        lower=np.floor((center-radius-origin)/h).astype(np.int64)
        upper=np.floor((center+radius-origin)/h).astype(np.int64)
        if np.any(lower<0) or np.any(upper>=shape):
            return np.array(new),2,pos,key,minimum,flip_margin,largest_remainder
        directions=np.zeros((15,3));support=np.zeros(15);directions[:3]=np.eye(3);k=3
        for i in range(3):
            for j in range(i+1,3):directions[k]=np.cross(J[:,i],J[:,j]);k+=1
        for i in range(3):
            for j in range(3):directions[k]=np.cross(J[:,i],np.eye(3)[j]);k+=1
        for k in range(15):
            length=np.linalg.norm(directions[k])
            if length>1e-15:directions[k]/=length
            support[k]=.5*h*np.sum(abs(directions[k]@J))+np.sum(abs(directions[k])*(remainder+.5*h))
        for ix in range(lower[0],upper[0]+1):
            for iy in range(lower[1],upper[1]+1):
                for iz in range(lower[2],upper[2]+1):
                    target=(ix*ny+iy)*nz+iz
                    if marked(mask,target):continue
                    displacement=origin+h*(np.array([ix,iy,iz])+.5)-center;intersects=True
                    for k in range(15):
                        if abs(np.dot(directions[k],displacement))>support[k]:intersects=False;break
                    if intersects:mark(mask,target);new.append(target)
        if len(new)>max_new:return np.array(new),3,pos,key,minimum,flip_margin,largest_remainder
    return np.array(new),0,len(frontier),-1,minimum,flip_margin,largest_remainder


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    for name in [Path(__file__).name,'search_soft_constellations.py','search_raar_constellations.py','propose_raar_trapping_bitset.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:z=f['states'];beta,tau=f['params'];points=f['points']
    # Case 29 split fitted before this proposal to the frozen long-run cloud.
    split=np.array([-.2375212291113019,-.9692369549104499,-.06452434391277834])
    split_center=np.array([.6946711958320783,-.2263318555502484,.49339841102498244])
    cut=-.038912546427851594+split@split_center
    h=args.step;origin=np.floor((z.min(0)-args.domain_pad)/h)*h
    shape=np.ceil((z.max(0)+args.domain_pad-origin)/h).astype(np.int64)
    assert int(np.prod(shape))<args.max_grid_cells
    ijk=np.floor((z-origin)/h).astype(np.int64);keys=np.unique((ijk[:,0]*shape[1]+ijk[:,1])*shape[2]+ijk[:,2])
    mask=initial_mask(keys,shape,args.initial_pad);frontier=occupied_keys(mask);initial=len(frontier);history=[];flag=6;key=-1
    for generation in range(args.generations):
        new,flag,pos,key,r,margin,remainder=advance(frontier,mask,shape,origin,h,beta,tau,points,args.max_new,args.slack,args.curvature,split,cut)
        row=dict(generation=generation,sources=len(frontier),processed=int(pos),added=len(new),status=int(flag),minimum_sampled_radius=float(r),minimum_sampled_enclosure_flip_margin=float(margin),maximum_sampled_remainder=float(remainder),seconds=time.time()-begin)
        history.append(row);print(json.dumps(row),flush=True)
        (args.out/'progress.json').write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8',newline='\n')
        if flag:break
        frontier=new
        if not len(frontier):break
    occupied=occupied_keys(mask)
    np.savez_compressed(args.out/'proposal.npz',params=[beta,tau],points=points,origin=origin,shape=shape,step=h,occupied=occupied,pending=frontier,split=split,cut=cut)
    result=dict(status='numerically_closed' if flag==0 and not len(frontier) else 'not_closed',reason={0:'Sampled closure or generation limit',1:'Undefined sampled map',2:'Image enclosure leaves proposed domain',3:'New-cell budget exceeded',4:'Source box meets split plane',5:'Sampled image enclosure fails lobe alternation',6:'No generations attempted'}[int(flag)],source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),step=h,slack=args.slack,curvature=args.curvature,initial_pad=args.initial_pad,domain_pad=args.domain_pad,shape=shape.tolist(),initial_occupied=initial,final_occupied=len(occupied),failing_key=int(key),history=history,scope='Sampled closure proposal only. Corner residual inflation is not a rigorous remainder bound. No proof of trapping, alternation, chaos, attraction, or novelty.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k!='history'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--step',type=float,required=True)
    p.add_argument('--initial-pad',type=int,default=1);p.add_argument('--domain-pad',type=float,default=.25);p.add_argument('--curvature',type=float,default=2.);p.add_argument('--slack',type=float,default=1e-8)
    p.add_argument('--max-new',type=int,default=1000000);p.add_argument('--max-grid-cells',type=int,default=300000000);p.add_argument('--generations',type=int,default=100);main(p.parse_args())
