"""Enumerate every exact-depth preimage of infinity numerically, then compare independent forward basin labels."""
import argparse
import hashlib
import json
import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import distance_transform_edt
from search import RAW,OUTPUT,roots,plane
from selected import D,K,CASE,APPROX,BOUNDS

def refined():
    from search import critical_iter
    a=APPROX
    for _ in range(6):
        f,v=critical_iter(a,D,K);a-=f/v
    return a

def sphere(z):
    r2=np.abs(z)**2
    return np.stack((2*z.real/(1+r2),2*z.imag/(1+r2),(r2-1)/(1+r2)),axis=1)

def inverse(a,depth):
    z=(a/D)**(1/(D-1))*np.exp(2j*np.pi*np.arange(D-1)/(D-1))
    stats=[]
    for level in range(2,depth+1):
        new=np.empty(len(z)*D,complex);normres=0.;native=0.
        # Batched companion eigenproblems enumerate all branches without using the forward classifier.
        for start in range(0,len(z),1600):
            w=z[start:start+1600];n=len(w)
            companion=np.zeros((n,D,D),dtype=complex)
            companion[:,np.arange(1,D),np.arange(D-1)]=1
            companion[:,0,0]=D*w/(D-1)
            companion[:,0,-1]=-a*(w-1)/(D-1)
            u=np.linalg.eigvals(companion)
            order=np.argsort(np.angle(u),axis=1)
            u=np.take_along_axis(u,order,axis=1)
            residual=(D-1)*u**D-D*w[:,None]*u**(D-1)+a*(w[:,None]-1)
            scale=(D-1)*abs(u)**D+D*abs(w[:,None])*abs(u)**(D-1)+abs(a*(w[:,None]-1))
            normres=max(normres,float(np.max(abs(residual)/scale)))
            f=((D-1)*u**D-a)/(D*u**(D-1)-a)
            native=max(native,float(np.max(abs(f-w[:,None]))))
            new[D*start:D*(start+n)]=u.ravel()
        stats.append(dict(first_hit_infinity=level,count=len(new),max_relative_equation_residual=normres,max_forward_residual=native))
        z=new
        print('inverse level',level,len(z),flush=True)
    return z,stats

def compare(z,labels,box):
    ny,nx=labels.shape;lo,hi,bot,top=box;dx=(hi-lo)/nx;dy=(top-bot)/ny
    mask=(labels>=0)
    boundary=np.zeros_like(mask)
    pair=(labels[:,1:]!=labels[:,:-1])&mask[:,1:]&mask[:,:-1]
    boundary[:,1:]|=pair;boundary[:,:-1]|=pair
    pair=(labels[1:,:]!=labels[:-1,:])&mask[1:,:]&mask[:-1,:]
    boundary[1:,:]|=pair;boundary[:-1,:]|=pair
    iy,ix=np.where(boundary)
    rng=np.random.default_rng(22712)
    if len(ix)>160000:
        selection=rng.choice(len(ix),160000,replace=False);ix=ix[selection];iy=iy[selection]
    sample=np.c_[ix+.5,iy+.5]
    inv=np.c_[(z.real-lo)/dx,(top-z.imag)/dy]
    sel=(inv[:,0]>=1)&(inv[:,0]<nx-1)&(inv[:,1]>=1)&(inv[:,1]<ny-1)
    candidates=inv[sel]
    if len(candidates)>200000:candidates=candidates[rng.choice(len(candidates),200000,replace=False)]
    distance=distance_transform_edt(~boundary)
    # Pixel-center residual, rounded to closest pixel for inverse-to-classified-boundary.
    d1=distance[np.floor(candidates[:,1]).astype(int),np.floor(candidates[:,0]).astype(int)]
    d2=cKDTree(inv).query(sample)[0]
    return dict(screen=list(box),shape=list(labels.shape),known_fraction=float(mask.mean()),pixel_width=dx,pixel_height=dy,
        forward_interface_pixels=int(boundary.sum()),inverse_to_forward_samples=len(d1),forward_to_inverse_samples=len(d2),
        quantile_levels=[.5,.9,.99,1],inverse_to_forward_grid_pixels=np.quantile(d1,[.5,.9,.99,1]).tolist(),
        forward_to_inverse_grid_pixels=np.quantile(d2,[.5,.9,.99,1]).tolist(),
        inverse_fraction_within_two_pixels=float((d1<=2).mean()),forward_fraction_within_two_pixels=float((d2<=2).mean()))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--depth',type=int,default=8);args=ap.parse_args()
    a=refined();z,stats=inverse(a,args.depth)
    points=sphere(z);f=RAW/'case-6902-inverse-sphere.npz';np.savez_compressed(f,points=points)
    planar=RAW/'case-6902-inverse-plane.npz';np.savez_compressed(planar,points=np.c_[z.real,z.imag])
    l,c=plane(D,a,K,roots(D,a),*BOUNDS,1480,960,330)
    np.savez_compressed(RAW/'case-6902-forward.npz',labels=l,steps=c,bounds=BOUNDS)
    control=compare(z,l,BOUNDS)
    def filedata(p,key,shape):
        return dict(cloud_relative='boundary/'+p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),key=key,shape=shape,bytes=p.stat().st_size)
    data=dict(case=CASE,degree=D,period=K,parameter_floating=[a.real,a.imag],critical_equation='N_a^9(0)=0',
        all_roots=[[p.real,p.imag] for p in roots(D,a)],first_infinity_depth=args.depth,levels=stats,
        sphere=filedata(f,'points',list(points.shape)),plane=filedata(planar,'points',[len(z),2]),
        forward=filedata(RAW/'case-6902-forward.npz','labels',list(l.shape)),
        sphere_max_unit_norm_defect=float(abs(np.linalg.norm(points,axis=1)-1).max()),screen_control=control,
        labels={'roots':'0 through 5 by argument in [0,2 pi)','phases':'6 plus the index relative to 0 at iterate zero','unclassified':-1},
        outcome_counts=np.bincount(l.ravel()+1,minlength=D+K+1).tolist())
    (OUTPUT/'inverse-check.json').write_text(json.dumps(data,indent=2)+'\n')
    print(control,flush=True)

if __name__=='__main__':main()
