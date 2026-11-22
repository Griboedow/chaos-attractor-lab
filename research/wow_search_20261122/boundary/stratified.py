"""Produce a spatially stratified inverse Julia sample, independent of all forward classifications."""
import hashlib,json
import numpy as np
from scipy.spatial import cKDTree
from search import ROOT,RAW,OUTPUT
from selected import D,K,BOUNDS
from inverse import refined,sphere,compare

def children(w,a):
    vals=np.empty((len(w),D),complex);largest_res=0.
    for s in range(0,len(w),2000):
        ww=w[s:s+2000];mat=np.zeros((len(ww),D,D),complex)
        mat[:,np.arange(1,D),np.arange(D-1)]=1;mat[:,0,0]=D*ww/(D-1);mat[:,0,-1]=-a*(ww-1)/(D-1)
        zz=np.linalg.eigvals(mat);idx=np.argsort(np.angle(zz),axis=1);zz=np.take_along_axis(zz,idx,axis=1);vals[s:s+len(ww)]=zz
        res=(D-1)*zz**D-D*ww[:,None]*zz**(D-1)+a*(ww[:,None]-1)
        denom=(D-1)*abs(zz)**D+D*abs(ww[:,None])*abs(zz)**(D-1)+abs(a*(ww[:,None]-1))
        largest_res=max(largest_res,np.max(abs(res)/denom))
    return vals.ravel(),float(largest_res)

def bins(z,spacing):
    pts=sphere(z)
    # A fixed XYZ voxel grid covers the whole sphere, including its poles, without a planar crop.
    ids=np.floor((pts+1)/spacing).astype('int64');n=round(2/spacing)+2
    return ids[:,0]*n*n+ids[:,1]*n+ids[:,2]

def run(depth=23,spacing=.004):
    a=refined();z=(a/D)**(1/(D-1))*np.exp(2j*np.pi*np.arange(D-1)/(D-1))
    allz=z.copy();depths=np.ones(len(z),np.uint8);stats=[]
    for d in range(2,depth+1):
        z,err=children(z,a);bid=bins(z,spacing);_,i=np.unique(bid,return_index=True);z=z[i]
        candidate=np.concatenate([allz,z]);ds=np.r_[depths,np.full(len(z),d,dtype=np.uint8)]
        ids=bins(candidate,spacing);_,g=np.unique(ids,return_index=True);allz=candidate[g];depths=ds[g]
        stats.append(dict(depth=d,retained=len(z),unique_total=len(allz),inverse_max_relative_residual=err))
        print(d,len(z),len(allz),flush=True)
    pts=sphere(allz)
    cloud=RAW/'case-6902-stratified-inverse-sphere.npz'
    np.savez_compressed(cloud,points=pts,first_infinity_depth=depths,plane=np.c_[allz.real,allz.imag])
    with np.load(RAW/'case-6902-forward.npz') as inp: control=compare(allz,inp['labels'],BOUNDS)
    report=dict(case=6902,cloud_relative='boundary/'+cloud.name,point_key='points',shape=list(pts.shape),dtype=str(pts.dtype),
                sha256=hashlib.sha256(cloud.read_bytes()).hexdigest(),voxel_spacing=spacing,final_depth=depth,levels=stats,
                screen_control=control,embedding='Unit Riemann sphere: (2 Re z, 2 Im z, |z|²−1)/(1+|z|²). No planar crop; no point interpolation. A representative is retained per occupied 3D voxel.')
    (OUTPUT/'stratified-inverse.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
if __name__=='__main__':run()
