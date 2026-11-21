"""One seventh of the angular itinerary in a fixed C7 chart, and a chart-composition test."""
import argparse
import json
import math
import numpy as np
from numba import njit
from native_section import RAW,packed,trace
from return_tube import CENTER,WIDTH,ZLO,ZHI,graph,grid

@njit(cache=True,nogil=True)
def reduced(xs,zs,h,p):
    phi=8*np.pi/7;c,s=np.cos(phi),np.sin(phi);n=len(zs)
    out=np.zeros((n,4))
    for i in range(n):
        hits=trace(np.array([xs[i]*c,xs[i]*s,zs[i]]),h,round(7.5/h),p)
        found=0
        for row in hits:
            if row[0]>5.7 and row[0]<7.5 and row[1]>.77 and row[1]<.85 and row[2]>.945 and row[2]<1.04:
                out[i,:3]=row
                found+=1
        out[i,3]=found
    return out

def main(n,levels,h):
    g=graph();z,d=np.meshgrid(np.linspace(ZLO,ZHI,n),np.linspace(-WIDTH,WIDTH,levels),indexing='ij')
    xs=g(z-CENTER)+d;p=packed()
    a=reduced(xs.ravel(),z.ravel(),h,p).reshape(n,levels,4)
    assert np.all(a[:,:,3]==1)
    zp=a[:,:,2];dp=a[:,:,1]-g(zp-CENTER)
    slopes=np.diff(zp[:,levels//2])/np.diff(z[:,0])
    base=np.arange(0,n,max(1,n//70));x=xs[base,0];y=z[base,0];full=grid(x,y,h,p)
    t=np.zeros(len(x))
    for k in range(7):
        ry=reduced(x,y,h,p);assert np.all(ry[:,3]==1);t+=ry[:,0];x=ry[:,1];y=ry[:,2]
    q=dict(n=n,levels=levels,h=h,window=[5.7,7.5],rotation_applied_to_input=8*math.pi/7,
           z_image=[zp.min(),zp.max()],delta_image=[dp.min(),dp.max()],
           times=[a[:,:,0].min(),a[:,:,0].max()],slopes=[slopes.min(),slopes.max()],
           centerline_turns=z[1:-1,0][slopes[1:]*slopes[:-1]<0].tolist(),
           seven_reduced_vs_fifth_native_max_error_t_x_z=np.max(np.abs(np.c_[t,x,y]-full[:,:3]),axis=0).tolist())
    print(q,flush=True)
    name=f'reduced-n{n}-l{levels}-h{h:g}'
    np.savez_compressed(RAW/(name+'.npz'),z=z,d=d,x=xs,returns=a,coefficients=g.coef)
    (RAW/(name+'.json')).write_text(json.dumps(q,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n',type=int,default=2049);ap.add_argument('--levels',type=int,default=17)
    ap.add_argument('--h',type=float,default=.00625)
    a=ap.parse_args();main(a.n,a.levels,a.h)
