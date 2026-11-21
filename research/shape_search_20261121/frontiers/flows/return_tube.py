"""Stress-test a two-dimensional tube for the fifth native Poincare return."""
import argparse
import json
import numpy as np
from numba import njit
from numpy.polynomial import Polynomial
from native_section import RAW,RESULT,P,packed,trace

ZLO=.98045
ZHI=.99945
WIDTH=3e-5
CENTER=.99

def graph():
    # Fit only the first of five observed bands, evaluated in centered coordinates.
    ev=np.load(RAW/'s000-h0.02-t90000.npz')['events']
    ev=ev[(ev[:,1]>.78)&(ev[:,1]<.83)&(ev[:,2]>.96)]
    return Polynomial(np.polynomial.polynomial.polyfit(ev[:,2]-CENTER,ev[:,1],5))

@njit(cache=True,nogil=True)
def grid(xs,zs,h,p):
    n=len(zs);out=np.zeros((n,5))
    for i in range(n):
        hits=trace(np.array([xs[i],0.,zs[i]]),h,round(51/h),p)
        if len(hits)<5:raise RuntimeError('Missing fifth hit')
        out[i,:3]=hits[4]
        # Positive normal speeds at the five crossings are independently assessed downstream.
        out[i,3]=np.min(np.diff(hits[:5,0].copy()))
        out[i,4]=len(hits)
    return out

def run(n=401,levels=13,h=.0125):
    g=graph()
    z,d=np.meshgrid(np.linspace(ZLO,ZHI,n),np.linspace(-WIDTH,WIDTH,levels),indexing='ij')
    xs=g(z-CENTER)+d
    ret=grid(xs.ravel(),z.ravel(),h,packed()).reshape((n,levels,5))
    zp=ret[:,:,2];dp=ret[:,:,1]-g(zp-CENTER)
    name=f'tube-n{n}-l{levels}-h{h:g}'
    np.savez_compressed(RAW/(name+'.npz'),z=z,d=d,x=xs,returns=ret,coefficients=g.coef)
    deriv=np.diff(zp[:,levels//2])/np.diff(z[:,levels//2])
    q=dict(n=n,levels=levels,h=h,z_domain=[ZLO,ZHI],transverse_halfwidth=WIDTH,
            graph_center=CENTER,graph_coefficients=g.coef.tolist(),
            z_image=[float(zp.min()),float(zp.max())],delta_image=[float(dp.min()),float(dp.max())],
            return_time=[float(ret[:,:,0].min()),float(ret[:,:,0].max())],
            least_time_between_native_hits=float(ret[:,:,3].min()),centerline_slopes=[float(deriv.min()),float(deriv.max())],
            centerline_turning_points=z[1:-1,0][deriv[1:]*deriv[:-1]<0].tolist())
    print(name,q,flush=True)
    (RAW/(name+'.json')).write_text(json.dumps(q,indent=2)+'\n')
    return q

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n',type=int,default=401);ap.add_argument('--levels',type=int,default=13)
    ap.add_argument('--h',type=float,default=.0125)
    a=ap.parse_args();run(a.n,a.levels,a.h)
