"""Independent Poincare multipliers, C7 section separation and perturbed basin tests."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.spatial import cKDTree
from branches import polar_ring,rotation
from field import integrate,ring_tangent


def section_matrix(row,p,scale):
    T=row['period'];v=np.array(row['nodes'][0]);p=p.copy();p[8]=row['e']
    f=lambda t,x:polar_ring(x,p)
    normal=f(0,v);normal/=np.linalg.norm(normal);one=np.cross(np.eye(3)[np.argmin(abs(normal))],normal);one/=np.linalg.norm(one)
    B=np.column_stack([one,np.cross(normal,one)])
    def event(t,s):return np.dot(s-v,normal)
    event.direction=1.
    records=[]
    for method,rtol in [('DOP853',2e-12),('RK45',2e-11)]:
        def ret(a):
            sol=solve_ivp(f,(0,1.01*T),a,method=method,rtol=rtol,atol=rtol/100,events=event)
            ts=sol.t_events[0];i=np.argmin(abs(ts-T))
            if abs(ts[i]-T)>.005*T:raise RuntimeError('Wrong return')
            return B.T@(sol.y_events[0][i]-v)
        center=ret(v)
        for eps in [3e-6,3e-7,1e-7]:
            h=eps*scale
            mat=np.column_stack([(ret(v+h*d)-ret(v-h*d))/(2*h) for d in B.T]);eig=np.linalg.eigvals(mat)
            records.append(dict(method=method,rtol=rtol,displacement_rms=eps,matrix=mat.tolist(),
                      dominant=float(eig[np.argmax(abs(eig))].real),center_residual_native=float(np.linalg.norm(center))))
    return records


def orbit_sections(row,p):
    p=p.copy();p[8]=row['e'];T=row['period'];v=row['nodes'][0]
    result=[]
    for method,rtol in [('DOP853',8e-13),('RK45',2e-12)]:
        sol=solve_ivp(lambda _,z:polar_ring(z,p),(0,T),v,method=method,rtol=rtol,atol=rtol/100,dense_output=True)
        t=np.linspace(0,T,8193);z=sol.sol(t).T;c=[]
        for k in range(7):
            R=rotation(2*np.pi*k/7);y=(z@R.T)[:,1];idx=np.flatnonzero((y[:-1]<0)&(y[1:]>=0));pts=[]
            for i in idx:
                crossing=brentq(lambda tt:(R@sol.sol(tt))[1],t[i],t[i+1])
                pts.append((R@sol.sol(crossing))[[0,2]])
            c.append(np.array(pts))
        values=[]
        for i,a in enumerate(c):
            for b in c[:i]:
                d=np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2)
                values.append([float(d.min()),float(max(d.min(0).max(),d.min(1).max()))])
        result.append(dict(method=method,minimum_distinct_section_point_distance=min(x[0] for x in values),
                        minimum_Hausdorff_distance=min(x[1] for x in values),counts=[len(x) for x in c]))
    return result


def basin(row,p):
    p=p.copy();p[8]=row['e'];T=row['period'];z=np.load(row['path'])['points']
    tree=cKDTree(z);rms=np.linalg.norm(np.std(z,axis=0));rng=np.random.default_rng(151927)
    def distances(trace):
        distance=[]
        for k in range(7):
            loc=trace@rotation(2*np.pi*k/7);_,near=tree.query(loc)
            possibilities=[]
            for shift in [-1,0]:
                a=z[(near+shift)%len(z)];v=z[(near+shift+1)%len(z)]-a
                s=np.clip(np.einsum('ij,ij->i',loc-a,v)/np.einsum('ij,ij->i',v,v),0,1)
                possibilities.append(np.linalg.norm(loc-a-s[:,None]*v,axis=1))
            distance.append(float(np.quantile(np.minimum(*possibilities),.95)/rms))
        return distance
    trials=[];broad=[]
    for k in range(7):
        for h in (.0125,.00625):
            for radius in [.0001,.001]:
                j=int(rng.integers(len(z)));move=rng.normal(size=3);start=rotation(2*np.pi*k/7)@z[j]+rms*radius*move/np.linalg.norm(move)
                _,tail=integrate(p,start,h,round(80*T/h),round(2*T/h),1,max(1,round(.07/h)),ring_tangent)
                d=distances(tail)
                trials.append(dict(copy=k,h=h,radius_rms=radius,initial=start.tolist(),nearest_copy=int(np.argmin(d)),distance_rms=min(d),distances=d))
    for i in range(12):
        start=rng.uniform(-4,4,3);h=.0125
        rates,tail=integrate(p,start,h,round(1200/h),round(350/h),6,6,ring_tangent)
        d=distances(tail)
        broad.append(dict(initial=start.tolist(),rates=rates.tolist(),nearest_copy=int(np.argmin(d)),distance_rms=min(d)))
    return dict(local=trials,broad=broad)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--part',choices=['basins','sections'],required=True)
    a=ap.parse_args();d=json.loads(a.data.read_text());cases=d['representatives'];p=np.array(cases['stable']['parameters'])
    a.out.mkdir(exist_ok=False,parents=True)
    if a.part=='basins':
        (a.out/'basins.json').write_text(json.dumps(basin(cases['stable'],p),indent=2)+'\n')
        return
    scale=.89145
    targets={**cases,**list(d['events'].values())[-1]}
    for name,row in targets.items():
        result=dict(poincare=section_matrix(row,p,scale),section=orbit_sections(row,p))
        (a.out/f'{name}.json').write_text(json.dumps(result,indent=2)+'\n')
        print(name,result['poincare'][-1]['dominant'],flush=True)


if __name__=='__main__':
    main()
