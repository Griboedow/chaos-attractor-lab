"""Locate the saddle-node and first flip and reconstruct full, unsymmetrized trajectories."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from branches import DELTA,E0,R,solve,equations,rotation,polar_ring


def loadx(row):
    return np.r_[np.array(row['nodes']).ravel(),row['period']/7,row['parameter_fraction']]


def stats(got,phase,level,p):
    x,J,mu,one,err=got;u=x[:-2].reshape(-1,3)
    return dict(e=float(E0+DELTA*x[-1]),parameter_fraction=float(x[-1]),period=float(7*x[-2]),multiplier=mu,
                neutral=one,residual=err,phase=phase,phase_level=level,nodes=u.tolist(),
                node_spread=float(np.max(np.linalg.norm(u-u.mean(0),axis=1))))


def arc_mid(left,right,phase,level,p,tol):
    normal=right-left;normal/=np.linalg.norm(normal);mid=(left+right)/2
    got=solve(mid,p,phase,level,normal,mid,tol)
    if got is None:raise RuntimeError('Pseudo-arclength midpoint failed')
    v=np.linalg.svd(got[1],full_matrices=True)[2][-1]
    if np.dot(v,normal)<0:v=-v
    return got,v


def event(rows,kind,p,tol):
    value=lambda row:row['de_ds'] if kind=='fold' else row['multiplier'][0]+1
    i=next(i for i in range(1,len(rows)) if value(rows[i-1])*value(rows[i])<0)
    phase=rows[0]['phase'];level=rows[0]['phase_level'];lo=loadx(rows[i-1]);hi=loadx(rows[i])
    sign=1 if value(rows[i-1])>0 else -1
    for _ in range(18):
        got,v=arc_mid(lo,hi,phase,level,p,tol)
        val=v[-1] if kind=='fold' else got[2][0]+1
        if val*sign>0:lo=got[0]
        else:hi=got[0]
    got,v=arc_mid(lo,hi,phase,level,p,tol)
    st=stats(got,phase,level,p);h=1e-5;other=[]
    for s in (-h,h):
        pred=got[0]+s*v
        a=solve(pred,p,phase,level,v,pred,tol);other.append(stats(a,phase,level,p))
    st.update(arclength_bracket=float(np.linalg.norm(hi-lo)),local_dmu_ds=(other[1]['multiplier'][0]-other[0]['multiplier'][0])/(2*h),
              local_d2e_ds2=(other[1]['e']-2*st['e']+other[0]['e'])/h**2)
    return st


def section_arc_at(rows,coordinate,wanted,p,tol,after=0):
    def value(r):return r[coordinate] if coordinate!='multiplier' else r[coordinate][0]
    i=next(i for i in range(max(1,after+1),len(rows)) if (value(rows[i-1])-wanted)*(value(rows[i])-wanted)<=0)
    lo=loadx(rows[i-1]);hi=loadx(rows[i]);phase=rows[0]['phase'];level=rows[0]['phase_level']
    if coordinate=='parameter_fraction':
        normal=np.zeros(len(lo));normal[-1]=1
        x=(lo+hi)/2;x[-1]=wanted;got=solve(x,p,phase,level,normal,x,tol)
        return stats(got,phase,level,p)
    sign=np.sign(value(rows[i-1])-wanted)
    for _ in range(16):
        got,v=arc_mid(lo,hi,phase,level,p,tol)
        if sign*(got[2][0]-wanted)>0:lo=got[0]
        else:hi=got[0]
    return stats(got,phase,level,p)


def reconstruct(row,p0,n,method='DOP853',tol=3e-12):
    p=p0.copy();p[8]=row['e'];nodes=np.asarray(row['nodes']);t=row['period']/7
    times=7*np.arange(n)/n;idx=np.floor(times).astype(int);z=np.empty((n,3));closures=[]
    if len(nodes)==1:nodes=np.repeat(nodes,7,axis=0)
    for j in range(7):
        sol=solve_ivp(lambda _,s:polar_ring(s,p),(0,t),nodes[j],method=method,rtol=tol,atol=tol/100,dense_output=True)
        assert sol.success
        closures.append(float(np.linalg.norm(sol.y[:,-1]-R@nodes[(j+1)%7])))
        z[idx==j]=sol.sol((times[idx==j]-j)*t).T@rotation(2*np.pi*4*j/7).T
    return z,closures


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ['old','branch','symmetric','out']:ap.add_argument('--'+key,type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    orig=json.loads((a.old/'final-link-crosscheck.json').read_text())['symmetric_companion'];p=np.array(orig['parameters'])
    rows=[json.loads(s) for s in a.branch.read_text().splitlines()];sym=[json.loads(s) for s in a.symmetric.read_text().splitlines()]
    result=dict(events={},representatives={});write=lambda:(a.out/'data.json').write_text(json.dumps(result,indent=2)+'\n')
    for tol in (8e-11,9e-13):
        marks={}
        for event_type in ('fold','flip'):
            marks[event_type]=event(rows,event_type,p,tol)
        result['events'][str(tol)]=marks;write()
    cases=dict(stable=section_arc_at(rows,'multiplier',-.45,p,2e-12,46),
               final_partner=section_arc_at(rows,'parameter_fraction',1.,p,2e-12,70))
    cases['unstable_at_stable_e']=section_arc_at(rows,'parameter_fraction',cases['stable']['parameter_fraction'],p,2e-12)
    cases['symmetric_at_stable_e']=section_arc_at(sym,'parameter_fraction',cases['stable']['parameter_fraction'],p,2e-12)
    for name,row in cases.items():
        n=32768;z,cl=reconstruct(row,p,n)
        zz,control=reconstruct(row,p,n,'RK45',2e-12)
        path=a.out/(name+'.npz');np.savez_compressed(path,points=z)
        row.update(case=180844,family='ring',identifier=name,parameters=[float(v) if i!=8 else row['e'] for i,v in enumerate(p)],
                   rotation_index=4,shape=list(z.shape),point_key='points',dtype='float64',path=str(path.resolve()),
                   cloud_relative=path.resolve().relative_to(a.out.resolve().parents[2]).as_posix(),
                   sha256=hashlib.sha256(path.read_bytes()).hexdigest(),solver_sector_closures=dict(DOP853=cl,RK45=control),
                   max_trajectory_discrepancy=float(np.max(np.linalg.norm(zz-z,axis=1))))
        result['representatives'][name]=row;write()
        print(name,row['e'],row['multiplier'],flush=True)


if __name__=='__main__':
    main()
