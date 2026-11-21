"""Follow post-flip attractors from the stable seven-copy interval using independent shooting seeds."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from branches import polar_ring
from field import integrate,ring_tangent
from orbits import refine
from adaptive import full_rhs


def guesses(z,dt,base,rms):
    out=[]
    for axis in range(3):
        level=np.median(z[:,axis]);idx=np.flatnonzero((z[:-1,axis]<level)&(z[1:,axis]>=level))
        a=(level-z[idx,axis])/(z[idx+1,axis]-z[idx,axis]);u=z[idx]+a[:,None]*(z[idx+1]-z[idx]);t=(idx+a)*dt
        for lag in range(1,min(len(t)-20,400)):
            period=np.mean(t[-20:]-t[-20-lag:-lag])
            if period>8.1*base:break
            mult=round(period/base)
            if mult not in [1,2,4,8] or abs(period/base-mult)>.015:continue
            error=np.max(np.linalg.norm(u[-20:]-u[-20-lag:-lag],axis=1))/rms
            if error<.0007:out.append((mult,error,axis,u[-1-lag],period))
    return sorted(out,key=lambda z:(z[0],z[1]))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--data',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    q=json.loads(a.data.read_text());stable=q['representatives']['stable'];ef=list(q['events'].values())[-1]['flip']['e']
    p=np.array(stable['parameters']);v=np.array(stable['nodes'][0]);base=stable['period'];z0=v+np.array([.00003,-.000011,.000023])
    report=[]
    for offset in [1e-7,3e-7,7e-7,1.2e-6,2e-6,3e-6,4e-6,6e-6,8e-6,1e-5]:
        p[8]=ef+offset
        h=.008;rates,z=integrate(p,z0,h,round(250*base/h),round(5*base/h),8,4,ring_tangent)
        found=guesses(z,.032,base,.89);row=dict(e=float(p[8]),offset_from_flip=offset,rates=rates.tolist(),guesses=len(found))
        for mult,err,axis,guess,T in found[:6]:
            got=refine(axis,guess,T,.89,p,1)
            if got is None:continue
            loc,per=got;sol=solve_ivp(full_rhs,(0,per),np.r_[loc,np.eye(3).ravel(),0],args=(p,1),method='DOP853',
                                     rtol=2e-12,atol=2e-14,dense_output=True)
            end=sol.y[:,-1];closure=np.linalg.norm(end[:3]-loc);vals=np.linalg.eigvals(end[3:12].reshape(3,3));k=np.argmin(abs(vals-1))
            mu=max(np.delete(vals,k),key=abs)
            if closure>2e-6 or abs(vals[k]-1)>.002 or abs(mu)>.995:continue
            tt=np.linspace(0,per,19,endpoint=False)
            floor=min(np.max(np.linalg.norm(sol.sol(tt)[:3]-sol.sol((tt+per/k)%per)[:3],axis=0))/.89 for k in range(2,65))
            if floor<.00015:continue
            row.update(order_relative_to_2T=mult,period=per,original_section_error=err,initial=loc.tolist(),
                       multiplier=float(mu.real),neutral=float(vals[k].real),independent_closure=float(closure),primitive_floor=floor)
            if mult>=2 and not any(x.get('cloud') and x.get('order_relative_to_2T')==mult for x in report):
                import hashlib
                pts=sol.sol(np.linspace(0,per,65536,endpoint=False))[:3].T
                path=a.out/f'daughter-{mult}.npz';np.savez_compressed(path,points=pts)
                row.update(cloud=path.name,cloud_relative=path.resolve().relative_to(a.out.resolve().parents[2]).as_posix(),
                           path=str(path.resolve()),shape=list(pts.shape),point_key='points',parameters=p.tolist(),
                           sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            break
        report.append(row);(a.out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(offset,row.get('order_relative_to_2T'),row.get('multiplier'),flush=True)


if __name__=='__main__':
    main()
