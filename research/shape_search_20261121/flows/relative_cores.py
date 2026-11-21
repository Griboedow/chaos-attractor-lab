"""Use C7-equivariant relative shooting to continue unstable closed curves over one sector."""
import argparse
import hashlib
import json
import os
from pathlib import Path

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from adaptive import full_rhs,polar_ring
from orbits import primary_rhs


def rotation(angle):
    c=np.cos(angle);s=np.sin(angle)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])


def solve_sector(seed,t,p,R):
    q=seed.copy();i=int(np.argmax(abs(polar_ring(q,p))));free=[j for j in range(3) if j!=i]
    for j in range(18):
        sol=solve_ivp(primary_rhs,(0,t),np.r_[q,np.eye(3).ravel()],args=(p,1),method='DOP853',rtol=8e-13,atol=8e-15)
        if not sol.success:return
        final=sol.y[:,-1];err=final[:3]-R@q
        if np.linalg.norm(err)<1e-11:return q,t
        m=final[3:].reshape(3,3)-R
        delta=np.linalg.solve(np.column_stack((m[:,free],polar_ring(final[:3],p))),-err)
        scale=max(1.,max(abs(delta[:2]))/.05,abs(delta[2])/.25)
        q[free]+=delta[:2]/scale;t+=delta[2]/scale
    return None


def sector_cloud(seed,t,p,k,n):
    R=rotation(2*np.pi*k/7)
    sol=solve_ivp(full_rhs,(0,t),np.r_[seed,np.eye(3).ravel(),0.],args=(p,1),
                  method='DOP853',rtol=3e-12,atol=3e-14,dense_output=True)
    end=sol.y[:,-1];μ=np.linalg.eigvals(R.T@end[3:12].reshape(3,3))
    li=np.argmin(abs(μ-1));other=np.delete(μ,li);mu=max(other,key=abs)
    times=7*np.arange(n)/n;sector=np.floor(times).astype(int);base=sol.sol((times-sector)*t)[:3].T
    pts=np.empty_like(base)
    for i in range(7):pts[sector==i]=base[sector==i]@rotation(2*np.pi*k*i/7).T
    check=dict(relative_rotation_index=k,sector_period=t,sector_return_error_into_rotated_center=float(np.linalg.norm(end[:3]-R@seed)),
               sector_neutral=[float(μ[li].real),float(μ[li].imag)],sector_unstable=[float(mu.real),float(mu.imag)],
               neutral_multiplier=[float((μ[li]**7).real),float((μ[li]**7).imag)],
               dominant_multiplier=[float((mu**7).real),float((mu**7).imag)],
               sector_integrated_divergence=float(end[12]))
    return pts,check


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,required=True);ap.add_argument('--branch',type=Path,required=True)
    ap.add_argument('--long-return',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    m=json.loads((a.packet/'manifest.json').read_text());anchor=next(r for r in m['cases'] if r['case']==180844)
    p=np.asarray(anchor['parameters']).copy();p0=next(r['parameters'][8] for r in m['stable_cycles'] if r['case']==28756);p1=p[8]
    cloud=np.load(anchor['path'])['points'];scale=np.linalg.norm(np.std(cloud,axis=0));tree=cKDTree(cloud)
    table={r['fraction']:r for r in json.loads(a.branch.read_text())['trials']}
    for order,start,k in [(4,.94,1),(8,.945,2),(16,.95,4)]:
        row=table[start].get('attracting_cycle') or json.loads(a.long_return.read_text())['accepted_cycles'][0]
        q=np.array(row['refined_start']);t=row['period']/7;R=rotation(2*np.pi*k/7);f=start;h=.002
        history=[]
        while f<1-1e-14:
            f1=min(1.,f+h);p[8]=(1-f1)*p0+f1*p1
            cand=solve_sector(q,t,p,R)
            if cand is None:
                h*=.5
                if h<2e-6:break
                continue
            q,t=cand;f=f1;history.append([f,t]);h=min(.003,h*1.25)
        record=dict(period_order=order,achieved_fraction=f,completed=False,history=history)
        if f>=1-1e-14:
            p[8]=p1;n=65536 if order>=8 else 32768
            pts,check=sector_cloud(q,t,p,k,n)
            u=check['dominant_multiplier'][0]
            valid=check['sector_return_error_into_rotated_center']<2e-7 and abs(check['neutral_multiplier'][0]-1)<.002 and abs(u)>1.01
            record.update(**check,completed=valid)
            if valid:
                path=a.out/f'orbit-{order}.npz';np.savez_compressed(path,points=pts)
                record.update(case=180844,family='ring',label=f'Unstable C7-relative orbit of order {order}',parameters=p.tolist(),
                    initial=q.tolist(),period=7*t,shape=[n,3],point_dtype='float64',point_key='points',sample_time_step=7*t/n,
                    method='C7 relative shooting',path=str(path.resolve()),cloud_relative=path.resolve().relative_to(a.out.resolve().parents[1]).as_posix(),
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest(),orbit_to_cloud_p95_rms=float(np.quantile(tree.query(pts[::5])[0],.95)/scale),
                    cloud_to_orbit_p95_rms=float(np.quantile(cKDTree(pts).query(cloud[::16])[0],.95)/scale))
        (a.out/f'order-{order}.json').write_text(json.dumps(record,indent=2)+'\n')
        print(order,f,record['completed'],record.get('dominant_multiplier'),flush=True)


if __name__=='__main__':
    main()
