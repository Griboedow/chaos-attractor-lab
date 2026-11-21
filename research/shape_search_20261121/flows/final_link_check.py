"""Compare the asymmetric 2T orbit and its C7-relative companion against the same 4T core."""
import argparse
import hashlib
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from adaptive import polar_ring
from dyadic_links import coordinates
from exact_knot import prepare
from relative_cores import rotation,solve_sector,sector_cloud
from symmetry_links import link_integer


def solve_state(x,t,k,p,method,tol):
    sol=solve_ivp(lambda _,z:polar_ring(z,p),(0,t),x,method=method,rtol=tol,atol=tol/100,dense_output=True)
    assert sol.success
    def at(time):
        phases=np.asarray(time)/t
        sector=np.floor(phases).astype(int)
        loc=sol.sol((phases-sector)*t).T
        if np.ndim(phases)==0:return rotation(2*np.pi*k*int(sector)/7)@loc
        out=np.empty_like(loc)
        for j in range(7):
            use=sector%7==j
            out[use]=loc[use]@rotation(2*np.pi*k*j/7).T
        return out
    return at,float(np.linalg.norm(sol.y[:,-1]-rotation(2*np.pi*k/7)@x))


def nearest(f,g,T,S,p):
    # Optimize smooth inter-curve separation from the best sampled neighborhoods,
    # restricted to one symmetry sector of f when T is a sector period.
    times=np.linspace(0,T,8192,endpoint=False);other=np.linspace(0,S,57344,endpoint=False)
    a=f(times);b=g(other);dist,ix=cKDTree(b).query(a)
    ranked=np.argsort(dist);seeds=[]
    for i in ranked:
        if all(min(abs(i-j),len(times)-abs(i-j))>len(times)/60 for j in seeds):seeds.append(i)
        if len(seeds)>=24:break
    answers=[]
    for i in seeds:
        def fun(q):return f(q[0])-g(q[1])
        def jac(q):return np.column_stack((polar_ring(f(q[0]),p),-polar_ring(g(q[1]),p)))
        opt=least_squares(fun,[times[i],other[ix[i]]],jac=jac,ftol=1e-14,xtol=1e-14,gtol=1e-14,max_nfev=200)
        answers.append(dict(native_distance=float(np.linalg.norm(opt.fun)),times=opt.x.tolist()))
    return min(answers,key=lambda s:s['native_distance'])


def minimum_segments(a,b):
    """Evaluate every segment pair that can improve the current vertex-distance bound."""
    au=np.roll(a,-1,axis=0)-a;bv=np.roll(b,-1,axis=0)-b
    al=np.linalg.norm(au,axis=1);bl=np.linalg.norm(bv,axis=1)
    mid_a=a+au/2;mid_b=b+bv/2;tree=cKDTree(mid_b)
    best=float(cKDTree(b).query(a)[0].min());tested=0
    def point_edge(x,v):
        t=np.clip(np.einsum('ij,ij->i',x,v)/np.einsum('ij,ij->i',v,v),0,1)
        d=x-t[:,None]*v
        return np.einsum('ij,ij->i',d,d)
    for i in range(0,len(a),1024):
        sl=slice(i,min(len(a),i+1024))
        near=tree.query_ball_point(mid_a[sl],best+(al[sl]+bl.max())/2+1e-11)
        counts=np.array([len(n) for n in near])
        if not counts.sum():continue
        ii=np.repeat(np.arange(i,i+len(near)),counts);jj=np.concatenate([q for q in near if q])
        u=au[ii];v=bv[jj];w=b[jj]-a[ii]
        d=np.minimum.reduce([point_edge(-w,v),point_edge(u-w,v),point_edge(w,u),point_edge(w+v,u)])
        n=np.cross(u,v);den=np.einsum('ij,ij->i',n,n)
        s=np.einsum('ij,ij->i',np.cross(w,v),n)/np.maximum(den,1e-300)
        t=np.einsum('ij,ij->i',np.cross(w,u),n)/np.maximum(den,1e-300)
        inside=(den>0)&(s>=0)&(s<=1)&(t>=0)&(t<=1)
        z=s[inside,None]*u[inside]-t[inside,None]*v[inside]-w[inside]
        d[inside]=np.minimum(d[inside],np.einsum('ij,ij->i',z,z))
        best=min(best,float(np.sqrt(d.min())));tested+=len(d)
    return dict(vertices_each=len(a),minimum_segment_distance_native=best,candidate_segment_pairs=tested)


def distance_audit(packet,report):
    rows={r['period_order']:r for r in json.loads((packet/'nested-skeletons.json').read_text())['cores']}
    two=report['symmetric_companion'];four=rows[4];p=np.array(two['parameters']);t=two['sector_period'];u=four['sector_period']
    F,_=solve_state(np.array(two['initial']),t,4,p,'DOP853',4e-13)
    G,_=solve_state(np.array(four['initial']),u,1,p,'DOP853',4e-13)
    old=rows[2];T=old['period']
    sol=solve_ivp(lambda _,z:polar_ring(z,p),(0,T),np.array(old['initial']),method='DOP853',rtol=4e-13,atol=4e-15,dense_output=True)
    checks=[]
    for n in [14336,28672,57344]:
        phases=np.arange(n)/n;z2=F(phases*7*t);z4=G(phases*7*u);old2=sol.sol(phases*T).T
        checks.append(dict(vertices=n,symmetric=minimum_segments(z2,z4),asymmetric=minimum_segments(old2,z4)))
    report['polygon_distance_method']='Full polygon meshes; conservative midpoint-ball exclusion followed by float64 endpoint-to-segment and interior line-to-line tests. These distances concern each polygon, without a rigorous smooth-curve error bound.'
    report['polygon_distance_checks']=checks


def symmetry_audit(packet,report):
    from scipy.optimize import brentq
    row=next(x for x in json.loads((packet/'nested-skeletons.json').read_text())['cores'] if x['period_order']==2)
    p=np.array(row['parameters']);v=np.array(row['initial']);T=row['period'];rots=[rotation(2*np.pi*k/7) for k in range(7)]
    results=[]
    for method,tol in [('DOP853',4e-13),('RK45',2e-12)]:
        sol=solve_ivp(lambda _,z:polar_ring(z,p),(0,T),v,method=method,rtol=tol,atol=tol/100,dense_output=True)
        ts=np.linspace(0,T,16385);xyz=sol.sol(ts).T;states=sol.sol(np.arange(1,7)*T/7).T
        phase=[[float(np.linalg.norm(z-R@v)) for R in rots[1:]] for z in states]
        cuts=[]
        for k,R in enumerate(rots):
            level=(xyz@R.T)[:,1];hits=np.flatnonzero((level[:-1]<0)&(level[1:]>=0));pts=[]
            for j in hits:
                root=brentq(lambda tt:(R@sol.sol(tt))[1],ts[j],ts[j+1],xtol=4e-15)
                pt=R@sol.sol(root);pts.append([float(pt[0]),float(pt[2])])
            cuts.append(np.array(pts))
        separate=[]
        for i,x in enumerate(cuts):
            for j,y in enumerate(cuts[:i]):
                sep=np.linalg.norm(x[:,None,:]-y[None,:,:],axis=2)
                # Setwise equality would force both directions of the Hausdorff comparison to vanish.
                h=max(float(sep.min(0).max()),float(sep.min(1).max()))
                separate.append(dict(components=[j,i],closest_section_points=float(sep.min()),section_Hausdorff=h))
        results.append(dict(method=method,rtol=tol,positive_y_zero_section=[z.tolist() for z in cuts],
                           all_nontrivial_seventh_period_to_rotations=phase,
                           minimum_nontrivial_phase_rotation_mismatch=float(np.min(phase)),
                           minimum_pairwise_section_distance=min(z['closest_section_points'] for z in separate),
                           minimum_setwise_section_Hausdorff=min(z['section_Hausdorff'] for z in separate),
                           pairwise_section_distances=separate))
    report['asymmetric_orbit_group']=dict(group='C7',expected_orbit_multiplicity=7,
        argument='A nontrivial stabilizer would equal C7 because seven is prime, forcing one of the six nonzero seventh-period shifts to equal the corresponding rotation of the initial state. Shared positive y=0 sections independently test the seven component sets.',
        numerical_periodic_orbits_at_final_parameters_including_all_copies=12,solvers=results)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--refine-distance',action='store_true',help='Add global polygon distance controls to an existing report')
    ap.add_argument('--symmetry',action='store_true',help='Check all rotated non-symmetric curves on a common section')
    ap.add_argument('--replay',type=Path,help='Rebuild the registered companion NPZ at this fresh path')
    a=ap.parse_args();file=a.packet/'final-link-crosscheck.json'
    if a.symmetry:
        record=json.loads(file.read_text());symmetry_audit(a.packet,record)
        file.write_text(json.dumps(record,indent=2)+'\n');return
    if a.refine_distance:
        record=json.loads(file.read_text());distance_audit(a.packet,record)
        file.write_text(json.dumps(record,indent=2)+'\n');return
    if a.replay:
        record=json.loads(file.read_text())['symmetric_companion'];z,_=sector_cloud(np.array(record['initial']),
             record['sector_period'],np.array(record['parameters']),4,record['shape'][0])
        np.savez_compressed(a.replay,points=z)
        assert hashlib.sha256(a.replay.read_bytes()).hexdigest()==record['sha256'];return
    a.out.mkdir(exist_ok=False,parents=True)
    rows={r['period_order']:r for r in json.loads((a.packet/'nested-skeletons.json').read_text())['cores']}
    p=np.array(rows[2]['parameters']);old=np.array(rows[2]['initial']);oldT=rows[2]['period']
    R=rotation(2*np.pi*4/7)
    q,t=solve_sector(old,oldT/7,p,R)
    candidate,stats=sector_cloud(q,t,p,4,65536)
    target=a.out/'symmetric-2T.npz';np.savez_compressed(target,points=candidate)
    m=json.loads((a.packet/'manifest.json').read_text());parent=next(r for r in m['cases'] if r['case']==180844)
    cloud=np.load(a.cache/parent['cloud_relative'])['points'];scale=np.linalg.norm(np.std(cloud,axis=0))
    from periodic_catalog import shifts
    gap,divisor=shifts(candidate)
    record=dict(case=180844,period_order=2,label='C7-symmetric unstable order-two companion',family='ring',period=7*t,
       **stats,initial=q.tolist(),parameters=p.tolist(),point_key='points',point_dtype='float64',shape=list(candidate.shape),
       path=str(target.resolve()),cloud_relative=target.resolve().relative_to(a.cache.resolve()).as_posix(),
       sha256=hashlib.sha256(target.read_bytes()).hexdigest(),minimum_phase_shift_rms=gap,worst_divisor=divisor,
       cloud_to_orbit_p95_rms=float(np.quantile(cKDTree(candidate).query(cloud[::16])[0],.95)/scale),
       orbit_to_cloud_p95_rms=float(np.quantile(cKDTree(cloud).query(candidate[::5])[0],.95)/scale))
    report=dict(case=180844,symmetric_companion=record,asymmetric_period=oldT,asymmetric_initial=old.tolist(),
                original_2T_full_time_polygons_linking=-90,
                method='Relative shooting for Phi(t,q)=R(8pi/7) q; independent state-only adaptive integration. All seven sectors are assembled before exact dyadic linking. Minimum distances are local smooth refinements from 24 separated nearest-neighbor seeds, not certified global bounds.',
                reproduction='python research/shape_search_20261121/flows/final_link_check.py --packet results/shape_search_20261121/flows --cache <campaign-cache> --out <scratch-dir> --replay <fresh.npz>',
                checks=[])
    write=lambda:file.write_text(json.dumps(report,indent=2)+'\n')
    write()
    four=np.array(rows[4]['initial']);u=rows[4]['sector_period']
    for method,tol,nlist in [('DOP853',4e-13,[2048,8192,16384]),('RK45',2e-12,[8192])]:
        F,rf=solve_state(q,t,4,p,method,tol);G,rg=solve_state(four,u,1,p,method,tol)
        orig=solve_ivp(lambda _,z:polar_ring(z,p),(0,oldT),old,method=method,rtol=tol,atol=tol/100,dense_output=True)
        def old_at(when):return orig.sol(np.asarray(when)%oldT).T
        chk=dict(method=method,rtol=tol,symmetric_two_relative_closure=rf,four_relative_closure=rg,
                 original_full_return=float(np.linalg.norm(orig.y[:,-1]-old)),original_2T_symmetry_defect=float(np.linalg.norm(orig.sol(oldT/7)-R@old)),
                 minimum_2T_symmetric_to_4T=nearest(F,G,t,7*u,p),minimum_2T_asymmetric_to_4T=nearest(old_at,G,oldT,7*u,p),
                 period_difference_native=7*t-oldT,meshes=[])
        report['checks'].append(chk);write()
        for per_sector in nlist:
            n=7*per_sector
            z1=F(np.arange(n)*t/per_sector);z2=G(np.arange(n)*u/per_sector)
            ints,den=coordinates([z1,z2]);views=[]
            for seed in (17,35):
                _,cam,_,_=prepare(z1[::64],seed);lk,cross=link_integer(*ints,cam)
                views.append(dict(camera=cam.tolist(),linking=lk,crossings=cross))
            chk['meshes'].append(dict(vertices=n,exact_common_dyadic_denominator=str(den),views=views));write()
            print(method,n,views[0]['linking'],flush=True)


if __name__=='__main__':
    main()
