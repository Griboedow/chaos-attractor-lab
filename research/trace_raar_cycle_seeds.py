"""Declared finite evolution of tangent seeds near the certified unstable seven-cycle.

These are straight tangent segments, not a validated invariant-manifold chart.
Every radius/resolution outcome is retained; no attraction/chaos theorem follows.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.spatial import cKDTree
from continue_raar_cycles import iterate
from render_twist import camera


@njit(cache=True)
def value(x,p):
    beta,r,t,cx,cy,cz=p;sn=np.sin(t);co=np.cos(t)
    dx=x[0]-cx;dy=x[1]-cy;dz=x[2]-cz;along=co*dy+sn*dz;length=np.sqrt(dx*dx+along*along)
    if length<1e-12:return x,1,length
    bx=cx+r*dx/length;by=cy+r*co*along/length;bz=cz+r*sn*along/length
    rx=2*bx-x[0];ry=2*by-x[1];lengtha=np.sqrt(rx*rx+ry*ry)
    if lengtha<1e-12:return x,1,min(length,lengtha)
    a=np.array([rx/lengtha,ry/lengtha,0.]);b=np.array([bx,by,bz]);return beta*x+beta*a+(1-2*beta)*b,0,min(length,lengtha)


@njit(cache=True,parallel=True)
def advance(states,p,steps):
    result=states.copy();status=np.zeros(len(states),np.int64);minimum=np.full(len(states),np.inf)
    for i in prange(len(states)):
        s=states[i].copy()
        for k in range(steps):
            s,flag,axis=value(s,p);minimum[i]=min(minimum[i],axis)
            if flag or not np.isfinite(s).all():status[i]=flag or 2;break
        result[i]=s
    return result,status,minimum


@njit(cache=True,parallel=True)
def phases(states,p):
    result=np.empty((7,len(states),3));status=np.zeros(len(states),np.int64);minimum=np.full(len(states),np.inf)
    for i in prange(len(states)):
        s=states[i].copy()
        for k in range(7):
            result[k,i]=s
            if k<6:
                s,flag,axis=value(s,p);minimum[i]=min(minimum[i],axis)
                if flag or not np.isfinite(s).all():status[i]=flag or 2;break
    return result,status,minimum


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(8);begin=time.time();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    cert=json.loads(args.certificate.read_text());p=np.array(cert['parameters']);x=np.array([float(v) for v in cert['center_decimal']]);_,M,_,orbit,status=iterate(x,p,7);assert status==0
    eigen,vectors=np.linalg.eig(M);which=np.argmax(np.abs(eigen));assert abs(eigen[which].imag)<1e-12 and eigen[which].real < -1
    direction=vectors[:,which].real;direction/=np.linalg.norm(direction)
    if direction[np.argmax(abs(direction))]<0:direction=-direction
    with np.load(args.reference) as d:reference=d['states'].copy()
    tree=cKDTree(reference);floor=cKDTree(reference[:len(reference)//2]).query(reference[len(reference)//2:],workers=8)[0]
    render=json.loads(args.camera.read_text())['runs'][0];center=np.array(render['center']);basis=np.array(render['basis']);cam=camera(render['azimuth'],render['elevation']);projection=basis@cam.T
    snapshots_at=np.array([0,4,8,12,16,20,24,32,48,64]);radii=[1e-5,1e-7,1e-9];sizes=[4096,16384];half=3.4;resolution=(900,700);rows=[]
    protocol=dict(certificate=str(args.certificate),certificate_sha256=sha(args.certificate),reference=str(args.reference),reference_sha256=sha(args.reference),camera_source=str(args.camera),camera_sha256=sha(args.camera),params=p.tolist(),cycle_seed=x.tolist(),unstable_eigenvector=direction.tolist(),numerical_multiplier=float(eigen[which].real),radii=radii,seed_counts=sizes,returns=snapshots_at.tolist(),phase_offsets=list(range(7)),maximum_main_iterations=448,initial_rule='Evenly spaced straight segment x + t*v, t in [-radius,+radius], even count omits t=0.',selection='R1051 selected after the prior cycle study; all six declared radius/resolution combinations retained.',scope='Finite tangent-seed evolution, not an exact or certified unstable manifold, basin/attraction proof, independent discovery trial or historical-priority claim.',histogram=dict(resolution=list(resolution),x_bounds=[-half,half],y_bounds=[-half*resolution[1]/resolution[0],half*resolution[1]/resolution[0]],camera_center=center.tolist(),projection=projection.tolist(),all_65_return_windows_for_case='radius 1e-7, 16384 seeds'))
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    for name in ['trace_raar_cycle_seeds.py','continue_raar_cycles.py','search_raar_circles.py','render_twist.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (args.out/'cycle-certificate.json').write_bytes(args.certificate.read_bytes());map_errors=[]
    for y in reference[::10000]:map_errors.append(np.max(np.abs(value(y,p)[0]-iterate(y,p,1)[0])))
    assert max(map_errors)<1e-12
    for radius in radii:
        for count in sizes:
            start=time.time();states=x+np.linspace(-radius,radius,count)[:,None]*direction;assert len(np.unique(states,axis=0))==count
            snapshots=[];distances=[];bounds=[];minaxis=np.inf;failed=np.zeros(count,dtype=int);histograms=[];clipped=[]
            for k in range(65):
                phase,flag,axis=phases(states,p);failed=np.maximum(failed,flag);minaxis=min(minaxis,float(axis.min()));bounds.append([phase.reshape(-1,3).min(0),phase.reshape(-1,3).max(0)])
                if k in snapshots_at:
                    snapshots.append(phase.copy());distances.append(tree.query(phase.reshape(-1,3),workers=8)[0])
                if radius==1e-7 and count==16384:
                    xy=(phase.reshape(-1,3)-center)@projection
                    h,_,_=np.histogram2d(xy[:,0],xy[:,1],bins=resolution,range=[[-half,half],[-half*resolution[1]/resolution[0],half*resolution[1]/resolution[0]]]);histograms.append(h.T[::-1].astype(np.uint32));clipped.append(int(len(xy)-h.sum()))
                if k<64:
                    states,flag,axis=advance(states,p,7);failed=np.maximum(failed,flag);minaxis=min(minaxis,float(axis.min()))
            key=f'radius-{radius:.0e}-n{count}';np.savez_compressed(args.out/(key+'.npz'),params=p,seed=x,direction=direction,returns=snapshots_at,snapshots=np.array(snapshots),reference_distances=np.array(distances),reference_split_distances=floor,bounds=np.array(bounds),status=failed)
            if histograms:
                np.savez_compressed(args.out/'evolution-histograms.npz',histograms=np.array(histograms),returns=np.arange(65),clipped=np.array(clipped));assert max(clipped)==0
            row=dict(key=key,radius=radius,seed_count=count,completed=int(np.sum(failed==0)),minimum_axis_distance=minaxis,maximum_state_norm_bound=float(np.max(np.linalg.norm(np.array(bounds),axis=2))),distance_quantiles=np.quantile(np.array(distances),[.5,.9,.99,1],axis=1).T.tolist(),reference_split_quantiles=np.quantile(floor,[.5,.9,.99,1]).tolist(),seconds=time.time()-start)
            rows.append(row);print(json.dumps(row),flush=True)
    report=dict(protocol,runs=rows,max_scalar_map_residual=float(max(map_errors)),seconds=time.time()-begin,source_sha256=sha(Path(__file__)))
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,default=Path('results/raar-cycle-structure/cycle-proof/base-unstable-period-7.json'));p.add_argument('--reference',type=Path,default=Path('results/raar-circles/render/case-1051-animation.npz'));p.add_argument('--camera',type=Path,default=Path('results/raar-circles/render/analysis.json'));p.add_argument('--out',type=Path,required=True);main(p.parse_args())
