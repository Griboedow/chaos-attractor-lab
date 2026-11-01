"""Follow a declared unit-ball ensemble toward the observed twist-map regimes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.spatial import cKDTree
from scipy.stats import qmc
from twist_map import step


@njit(parallel=True)
def evolve(p,initial,times):
    count=len(initial);snapshots=np.empty((len(times),count,3));tails=np.empty((count,24,3))
    for lane in prange(count):
        z=initial[lane].copy();next_time=0
        for t in range(times[-1]+1):
            if t==times[next_time]:
                snapshots[next_time,lane]=z
                if next_time<len(times)-1:next_time+=1
            if t>times[-1]-24:tails[lane,t-(times[-1]-23)]=z
            if t<times[-1]:z=step(z,p)
    return snapshots,tails


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    seed=2026110118;u=qmc.Sobol(3,scramble=True,seed=seed).random_base2(15)
    rad=u[:,0]**(1/3);cost=2*u[:,1]-1;angle=2*np.pi*u[:,2]
    initial=np.c_[rad*np.sqrt(1-cost**2)*np.cos(angle),rad*np.sqrt(1-cost**2)*np.sin(angle),rad*cost]
    times=np.array([0,1,2,4,8,16,32,64,256,1024,4096],dtype=np.int64)
    long=np.load(args.source/'long-runs.npz');ids=[1709,3254]
    protocol=dict(source=str(args.source),sha256=hashlib.sha256((args.source/'long-runs.npz').read_bytes()).hexdigest(),
        seed=seed,initial='32768 scrambled Sobol points transformed to uniform unit-ball volume using spherical coordinates.',
        candidates=ids,snapshot_times=times.tolist(),endpoint_test='Maximum full-state difference at lag3 across final24states<1e-8. Also test lag1 and lag2.',
        support_distance='Euclidean nearest-reference-tail distance; retain quantiles for every snapshot. It is not distance to the exact invariant set.',
        scope='Ensemble attraction and basin sample visualization; no Lyapunov calculation or asymptotic classification for these new starts. All declared points are retained.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_twist_approach.py','twist_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    for ident in ids:
        j=list(long['ids']).index(ident);p=long['parameters'][j];reference=long['tail'][j];tree=cKDTree(reference)
        snapshots,tails=evolve(p,initial,times);assert np.max(np.linalg.norm(snapshots,axis=2))<=1+1e-12
        residual=np.array([np.max(np.linalg.norm(tails[:,lag:]-tails[:,:-lag],axis=2),axis=1) for lag in [1,2,3]]).T
        distance=np.array([tree.query(x,workers=4)[0] for x in snapshots])
        path=out/f'case-{ident}.npz';np.savez_compressed(path,initial=initial,snapshot_times=times,snapshots=snapshots,final24= tails,
            nearest_reference_distance=distance,return_residual=residual,parameters=p)
        row=dict(id=ident,file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            period1_count=int(np.count_nonzero(residual[:,0]<1e-8)),period2_count=int(np.count_nonzero(residual[:,1]<1e-8)),period3_count=int(np.count_nonzero(residual[:,2]<1e-8)),
            support_distance_quantiles=[0,.5,.95,.99,1],distance_by_snapshot=np.quantile(distance,[0,.5,.95,.99,1],axis=1).T.tolist())
        records.append(row);(out/'analysis.json').write_text(json.dumps(dict(runs=records),indent=2)+'\n')
        print(ident,'period3',row['period3_count'],'final distance',row['distance_by_snapshot'][-1],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
