"""Finite volumes approaching the selected accepted-step clouds."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc
from scipy.spatial import cKDTree
from search_rk45_oscillator import matrices,oscillator_step

@njit(cache=True,parallel=True)
def evolve(params,starts,times):
    ms=matrices(params[0]);n=len(starts);snapshots=np.full((len(times),n,3),np.nan);status=np.zeros(n,dtype=np.int64)
    for i in prange(n):
        state=starts[i].copy();j=0
        for t in range(times[-1]+1):
            if t==times[j]:
                snapshots[j,i]=np.array([state[0],state[1],np.exp(state[2])]);j+=1
                if j==len(times):break
            state,_,_,code=oscillator_step(state,params[1],*ms)
            if code or not np.all(np.isfinite(state)) or abs(state[2])>50 or max(abs(state[0]),abs(state[1]))>1e10:
                status[i]=code+1;break
    return snapshots,status

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    with np.load(args.source) as f:params=f['params']
    times=np.array([0,4,16,64,256,1024],dtype=np.int64);records=[]
    for ident in [335,448]:
        raw=qmc.Sobol(3,scramble=True,seed=2026110220+ident).random_base2(15)
        starts=np.column_stack([-4+8*raw[:,0],-4+8*raw[:,1],np.log(.25+7.75*raw[:,2])])
        begin=time.time();snapshots,status=evolve(params[ident],starts,times)
        with np.load(args.gallery/f'case-{ident}.npz') as f:reference=f['states']
        tree=cKDTree(reference);distances=np.full((len(times),len(starts)),np.nan)
        for j in range(len(times)):
            valid=np.all(np.isfinite(snapshots[j]),axis=1);distances[j,valid]=tree.query(snapshots[j,valid],workers=4)[0]
        # Reference split establishes the finite sampling distance floor.
        baseline=cKDTree(reference[:100000]).query(reference[100000:],workers=4)[0]
        np.savez_compressed(args.out/f'case-{ident}.npz',params=params[ident],starts=starts,times=times,snapshots=snapshots,status=status,distances=distances,reference_split_distances=baseline)
        row=dict(id=ident,params=params[ident].tolist(),count=len(starts),times=times.tolist(),initial_box=[[-4.,4.],[-4.,4.],[.25,8.]],initial_distribution='Scrambled Sobol uniform in actual (q/atol,p/atol,h), converted to log h only for computation',
                 completed=int(np.sum(status==0)),distance_quantiles={str(t):np.nanquantile(d,[.1,.5,.9,.99]).tolist() for t,d in zip(times,distances)},reference_split_distance_quantiles=np.quantile(baseline,[.1,.5,.9,.99]).tolist(),reference_points=len(reference),seconds=time.time()-begin,
                 scope='Distance to a finite sampled reference cloud in equal original-state units; not a certified basin measure or invariant-set distance.')
        records.append(row);print(json.dumps(row),flush=True)
    report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),runs=records)
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for n in ['approach_rk45_oscillator.py','search_rk45_oscillator.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--gallery',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
