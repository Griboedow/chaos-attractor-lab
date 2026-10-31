"""Fresh starts and a declared local parameter box around the motor ribbon."""
import argparse
from pathlib import Path
import hashlib
import numpy as np
from scipy.stats import qmc
from motor_rows_gpu import MotorRowsGPU,save
from search_motor_asymmetry import write


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    reference=np.load(args.source);anchor=reference['parameters']
    d=reference['states'].shape[1];n=(d+1)//3
    assert d==3*n-1 and n>=3
    u=qmc.Sobol(3,scramble=True,seed=2026103145).random_base2(6)
    neighbors=anchor*np.exp(np.log(.95)+u*np.log(1.05/.95))
    params=np.vstack([np.tile(anchor,(32,1)),np.repeat(neighbors,2,axis=0)])
    rng=np.random.default_rng(2026103146);positions=rng.uniform(-.1,.1,(len(params),n));positions-=positions.mean(axis=1)[:,None]
    starts=np.column_stack([np.ones((len(params),n)),np.repeat(params[:,2,None],n,axis=1),positions[:,:n-1]])
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        anchor=anchor.tolist(),N=n,dimension=d,anchor_fresh_starts=32,neighbor_tuples=64,starts_per_neighbor=2,
        box='Each of rho,r,b lies independently in [.95,1.05] times the anchor; logarithmic scrambled Sobol.',
        parameter_seed=2026103145,start_seed=2026103146,
        initial='c_j=1,s_j=b; independent uniform[-.1,.1] displacement N-vectors centered to zero sum. No saved attractor endpoints.',
        warm=500.,observation=6000.,rtol=1e-9,step_factor=.5,time='original tau',
        limitation='A local finite numerical survey, not a rigorous open basin, parameter robustness theorem, or probability over all physical starts.'))
    np.savez_compressed(out/'declared-inputs.npz',parameters=params,initial=starts,positions=positions,
                        tuple_id=np.r_[np.zeros(32,dtype=int),np.repeat(np.arange(1,65),2)])
    for name in ['search_motor_ribbon_neighborhood.py','search_motor_asymmetry.py','motor_rows_gpu.py','motor_rows_gpu.cl','motor_rows.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    result=MotorRowsGPU(n).scan(params,starts,transient=500,duration=6000,rtol=1e-9,time_mode='original',
                               step_factor=.5,max_steps=15000000,max_wall=1200.)
    result['tuple_id']=np.r_[np.zeros(32,dtype=int),np.repeat(np.arange(1,65),2)]
    save(out/'complete-or-budgeted.npz',result)
    keep=(result['status']==0)&(result['lle']>.01)&(result['lle_windows'].min(axis=1)>.005)
    write(out/'local-analysis.json',dict(lanes=len(params),completed=int((result['status']==0).sum()),
        retained_anchor_ids=np.flatnonzero(keep[:32]).tolist(),retained_neighbor_lane_ids=(32+np.flatnonzero(keep[32:])).tolist(),
        retained_neighbor_tuple_ids=np.unique(result['tuple_id'][keep&(np.arange(len(params))>=32)]).tolist(),
        retention='mean leading rate>.01 and all four1500tau windows>.005; complete lanes only',
        warning='No attribution of the same attractor based on rates alone; histories and symmetry need separate comparison.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
