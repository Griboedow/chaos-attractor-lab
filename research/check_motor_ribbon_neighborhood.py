"""Independent spectra and geometry for a declared subset of the local survey."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
from pathlib import Path
import numpy as np
from check_motor_asymmetry import run,write


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);data=np.load(args.source)
    p=data['parameters'];complete=data['status']==0
    retained=complete&(data['lle']>.01)&(data['lle_windows'].min(axis=1)>.005)
    # First lane per retained neighbor; choose farthest-point coverage in the
    # three original logarithmic parameters. This selection is not based on images.
    candidates=np.array([i for i in range(32,len(p),2) if retained[i]],dtype=int)
    selected=[0];pool=candidates.tolist();logp=np.log(p)
    for _ in range(min(4,len(pool))):
        distances=np.min(np.linalg.norm(logp[pool,None,:]-logp[np.array(selected)[None,:],:],axis=2),axis=1)
        chosen=pool[int(np.argmax(distances))];selected.append(chosen);pool.remove(chosen)
    nonretained=np.array([i for i in range(32,len(p),2) if complete[i] and not retained[i]],dtype=int)
    if len(nonretained):
        order=nonretained[np.argsort(data['lle'][nonretained])]
        selected.extend(order[[0,-1]].tolist())
    selected=sorted(set(selected))
    jobs=[dict(id=i,parameters=p[i].tolist(),initial=data['final_states'][i].tolist(),warm=200.,duration=args.duration,
               output=str(out),rtol=2e-10,max_step=.02) for i in selected]
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        selection='Anchor lane0; four farthest-point retained first-lane neighbors in logarithmic parameter distance; lowest/highest-rate nonretained first-lane neighbors.',
        selected=selected,jobs=jobs,method='Independent DOP853 full8D QR.5tau;rtol2e-10,atol2e-12,maxstep.02,samples.02.',
        limitation='Finite numerical checks, no proof of an open basin or a common attractor.'))
    for name in ['check_motor_ribbon_neighborhood.py','check_motor_asymmetry.py','motor_rows.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            result=future.result();records.append(result)
            write(out/'analysis.json',dict(completed=len(records),declared=len(jobs),runs=sorted(records,key=lambda r:r['id'])))
            print(result['id'],result['spectrum'][:3],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--duration',type=float,default=3000.)
    main(p.parse_args())
