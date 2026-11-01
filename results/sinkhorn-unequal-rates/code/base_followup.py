"""Continue every positive pilot flag; retain complete final windows explicitly."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads
import model

@njit(cache=True, parallel=True)
def follow(matrices, starts, omegas, block, blocks, keep):
    count,d=starts.shape
    stats=np.zeros((count,blocks,4+2*d))
    checkpoints=np.zeros((count,blocks,d))
    tails=np.empty((count,keep,d))
    for lane in prange(count):
        L=matrices[lane]; la,lb=model.targets(L); w=omegas[lane]
        z=starts[lane].copy();v=np.ones(d)/np.sqrt(d)
        for b in range(blocks):
            rate=0.;s=np.zeros(d);ss=np.zeros(d);maximum=0.;minimum=1e300;residual=0.
            for k in range(block):
                znew,v=model.step(z,v,L,la,lb,w)
                residual=np.max(np.abs(znew-z));z=znew
                norm=np.linalg.norm(v);rate+=np.log(max(norm,1e-300))
                v/=max(norm,1e-300)
                s+=z;ss+=z*z
                maximum=max(maximum,np.max(np.abs(z)));minimum=min(minimum,np.max(np.abs(z)))
                index=b*block+k-(blocks*block-keep)
                if index>=0:tails[lane,index]=z
            stats[lane,b,0]=rate/block
            stats[lane,b,1]=minimum;stats[lane,b,2]=maximum;stats[lane,b,3]=residual
            stats[lane,b,4:4+d]=s/block
            stats[lane,b,4+d:]=np.sqrt(np.maximum(0,ss/block-(s/block)**2))
            checkpoints[lane,b]=z
    return stats,checkpoints,tails

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pilot',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--block',type=int,default=25000);ap.add_argument('--blocks',type=int,default=20)
    ap.add_argument('--keep',type=int,default=25000)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    for n in (3,4):
        p=json.loads((args.pilot/f'n{n}.json').read_text())
        rows=[r for r in p['runs'] if r['capture_step']==0 and r['finite_rate']>1e-5]
        endpoints=np.load(args.pilot/f'n{n}-endpoints-and-last512.npz')['endpoints']
        L=np.array([r['log_kernel'] for r in rows]);z=np.array([endpoints[r['case']] for r in rows])
        w=np.array([r['omega'] for r in rows]);t=time.perf_counter()
        stats,checkpoints,tails=follow(L,z,w,args.block,args.blocks,args.keep)
        output=[]
        for i,r in enumerate(rows):
            output.append(dict(case=r['case'],shape=[n,n],omega=r['omega'],start_step=30000,
                starting_state=z[i].tolist(),log_kernel=r['log_kernel'],
                blocks=[dict(stop_step=30000+(j+1)*args.block,rate=float(s[0]),min_abs_state=float(s[1]),max_abs_state=float(s[2]),
                    final_step_residual=float(s[3]),mean=s[4:4+2*n-2].tolist(),std=s[4+2*n-2:].tolist(),checkpoint=checkpoints[i,j].tolist()) for j,s in enumerate(stats[i])]))
        record=dict(block_length=args.block,blocks=args.blocks,retained_final_window_steps=args.keep,
            elapsed_seconds=time.perf_counter()-t,runs=output,
            sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(model.__file__))})
        (args.output/f'n{n}.json').write_text(json.dumps(record,indent=2)+'\n')
        np.savez_compressed(args.output/f'n{n}-final-{args.keep}-steps.npz',case=np.array([r['case'] for r in rows]),states=tails)
        print(json.dumps(dict(size=n,elapsed=record['elapsed_seconds'],cases=[dict(case=r['case'],first=r['blocks'][0]['rate'],last=r['blocks'][-1]['rate'],norm=r['blocks'][-1]['max_abs_state']) for r in output])),flush=True)
