"""Continue every positive unequal-rate flag, retaining escape outcomes."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numba import njit,prange,set_num_threads
import model,model_dual

@njit(cache=True,parallel=True)
def follow(Ls,starts,ws,block,blocks,keep):
    count,d=starts.shape;stats=np.full((count,blocks,4+2*d),np.nan);checkpoints=np.full((count,blocks,d),np.nan);tails=np.full((count,keep,d),np.nan);completed=np.zeros(count,np.int64)
    for lane in prange(count):
        L=Ls[lane];la,lb=model.targets(L);z=starts[lane].copy();v=np.ones(d)/np.sqrt(d);wx,wy=ws[lane];escaped=False
        for b in range(blocks):
            total=0.;s=np.zeros(d);ss=np.zeros(d);maximum=0.;minimum=1e300;residual=0.
            for k in range(block):
                zn,v=model_dual.step(z,v,L,la,lb,wx,wy);residual=np.max(np.abs(zn-z));z=zn;norm=np.linalg.norm(v)
                if norm>1e-250:v/=norm
                else:v=np.ones(d)/np.sqrt(d)
                if not np.isfinite(z).all() or np.max(np.abs(z))>1e5:escaped=True;break
                total+=np.log(max(norm,1e-300));s+=z;ss+=z*z;maximum=max(maximum,np.max(np.abs(z)));minimum=min(minimum,np.max(np.abs(z)))
                completed[lane]+=1;index=b*block+k-(block*blocks-keep)
                if index>=0:tails[lane,index]=z
            if escaped:break
            stats[lane,b,0]=total/block;stats[lane,b,1]=minimum;stats[lane,b,2]=maximum;stats[lane,b,3]=residual
            stats[lane,b,4:4+d]=s/block;stats[lane,b,4+d:]=np.sqrt(np.maximum(0,ss/block-(s/block)**2));checkpoints[lane,b]=z
    return stats,checkpoints,tails,completed

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--continuation',action='store_true')
    ap.add_argument('--block',type=int,default=50000);ap.add_argument('--blocks',type=int,default=20);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    for m,n in [(2,2),(2,3),(3,2)]:
        parent=args.input/f'{m}x{n}.json';data=json.loads(parent.read_text())
        if args.continuation:
            rows=[r for r in data['runs'] if r['completed_steps']==data['block_length']*data['blocks']];starts=np.array([r['blocks'][-1]['checkpoint'] for r in rows]);origin=data['stop_step']
        else:
            rows=[r for r in data['runs'] if r['rates'] and r['rates']['mean']>.005]
            endpoints=np.load(args.input/f'{m}x{n}-endpoints-and-last2048.npz')['endpoints'];starts=np.array([endpoints[r['case']] for r in rows]);origin=30000
        L=np.array([r['log_kernel'] for r in rows]);ws=np.array([[r['wx'],r['wy']] for r in rows]);t=time.perf_counter();d=m+n-2
        stats,checkpoints,tails,completed=follow(L,starts,ws,args.block,args.blocks,50000);out=[]
        for i,r in enumerate(rows):
            out.append(dict(case=r['case'],shape=[m,n],wx=r['wx'],wy=r['wy'],log_kernel=r['log_kernel'],starting_state=starts[i].tolist(),completed_steps=int(completed[i]),
                blocks=[dict(stop_step=origin+(j+1)*args.block,rate=float(s[0]),min_abs_state=float(s[1]),max_abs_state=float(s[2]),final_step_residual=float(s[3]),mean=s[4:4+d].tolist(),std=s[4+d:].tolist(),checkpoint=checkpoints[i,j].tolist()) for j,s in enumerate(stats[i]) if np.isfinite(s).all()]))
        result=dict(runs=out,start_step=origin,stop_step=origin+args.block*args.blocks,block_length=args.block,blocks=args.blocks,retained_final_window_steps=50000,parent_sha256=hashlib.sha256(parent.read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-t,
            sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'))})
        (args.output/f'{m}x{n}.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        np.savez_compressed(args.output/f'{m}x{n}-final-50000-states.npz',case=np.array([r['case'] for r in rows]),states=tails)
        print(json.dumps(dict(shape=[m,n],elapsed=result['elapsed_seconds'],cases=[dict(case=r['case'],completed=r['completed_steps'],last=r['blocks'][-1]['rate'] if r['blocks'] else None,norm=r['blocks'][-1]['max_abs_state'] if r['blocks'] else None) for r in out])),flush=True)
