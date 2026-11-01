"""Long full-spectrum follow-up of every declared tanh-Yogi scout flag."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from yogi_tanh_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    with np.load(args.source) as a:d={k:a[k] for k in a.files}
    mask=(d['status']==0)&(d['windows'].min(1)>.005)&(np.ptp(d['tails'][:,:,0],axis=1)>1e-4)&(d['periods']==0)
    ids=np.flatnonzero(mask) if args.ids is None else np.array(args.ids)
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),
        selection='All primary scout flags, or explicit IDs supplied on the command line.',warmup=100000,observed=args.observed,window=50000,tail=4096,
        initialization='Fresh original scout starts, retaining each prescribed bias-correction setting throughout. No continuation from the scout endpoint.',
        scope='Finite smooth-map full derivative spectra. No theorem of chaos, infinite-time positivity, attractor existence or historical uniqueness.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['follow_yogi_tanh.py','yogi_tanh_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    Q=np.linalg.qr(np.array([[1.,2.,3.],[4.,5.,7.],[2.,6.,8.]]))[0]
    before=time.perf_counter();blocks=[];rows=[]
    for lo in range(0,len(ids),32):
        chosen=ids[lo:lo+32];p=d['parameters'][chosen];z=d['starts'][chosen];bias=d['bias'][chosen]
        w,det,tail,end,status,negative=spectra(p,z,bias,100000,args.observed,50000,4096,Q)
        block=dict(ids=chosen,parameters=p,starts=z,bias=bias,windows=w,determinant_windows=det,tails=tail,ends=end,status=status,negative_determinant_counts=negative)
        np.savez_compressed(out/f'block-{lo:04}.npz',**block);blocks.append(block)
        for i,ident in enumerate(chosen):
            if status[i]:rows.append(dict(id=int(ident),status=int(status[i])));continue
            tol=1e-9*(1+np.max(abs(tail[i])))
            period=next((lag for lag in range(1,129) if np.max(abs(tail[i,lag:]-tail[i,:-lag]))<tol),None)
            residual=float(np.max(abs(w[i].sum(1)-det[i])));assert residual<1e-7
            leading=bool(w[i,:,0].min()>.005 and np.ptp(tail[i,:,0])>1e-4 and period is None)
            rows.append(dict(id=int(ident),status=0,spectrum=w[i].mean(0).tolist(),window_minima=w[i].min(0).tolist(),period_up_to128=period,
                leading=leading,two=bool(leading and w[i,:,1].min()>.001),determinant_sum_max_residual=residual,negative_determinant_fraction=float(negative[i]/args.observed)))
        print('followed',min(lo+32,len(ids)),'/',len(ids),'seconds',round(time.perf_counter()-before,1),flush=True)
    joined={k:np.concatenate([r[k] for r in blocks]) for k in blocks[0]};np.savez_compressed(out/'histories.npz',**joined)
    result=dict(declared=len(ids),completed=sum(r['status']==0 for r in rows),leading=sum(r.get('leading',False) for r in rows),two=sum(r.get('two',False) for r in rows),runs=rows,seconds=time.perf_counter()-before)
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='runs'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--observed',type=int,default=1000000);p.add_argument('--ids',type=int,nargs='+');main(p.parse_args())
