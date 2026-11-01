"""Declared cold-start Yogi parameter slices, retaining all three branch rates."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from yogi_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    designs=[dict(name='two-direction-neighborhood',b=np.linspace(.55,.995,32),c=np.linspace(.02,.92,32),e=9.301461332493893e-5),
        dict(name='fan-region',b=np.linspace(.005,.5,32),c=np.linspace(.05,.95,32),e=.1)]
    rows=[];p=[];slice_ids=[];grid=[]
    for k,design in enumerate(designs):
        rows.append(dict(name=design['name'],beta1=design['b'].tolist(),beta2=design['c'].tolist(),epsilon_normalized=design['e']))
        for j,c in enumerate(design['c']):
            for i,b in enumerate(design['b']):p.append([b,c,design['e']]);slice_ids.append(k);grid.append([j,i])
    p=np.array(p);starts=np.tile([.73,0.,0.],(len(p),1));protocol=dict(slices=rows,declared=len(p),start=[.73,0.,0.],
        warmup=100000,observed=200000,window=50000,tail=2048,bias_correction=False,
        design='Two explicitly declared32x32rectangular grids, fixed normalizedepsilon per slice; independent cold start at every point, no parameter continuation.',
        scope='Finite parameter slices do not prove open hyperchaotic regions or bifurcation mechanisms. The second rate is reported with all four window extrema, not just its mean.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['scan_yogi_plane.py','yogi_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    before=time.perf_counter();blocks=[]
    for offset in range(0,len(p),128):
        result=spectra(p[offset:offset+128],starts[offset:offset+128],100000,200000,50000,2048)
        np.savez_compressed(out/f'block-{offset:04d}.npz',parameters=p[offset:offset+128],starts=starts[offset:offset+128],windows=result[0],tails=result[1],ends=result[2],minimum_margins=result[3])
        assert all(np.isfinite(v).all() for v in result)
        assert np.max(abs(result[0].sum(2)-np.log(p[offset:offset+128,0,None])))<1e-9
        blocks.append(result);print('Completed',offset+len(result[0]),'/',len(p),'seconds',time.perf_counter()-before,flush=True)
    data=[np.concatenate([r[k] for r in blocks]) for k in range(4)];w,tail,ends,margin=data
    np.savez_compressed(out/'histories.npz',parameters=p,starts=starts,slice_ids=slice_ids,grid_indices=grid,windows=w,tails=tail,ends=ends,minimum_margins=margin)
    periods=[]
    for z in tail:
        tol=1e-9*(1+np.max(abs(z)));periods.append(next((k for k in range(1,129) if np.max(np.linalg.norm(z[k:]-z[:-k],axis=1))<tol),0))
    periods=np.array(periods);leading=(w[:,:,0].min(1)>.005)&(np.ptp(tail[:,:,0],axis=1)>1e-4)&(periods==0)
    second=leading&(w[:,:,1].min(1)>.001);summaries=[]
    for k,design in enumerate(designs):
        mask=np.array(slice_ids)==k;summaries.append(dict(name=design['name'],count=int(mask.sum()),leading_pass=int(np.sum(leading&mask)),two_pass=int(np.sum(second&mask)),periodic=int(np.sum((periods>0)&mask))))
    (out/'analysis.json').write_text(json.dumps(dict(completed=len(p),declared=len(p),slices=summaries,leading_ids=np.flatnonzero(leading).tolist(),two_ids=np.flatnonzero(second).tolist(),periods=periods.tolist(),wall_seconds=time.perf_counter()-before,
        histories_sha256=hashlib.sha256((out/'histories.npz').read_bytes()).hexdigest()),indent=2)+'\n');print(summaries,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
