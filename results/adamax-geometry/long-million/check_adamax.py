"""Long two-direction tangent checks of explicitly selected Adamax geometries."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from adamax_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    d=np.load(args.source);ids=np.array(args.ids,dtype=int);parameters=d['parameters'][ids];starts=d['starts'][ids]
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),
        selection='Four positive explicit controls and visually selected gallery02cases. Selection precedes these long calculations.',
        starts='Fresh zero-memory cold starts from original scout design; continuously bias-corrected warmup.',
        warmup=100000,observed=args.observed,window=50000,tail=100000,
        tangent='Two surviving directions, each-step repeated modified Gram-Schmidt. Any reset erases the input memory direction exactly; third derivative-product exponent is -infinity after a reset, not estimated by floating-point QR.',
        scope='Finite-time numerical expansion and geometry; neither novelty nor asymptotic chaotic attraction is established.')
    assert args.observed>=100000 and args.observed%50000==0
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_adamax.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    before=time.perf_counter();windows,tails,ends,resets,margins=spectra(parameters,starts,100000,args.observed,50000,100000,True)
    assert np.isfinite(windows).all() and np.isfinite(tails).all()
    file=out/'histories.npz';np.savez_compressed(file,ids=ids,parameters=parameters,starts=starts,windows=windows,tails=tails,ends=ends,reset_counts=resets,minimum_margins=margins)
    rows=[]
    for k,ident in enumerate(ids):
        tail=tails[k,-2048:];returns=np.array([np.max(np.linalg.norm(tail[lag:]-tail[:-lag],axis=1)) for lag in range(1,65)])
        period=next((i+1 for i,v in enumerate(returns) if v<1e-10*(1+np.linalg.norm(np.ptp(tail,axis=0)))),None)
        row=dict(id=int(ident),parameters=parameters[k].tolist(),spectrum_two=windows[k].mean(0).tolist(),
            leading_window_range=[float(windows[k,:,0].min()),float(windows[k,:,0].max())],
            tail_extrema=[tails[k].min(0).tolist(),tails[k].max(0).tolist()],tail_period_up_to64=period,
            reset_count=int(resets[k]),reset_fraction=float(resets[k]/args.observed),minimum_scaled_switch_margin=float(margins[k]),
            return_residuals=returns.tolist())
        rows.append(row);print(ident,row['spectrum_two'],'minleading',row['leading_window_range'][0],'period',period,flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,completed=len(rows),declared=len(ids),seconds=time.perf_counter()-before,sha256=hashlib.sha256(file.read_bytes()).hexdigest()),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--observed',type=int,default=1000000);main(p.parse_args())
