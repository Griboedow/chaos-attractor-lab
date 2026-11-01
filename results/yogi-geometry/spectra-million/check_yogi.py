"""Full three-direction branch spectra for selected scalar-quadratic Yogi runs."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from yogi_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    d=np.load(args.source);ids=np.array(args.ids,dtype=int);parameters=d['parameters'][ids];starts=d['starts'][ids].copy()
    if args.initial_v is not None:starts[:,2]=args.initial_v
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),
        selection='Visually selected cases44,2220,7591,5174,4214,7158; known later-collapse control4454 and default-pair controls8192/8194.',
        starts=starts.tolist(),bias_corrected_warmup=args.bias_warm,warmup=100000,observed=args.observed,window=50000,tail=100000,
        tangent='Full3D branch derivative, repeated modified Gram-Schmidt each step; determinant beta1 supplies an independent sum check.',
        scope='The map is discontinuous at v=x^2. Finite branch derivative estimates do not prove asymptotic chaotic attraction or historical novelty.')
    assert args.observed>=100000 and args.observed%50000==0
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_yogi.py','yogi_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    before=time.perf_counter();windows,tails,ends,margins=spectra(parameters,starts,100000,args.observed,50000,100000,args.bias_warm)
    assert np.isfinite(windows).all() and np.isfinite(tails).all()
    sums=windows.sum(2)-np.log(parameters[:,0,None]);assert np.max(abs(sums))<1e-9
    file=out/'histories.npz';np.savez_compressed(file,ids=ids,parameters=parameters,starts=starts,windows=windows,tails=tails,ends=ends,minimum_margins=margins,determinant_sum_residuals=sums)
    rows=[]
    for k,ident in enumerate(ids):
        tail=tails[k,-4096:];returns=np.array([np.max(np.linalg.norm(tail[lag:]-tail[:-lag],axis=1)) for lag in range(1,129)])
        period=next((i+1 for i,v in enumerate(returns) if v<1e-9*(1+np.max(abs(tail)))),None)
        row=dict(id=int(ident),parameters=parameters[k].tolist(),spectrum=windows[k].mean(0).tolist(),
            leading_window_range=[float(windows[k,:,0].min()),float(windows[k,:,0].max())],
            second_window_range=[float(windows[k,:,1].min()),float(windows[k,:,1].max())],
            tail_extrema=[tails[k].min(0).tolist(),tails[k].max(0).tolist()],tail_period_up_to128=period,
            determinant_sum_max_residual=float(np.max(abs(sums[k]))),minimum_scaled_switch_margin=float(margins[k]),return_residuals=returns.tolist())
        rows.append(row);print(ident,row['spectrum'],'minleading',row['leading_window_range'][0],'period',period,flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,completed=len(rows),declared=len(ids),seconds=time.perf_counter()-before,sha256=hashlib.sha256(file.read_bytes()).hexdigest()),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--observed',type=int,default=1000000)
    p.add_argument('--bias-warm',action='store_true');p.add_argument('--initial-v',type=float);main(p.parse_args())
