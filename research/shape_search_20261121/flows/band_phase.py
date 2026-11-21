"""Scan coefficients adjacent to stable knotted cycles for chaotic bands."""
import argparse
import csv
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.spatial import cKDTree
from field import integrate,ring_tangent


def phase_measurement(points,anchor):
    scale=np.sqrt(np.mean(np.sum((anchor-anchor.mean(0))**2,axis=1)))
    tree=cKDTree(anchor)
    dist,idx=tree.query(points,workers=1)
    phase=2*np.pi*idx/len(anchor)
    increments=np.diff(np.unwrap(phase))/(2*np.pi)
    counts,_=np.histogram(idx,np.linspace(0,len(anchor),33))
    return dict(distance95=float(np.quantile(dist,.95)/scale),median_distance=float(np.median(dist)/scale),
                phase_turns=float(increments.sum()),strong_phase_jumps=int(np.sum(abs(increments)>.08)),
                phase_coverage=float(np.mean(counts>0)),phase_increment95=float(np.quantile(abs(increments),.95)))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--ids',nargs='+',type=int,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--points',type=int,default=65)
    ap.add_argument('--span',type=float,default=.45)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    rows=[];inputs={}
    with (args.out/'scan.csv').open('w',newline='') as out:
        writer=csv.writer(out)
        writer.writerow(['case','base','parameter','fraction','minimum_rate','mean_rate','distance95','jumps','coverage'])
        for group,base in enumerate(args.ids):
            ref=json.loads((args.source/f'case-{base}.json').read_text())
            anchor=np.load(args.source/f'case-{base}.npz')['points'][::4]
            original=np.array(ref['parameters'])
            initial=np.array(ref['refined_start'])
            for par in (2,3,4,5,6,7,8,9,10):
                for step,fraction in enumerate(np.linspace(-args.span,args.span,args.points)):
                    if abs(fraction)<1e-12:
                        continue
                    case=160000+group*10000+par*100+step
                    p=original.copy();p[par]+=fraction*max(.45,abs(p[par]))
                    if p[7]<.02:
                        continue
                    rates,z=integrate(p,initial,.02,10000,6000,4,4,ring_tangent)
                    q95=None;jumps=None;coverage=None
                    if min(rates)>.008:
                        stats=phase_measurement(z,anchor)
                        q95=stats['distance95'];jumps=stats['strong_phase_jumps'];coverage=stats['phase_coverage']
                        if q95<.5 and coverage==1:
                            rows.append(dict(case=case,base=base,parameter_index=par,fraction=float(fraction),
                                             windows=rates.tolist(),**stats))
                            inputs[str(case)]=dict(parameters=p.tolist(),initial=initial.tolist())
                    writer.writerow([case,base,par,fraction,min(rates),rates.mean(),q95,jumps,coverage])
                out.flush()
            (args.out/'leads.json').write_text(json.dumps(rows,indent=2)+'\n')
            (args.out/'inputs.json').write_text(json.dumps(inputs,separators=(',',':'))+'\n')
            print(base,len(rows),flush=True)
    (args.out/'settings.json').write_text(json.dumps({**vars(args),'out':str(args.out),'source':str(args.source)},indent=2)+'\n')


if __name__=='__main__':
    main()
