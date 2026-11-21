"""Full binary64 QR replay of every GPU set with two positive tangent windows."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from model import replay


def main(args):
    args.cache.mkdir(parents=True,exist_ok=True)
    rows=[]
    for path in sorted(args.source.glob('deep-*.npz')):
        with np.load(path) as f:
            items=f['data'].transpose(1,0,2).astype(float)
            ids=f['ids'];kinds=f['kinds'];gpu=f['metrics']
        for j,(number,raw) in enumerate(zip(ids,items)):
            case=int(number)
            z,rates,m,rad,branches=replay(raw,raw[1,:3].copy(),args.burn,args.window,4,24576,True)
            passed=bool(np.min(rates[:,:2])>.015)
            row=dict(case=case,passes=passed,gpu_second_mean=float(np.mean(gpu[-4:,j])),cpu_qr=rates.tolist(),branch_visits=branches.tolist())
            rows.append(row)
            if passed:
                target=args.cache/f'case-{case}.npz'
                np.savez(target,points=z.astype(np.float32),data=raw,kind=kinds[j],case=case,rates=rates,minima=m,max_norm=rad,branches=branches)
    result=dict(burn=args.burn,window_length=args.window,windows=4,gpu_flags=len(rows),cpu_passes=sum(v['passes'] for v in rows),cases=rows)
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(result['gpu_flags'],result['cpu_passes'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--burn',type=int,default=32768)
    p.add_argument('--window',type=int,default=32768)
    main(p.parse_args())
