"""Longer binary64 confirmation of selected primary and companion portraits."""
import argparse
import json
from pathlib import Path
import numpy as np

from confirm import periods
from model import replay
from reproduce import raw_parameters


def main(args):
    d=json.loads(args.record.read_text());raw=raw_parameters(d)
    tail,rates,minimum,radius,visits=replay(raw,raw[1,:3].copy(),args.burn,args.window,args.windows,65536,True)
    p,res=periods(tail)
    result=dict(case=d['case'],record=args.record.name,burn=args.burn,window_length=args.window,windows=args.windows,qr_windows=rates.tolist(),mean=rates.mean(0).tolist(),minimum=rates.min(0).tolist(),maximum=rates.max(0).tolist(),period=p,minimum_period_residual=res,sampled_minimum_axis=float(minimum[0]),sampled_minimum_centerline=float(minimum[1]),sampled_maximum_norm=radius,initial=d['initial'],source_record_spectrum=d['rates'])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(d['case'],args.window*args.windows,result['mean'],result['minimum'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--record',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--burn',type=int,default=100000);p.add_argument('--window',type=int,default=5000000);p.add_argument('--windows',type=int,default=20)
    main(p.parse_args())
