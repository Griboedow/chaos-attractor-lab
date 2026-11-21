"""Recompute one stored configuration from its compact JSON record."""
import argparse
import json
from pathlib import Path
import numpy as np
from model import replay


def raw_parameters(record):
    p=record['parameters']
    sites=np.asarray(p['sites'],dtype=np.float64)
    raw=np.zeros((10,4),np.float64)
    raw[0]=[p['beta'],p['inv_tau'],p['tube_radius'],len(sites)]
    raw[1,:3]=record['initial']
    raw[1,3]=2 if 'mix_temperature' in p else int(p['second_ring'])
    raw[2:2+len(sites),:3]=sites
    raw[2,3]=p.get('mix_temperature',0)
    return raw


def main(args):
    rec=json.loads(args.record.read_text())
    raw=raw_parameters(rec)
    z,rates,_,_,_=replay(raw,raw[1,:3].copy(),rec['burn'],rec['window_length'],rec['windows'],rec['points_saved'],True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.out,points=z.astype(np.float32),case=rec['case'],initial=raw[1,:3])
    details=dict(case=rec['case'],count=len(z),rates=rates.mean(axis=0).tolist(),max_window_abs_difference=float(np.max(abs(rates-np.asarray(rec['qr_windows'])))))
    if args.reference:
        with np.load(args.reference) as f: existing=f['points']
        details['identical_float32_cloud']=bool(np.array_equal(z.astype('f4'),existing))
    print(json.dumps(details))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--record',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--reference',type=Path)
    main(p.parse_args())
