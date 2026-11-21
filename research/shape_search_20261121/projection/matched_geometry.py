"""Replay both linked-surface operators at identical parameters and initial states."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from model import replay
from screen import parameters
from reproduce import raw_parameters
from confirm import periods

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from compare_attractor_shapes import fingerprint,distance


def main(args):
    candidates={}
    for p in args.source.glob('case-*.json'):
        rec=json.loads(p.read_text())
        if rec['passes'] and rec['parameters']['second_ring']:
            candidates[rec['case']%1000000]=rec
    rows=[]
    for base,rec in sorted(candidates.items()):
        batch,ix=divmod(base,4096)
        data,_=parameters(batch,4096,'smooth')
        soft=data[:,ix,:].astype(float)
        old=raw_parameters(rec)
        assert np.array_equal(old[0],soft[0]) and np.array_equal(old[2:,:3],soft[2:,:3])
        outcomes=[];fps=[]
        for mode in (1,2):
            raw=soft.copy();raw[1,3]=mode
            z,rates,_,_,visits=replay(raw,raw[1,:3].copy(),args.burn,args.window,4,65536,True)
            period,_=periods(z)
            spread=float(np.sqrt(np.mean(np.sum((z-z.mean(0))**2,axis=1))))
            fp=fingerprint(z) if spread>1e-9 else None
            halves=float(distance(fingerprint(z[:len(z)//2]),fingerprint(z[len(z)//2:],1))) if fp is not None else None
            outcomes.append(dict(mode='hard' if mode==1 else 'blended',qr=rates.tolist(),period=period,split_distance=halves,rms_radius=spread,nearest_second_fraction=float(visits[1]/visits.sum())))
            fps.append(fp)
        rows.append(dict(base=base,hard_case=1000000+base,blended_case=2000000+base,mix_temperature=float(soft[2,3]),fingerprint_distance=float(distance(*fps)) if all(v is not None for v in fps) else None,outcomes=outcomes))
    args.out.write_text(json.dumps(dict(burn=args.burn,window_length=args.window,windows=4,tail_length=65536,scope='Matched deterministic configurations. Fingerprint distance is an empirical equal-unit geometric score.',cases=rows),indent=2)+'\n',encoding='utf8')
    print(len(rows),'matched pairs',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--burn',type=int,default=65536)
    p.add_argument('--window',type=int,default=65536)
    main(p.parse_args())
