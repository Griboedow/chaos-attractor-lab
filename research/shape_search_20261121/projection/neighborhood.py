"""Local sensitivity in tube radius and posterior temperature around recorded cases."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from model import replay


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    for path in sorted(args.source.glob('case-*.json')):
        record=json.loads(path.read_text())
        if not record.get('three_positive_passes',record['passes']): continue
        target=args.out/f'near-{record["case"]}.json'
        if target.exists(): continue
        params=record['parameters']
        raw=np.zeros((10,4))
        sites=np.array(params['sites']);n=len(sites)
        raw[0]=[params['beta'],params['inv_tau'],params['tube_radius'],n]
        raw[1,:3]=record['initial']; raw[1,3]=params['second_ring']
        raw[2:2+n,:3]=sites
        if 'mix_temperature' in params:
            raw[1,3]=2
            raw[2,3]=params['mix_temperature']
        rows=[]
        for temperature in args.factors:
            for radius in args.factors:
                case=raw.copy()
                case[0,1]/=temperature
                case[0,2]*=radius
                if raw[1,3] and case[0,2]>=.5: continue
                z,rates,_,_,visits=replay(case,raw[1,:3].copy(),args.burn,args.window,4,8192,True)
                eig=np.linalg.eigvalsh(np.cov(z.T))
                eigenfractions=(eig/eig.sum()).tolist() if eig.sum()>1e-20 else None
                rows.append(dict(tau_factor=temperature,tube_factor=radius,rates=rates.mean(axis=0).tolist(),min_leading=float(rates[:,0].min()),min_second=float(rates[:,1].min()),min_third=float(rates[:,2].min()),branch_fraction=float(visits[1]/visits.sum()),pca_fractions=eigenfractions))
        out=dict(case=record['case'],burn=args.burn,window_length=args.window,windows=4,initial=record['initial'],factors=args.factors,scope='One initial condition per grid point; positive finite-time rates and geometric variation only.',rows=rows)
        tmp=target.with_suffix('.tmp.json')
        tmp.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
        os.replace(tmp,target)
        print(record['case'],sum(r['min_leading']>.015 for r in rows),sum(r['min_second']>.015 for r in rows),'/',len(rows),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--factors',nargs='+',type=float,default=[.94,.98,1.,1.02,1.06])
    p.add_argument('--burn',type=int,default=16384)
    p.add_argument('--window',type=int,default=16384)
    main(p.parse_args())
