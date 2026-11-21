"""Long forward portraits of contrasting chaotic basins for identical parameters."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from model import replay
from reproduce import raw_parameters
from confirm import periods


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    args.cache.mkdir(parents=True,exist_ok=True)
    for path in sorted(args.source.glob('case-*.json')):
        rec=json.loads(path.read_text())
        if not rec['three_positive_passes']:continue
        case=rec['case'];dest=args.out/f'companion-{case}.json'
        if dest.exists():continue
        choices=[]
        for i,row in enumerate(rec['cohort']):
            spec=np.array(row['qr_windows'])
            if row['period'] is None and spec[:,0].min()>.03 and spec[:,1].max()<-.015:
                choices.append(i)
        if not choices:continue
        state=np.array(rec['cohort'][choices[0]]['initial'])
        raw=raw_parameters(rec)
        cloud,rates,_,_,_=replay(raw,state,100000,500000,10,250000,True)
        p,res=periods(cloud)
        if p or rates[:,0].min()<.03 or rates[:,1].max()>-.015:continue
        cache=args.cache/f'companion-{case}.npz'
        np.savez_compressed(cache,points=cloud.astype('f4'),case=case,initial=state)
        with np.load(rec['cloud']) as f:primary=f['points'].astype(float)
        radius=np.sqrt(np.var(np.vstack((primary,cloud)),axis=0).sum())
        def distances(a,b):
            q=cKDTree(b).query(a[::7])[0]/radius
            return np.quantile(q,[0,.01,.5,.99,1]).tolist()
        doc=dict(case=case,parameters=rec['parameters'],initial=state.tolist(),cohort_indices=choices,cloud=str(cache.resolve()),cloud_relative='projection/three-companions/'+cache.name,cloud_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),rates=rates.mean(0).tolist(),qr_windows=rates.tolist(),burn=100000,window_length=500000,windows=10,points_saved=len(cloud),detected_period=p,minimum_period_residual=res,nearest_distance_quantiles=[0,.01,.5,.99,1],distances_shared_rms=dict(primary_to_companion=distances(primary,cloud),companion_to_primary=distances(cloud,primary),primary_temporal=distances(primary[len(primary)//2:],primary[:len(primary)//2]),companion_temporal=distances(cloud[len(cloud)//2:],cloud[:len(cloud)//2])))
        dest.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
        print(case,choices,doc['rates'],doc['distances_shared_rms'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--cache',required=True,type=Path)
    main(p.parse_args())
