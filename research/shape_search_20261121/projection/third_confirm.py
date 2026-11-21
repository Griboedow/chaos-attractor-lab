"""Run the long three-rate controls for configurations retained by the CPU rescreen."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from confirm import main as confirm,save_json


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    inspection=args.cache/'third-inspection'
    inspection.mkdir(parents=True,exist_ok=True)
    for geometry in args.geometries:
        report=json.loads((args.source/f'three-rate-{geometry}.json').read_text())
        with np.load(args.cache/'three-deep'/f'deep-{geometry}.npz') as f:
            ids=f['ids'];raw=f['data'];kinds=f['kind']
        for row in report['cases']:
            if row['min_cpu_third']<=.001:continue
            i=int(np.flatnonzero(ids==row['index'])[0])
            offset={'single':3000000,'linked':4000000,'smooth':5000000,'evolved':6000000}[geometry]
            case=offset+int(ids[i]);dest=args.out/f'case-{case}.json'
            if dest.exists():continue
            np.savez(inspection/f'case-{case}.npz',data=raw[:,i,:],kind=kinds[i],rates=np.asarray(row['cpu_windows']))
            # Isolate each normal forward-run result from aggregate indexes.
            tmp_out=inspection/str(case)
            confirm(SimpleNamespace(out=tmp_out,cache=args.cache/'final',inspection=inspection,cases=[case],burn=100000,window=500000,windows=10,tail=250000,starts=17,start_burn=25000,start_window=16384,cohort_qr=True))
            record=json.loads((tmp_out/f'case-{case}.json').read_text())
            starts=sum(min(v[2] for v in m['qr_windows'])>.001 and m['period'] is None for m in record['cohort'])
            record.update(kind='Perturbed '+record['kind'],parent_case=row['parent'],third_growth_stream=geometry,stream_index=int(ids[i]),cohort_three_positive=starts,minimum_third_window=min(v[2] for v in record['qr_windows']))
            if geometry=='evolved':record['added_posterior_sites']=row['added']
            record['three_positive_passes']=bool(record['passes'] and record['minimum_third_window']>.001 and starts>=9)
            save_json(dest,record)
            print('three long',case,record['three_positive_passes'],starts,record['rates'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--geometries',nargs='+',default=['single','linked','smooth'])
    main(p.parse_args())
