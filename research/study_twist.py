"""Reproduce long full spectra, fresh starts and a local parameter cohort."""
import argparse
import json
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.stats import qmc
from twist_map import spectra


def save_run(path, ids, parameters, starts, warmup, observed, window, keep):
    windows, tail, end, means, radii = spectra(parameters, starts, warmup, observed, window, keep)
    np.savez_compressed(path, ids=ids, parameters=parameters, starts=starts,
                        warmup=warmup, observed=observed, window=window,
                        windows=windows, tail=tail, end=end, means=means, radii=radii)
    print(json.dumps({'file':str(path),'lanes':len(ids),
                      'positive_all_windows':int((windows[:,:,0].min(1)>0).sum())}),flush=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--directory',type=Path,default=Path('results/twist-geometry'))
    ap.add_argument('--stage',choices=['long','basins','neighbors','all'],default='all')
    args=ap.parse_args()
    set_num_threads(4)
    candidates=json.loads((args.directory/'candidates.json').read_text())['candidates']
    ids=np.array([c['id'] for c in candidates])
    p=np.array([c['parameters'] for c in candidates])
    if args.stage in ['long','all']:
        save_run(args.directory/'long-runs.npz',ids,p,np.tile([.1,.2,.3],(3,1)),100000,4000000,50000,200000)
    if args.stage in ['basins','all']:
        starts=np.random.default_rng(20261105).uniform(-2,2,(96,3))
        save_run(args.directory/'basins.npz',np.repeat(ids,32),np.repeat(p,32,axis=0),starts,100000,500000,25000,8192)
    if args.stage in ['neighbors','all']:
        # Perturb each of the seven/ten active parameters by <= 0.001.
        offsets=(qmc.Sobol(10,scramble=True,seed=20261107).random_base2(5)*2-1)*.001
        rows=[]
        for row in p[:2]:
            change=offsets.copy()
            if np.all(row[7:]==0):
                change[:,7:]=0
            rows.append(np.vstack([row+change,np.round(row,3)]))
        pp=np.vstack(rows)
        save_run(args.directory/'neighbors.npz',np.repeat(ids[:2],33),pp,np.tile([.1,.2,.3],(66,1)),100000,500000,25000,4096)


if __name__=='__main__':
    main()
