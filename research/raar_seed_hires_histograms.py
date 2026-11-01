"""Recompute a denser image grid, checking every state snapshot and old pixel count."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import set_num_threads
from trace_raar_cycle_seeds import advance,phases


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(8);report=json.loads((args.source/'report.json').read_text());p=np.array(report['params']);x=np.array(report['cycle_seed']);direction=np.array(report['unstable_eigenvector']);states=x+np.linspace(-1e-7,1e-7,16384)[:,None]*direction
    with np.load(args.source/'radius-1e-07-n16384.npz') as d:stored=d['snapshots'];returns=d['returns']
    with np.load(args.source/'evolution-histograms.npz') as d:old=d['histograms']
    meta=report['histogram'];center=np.array(meta['camera_center']);projection=np.array(meta['projection']);histograms=[];checks=[]
    for k in range(65):
        phase,status,_=phases(states,p);assert not status.any()
        if k in returns:assert np.array_equal(phase,stored[list(returns).index(k)]);checks.append(k)
        xy=(phase.reshape(-1,3)-center)@projection;h,_,_=np.histogram2d(xy[:,0],xy[:,1],bins=(1800,1400),range=[meta['x_bounds'],meta['y_bounds']]);h=h.T[::-1].astype(np.uint32);assert h.sum()==7*16384
        assert np.array_equal(h.reshape(700,2,900,2).sum((1,3)),old[k]);histograms.append(h)
        if k<64:states,status,_=advance(states,p,7);assert not status.any()
    np.savez_compressed(args.out/'evolution-histograms.npz',histograms=np.array(histograms),returns=np.arange(65),clipped=np.zeros(65,dtype=int))
    for name in ['raar_seed_hires_histograms.py','trace_raar_cycle_seeds.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    row=dict(status='passed',source=str(args.source),source_report_sha256=hashlib.sha256((args.source/'report.json').read_bytes()).hexdigest(),exact_snapshot_checks=checks,all_65_hires_histograms_sum_exactly_to_original_pixels=True,resolution=[1800,1400],scope='Rendering resolution update only. All ten stored state snapshots are bit-identical and every 2x2 pixel aggregation matches the original complete histogram.')
    (args.out/'report.json').write_text(json.dumps(row,indent=2),encoding='utf-8');print(json.dumps(row),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
