"""Compare occupation across independent starts in the same physical coordinates."""
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np


def cloud_stats(z, directions, shift, scale):
    # Quantile Wasserstein approximation with fixed axes and the anchor's scale.
    z = z[::max(1,len(z)//8000)]
    units = (z-shift)/scale
    qs = np.linspace(.01, .99, 99)
    return np.quantile(units @ directions.T,qs,axis=0)


def rotation(z, n, dt):
    theta = np.arctan2(z[:,1],z[:,0])
    increments = np.diff(np.unwrap(theta))
    dt_total = (len(theta)-1)*dt
    turns = increments.sum()/2/np.pi
    counts,_ = np.histogram(np.mod(theta,2*np.pi),np.linspace(0,2*np.pi,n+1))
    return dict(symmetry_order=n, total_net_turns=float(turns), net_turns_per_time=float(turns/dt_total),
                positive_increment_fraction=float(np.mean(increments>0)),
                angular_distance_signed_over_total=float(increments.sum()/np.abs(increments).sum()),
                fraction_per_equal_sector=(counts/counts.sum()).tolist(),
                max_sampled_angular_increment=float(np.max(np.abs(increments))))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--selections',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    settings=json.loads(args.selections.read_text(encoding='utf-8-sig'))
    rng=np.random.default_rng(48293)
    directions=rng.normal(size=(32,3))
    directions/=np.linalg.norm(directions,axis=1,keepdims=True)
    report={}
    for bundle in settings:
        for case,_ in bundle['cases']:
            with np.load(args.cache/bundle['dense']/f'case-{case}.npz') as f:
                z=f['points']
                p=f['parameters']
                dt=float(f['step_size']*f['sample_stride'])
            mean=z.mean(0)
            scale=np.sqrt(np.mean(np.sum((z-mean)**2,axis=1)))
            anchor=cloud_stats(z,directions,mean,scale)
            halves=[cloud_stats(v,directions,mean,scale) for v in np.array_split(z,2)]
            result=dict(anchor_mean=mean.tolist(),anchor_rms_radius=float(scale),
                        split_halves_mean_sliced_distance=float(np.abs(halves[0]-halves[1]).mean()),
                        control_mean_sliced_distance={})
            for name in bundle['controls']:
                with np.load(args.cache/name/f'case-{case}.npz') as f:
                    other=cloud_stats(f['points'],directions,mean,scale)
                result['control_mean_sliced_distance'][name]=float(np.abs(other-anchor).mean())
            if bundle['family']=='ring':
                result['rotation']=rotation(z,int(p[0]),dt)
            report[str(case)]=result
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
