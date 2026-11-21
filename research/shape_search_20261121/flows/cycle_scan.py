"""Find stable-cycle candidates whose sampled native curves have nontrivial determinant."""
import argparse
import csv
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from field import integrate,make_case,make_ring,make_rotor,make_planes,make_crossed
from field import velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent
from knot_diagram import alexander


def resample(cloud, n):
    polygon=np.vstack([cloud,cloud[0]])
    distance=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(polygon,axis=0),axis=1))]
    grid=np.linspace(0,distance[-1],n,endpoint=False)
    return np.column_stack([np.interp(grid,distance,a) for a in polygon.T])


def period_guess(z,dt):
    scale=np.sqrt(np.mean(np.sum((z-z.mean(0))**2,axis=1)))
    if scale<1e-4:
        return
    for axis in np.argsort(np.ptp(z,axis=0))[::-1]:
        level=np.median(z[:,axis])
        hits=np.flatnonzero((z[:-1,axis]<level)&(z[1:,axis]>=level))
        if len(hits)<4:
            continue
        fractions=(level-z[hits,axis])/(z[hits+1,axis]-z[hits,axis])
        cross=z[hits]+fractions[:,None]*(z[hits+1]-z[hits])
        times=(hits+fractions)*dt
        for lag in range(1,min(25,len(hits)-2)):
            comparison=np.linalg.norm(cross[-3:]-cross[-3-lag:-lag],axis=1)/scale
            periods=times[-3:]-times[-3-lag:-lag]
            if np.max(comparison)>.004 or np.ptp(periods)>.004*np.mean(periods):
                continue
            i,j=hits[-1-lag],hits[-1]
            if j-i<24:
                continue
            segment=np.vstack([cross[-1-lag],z[i+1:j+1]])
            return axis,cross[-1-lag],float(periods[-1]),resample(segment,220),float(np.max(comparison)),scale
    return None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--in',dest='sources',nargs='+',type=Path,default=[])
    ap.add_argument('--first',type=int)
    ap.add_argument('--count',type=int,default=0)
    ap.add_argument('--variants',nargs='+',type=int,choices=range(6))
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--family',choices=('cartesian','ring','rotor','planes','crossed'),required=True)
    ap.add_argument('--max',type=int)
    ap.add_argument('--seed',type=int,default=2026112123)
    ap.add_argument('--shuffle',type=int,default=294179)
    ap.add_argument('--dt',type=float,default=.025)
    ap.add_argument('--burn',type=float,default=360.)
    ap.add_argument('--time',type=float,default=320.)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    mode=('cartesian','ring','rotor','planes','crossed').index(args.family)
    gen=(make_case,make_ring,make_rotor,make_planes,make_crossed)[mode]
    rhs=(velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent)[mode]
    eligible=[]
    if args.first is not None:
        eligible.extend(i for i in range(args.first,args.first+args.count)
                        if args.variants is None or i%6 in args.variants)
    for path in args.sources:
        with (path/'screen.csv').open() as f:
            for row in csv.DictReader(f):
                keys=[k for k in row if k.startswith('window_')]
                growth=np.mean([float(row[k]) for k in keys])
                if float(row['variance_ratio_min'])>.003 and float(row['span'])>.25 and abs(growth)<.04:
                    eligible.append(int(row['case']))
    eligible=np.unique(eligible)
    np.random.default_rng(args.shuffle).shuffle(eligible)
    chosen=eligible[:args.max]
    spec={**vars(args),'sources':[str(p) for p in args.sources],'out':str(args.out),'eligible':len(eligible),'sampled':len(chosen)}
    (args.out/'run.json').write_text(json.dumps(spec,indent=2)+'\n')
    with (args.out/'screen.csv').open('w',newline='') as handle, (args.out/'hits.jsonl').open('w') as hit_file:
        writer=csv.writer(handle)
        writer.writerow(['case','status','period','determinant','crossings','closure_scale'])
        for i,ident in enumerate(chosen):
            case=int(ident);p,start=gen(args.seed,case)
            _,z=integrate(p,start,args.dt,round(args.burn/args.dt),round(args.time/args.dt),1,2,rhs)
            found=period_guess(z,2*args.dt) if args.first is None or np.linalg.norm(np.ptp(z,axis=0))>.25 else None
            row=[case,'unresolved','','','','']
            if found:
                axis,seed,period,closed,error,rms=found
                pd=alexander(closed,17)
                if pd is not None:
                    row=[case,'trivial-determinant',period,pd['determinant'],pd['crossings'],error]
                    if pd['determinant']>1:
                        second=alexander(closed,81)
                        if second is not None and second['determinant']==pd['determinant']:
                            row[1]='nontrivial-determinant'
                            record=dict(case=case,family=args.family,axis=int(axis),start=seed.tolist(),period=period,
                                        rms=rms,determinant=pd['determinant'],two_crossing_counts=[pd['crossings'],second['crossings']])
                            hit_file.write(json.dumps(record)+'\n');hit_file.flush()
            writer.writerow(row)
            if i%32==0:
                handle.flush()
            if (i+1)%256==0:
                print(i+1,len(chosen),flush=True)


if __name__=='__main__':
    main()
