"""Continue period-doubled orbits into the final chaotic band at identical coefficients."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.spatial import cKDTree
from adaptive import polar_ring
from band_continuation import dynamics
from orbits import refine


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,required=True)
    ap.add_argument('--branch',type=Path,required=True)
    ap.add_argument('--long-return',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--orders',nargs='+',type=int,default=[2,4,8,16])
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    m=json.loads((a.packet/'manifest.json').read_text())
    anchor=next(r for r in m['cases'] if r['case']==180844)
    p0=np.asarray(next(r for r in m['stable_cycles'] if r['case']==28756)['parameters'])
    p1=np.asarray(anchor['parameters'])
    with np.load(anchor['path']) as f:cloud=f['points']
    scale=np.linalg.norm(np.std(cloud,axis=0));tree=cKDTree(cloud)
    table={r['fraction']:r for r in json.loads(a.branch.read_text())['trials']}
    for order in a.orders:
        frac={2:.90,4:.94,8:.945,16:.95}[order]
        start=table[frac].get('attracting_cycle')
        if order==16:start=json.loads(a.long_return.read_text())['accepted_cycles'][0]
        guess=np.array(start['refined_start']);T=start['period'];f=frac;delta=.002;steps=[]
        while f<1-1e-14:
            dest=min(1.,f+delta)
            p=p0.copy();p[8]=(1-dest)*p0[8]+dest*p1[8]
            axis=int(np.argmax(abs(polar_ring(guess,p))))
            got=refine(axis,guess,T,scale,p,1)
            accepted=False
            if got is not None:
                proposed,period=got;d,sol=dynamics(proposed,period,p,False)
                accepted=(d['closure']/scale<3e-6 and abs(d['neutral'][0]-1)<.025 and abs(period/T-1)<.005)
            if not accepted:
                delta*=.5
                if delta<.00004:break
                continue
            f=dest;guess=proposed;T=period;steps.append(dict(fraction=f,e=float(p[8]),**d))
            delta=min(.003,delta*1.25)
        report=dict(order=order,starting_fraction=frac,reached_fraction=f,completed=f>=1-1e-14,steps=steps)
        if report['completed']:
            d,sol=dynamics(guess,T,p1,True)
            n=32768 if order<=4 else 65536
            points=sol.sol(np.linspace(0,T,n,endpoint=False))[:3].T
            gaps=[]
            for den in range(2,33):
                x=np.linspace(0,T,19,endpoint=False)
                gaps.append(np.max(np.linalg.norm(sol.sol(x)[:3]-sol.sol((x+T/den)%T)[:3],axis=0))/scale)
            if min(gaps)<1e-3 or d['closure']/scale>2e-5 or abs(d['neutral'][0]-1)>.025 or d['modulus']<1.01:
                report['completed']=False;report['final_rejection']=d
            else:
                fpath=a.out/f'orbit-{order}.npz';np.savez_compressed(fpath,points=points)
                report.update(final=d,refined_start=guess.tolist(),period=T,period_order=order,
                    minimum_phase_divisor_rms=float(min(gaps)),cloud_distance_p95=float(np.quantile(tree.query(points[::5])[0],.95)/scale),
                    reverse_distance_p95=float(np.quantile(cKDTree(points).query(cloud[::16])[0],.95)/scale),
                    cloud_relative=fpath.resolve().relative_to(a.out.resolve().parents[1]).as_posix())
        (a.out/f'order-{order}.json').write_text(json.dumps(report,indent=2)+'\n')
        print(order,report['reached_fraction'],report['completed'],report.get('final',{}).get('dominant'),flush=True)


if __name__=='__main__':
    main()
