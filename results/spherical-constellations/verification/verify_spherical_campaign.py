"""Check archive bytes, paired inputs, long-run replay and sphere portraits."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_spherical_constellations import follow
from render_spherical_constellations import orbit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from data_cache import data_path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/spherical-constellations'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    inventory=json.loads((BASE/'inventory.json').read_text(encoding='utf-8'))
    for row in inventory['files']:
        p=data_path(BASE/'runs'/row['run']/row['file'])
        assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256']
    sphere=data_path(BASE/'runs/spherical-constellation-scout01/screen.npz')
    circle=data_path(ROOT/'results/raar-constellations/runs/raar-constellation-scout01/screen.npz')
    with np.load(sphere) as a,np.load(circle) as b:
        for key in ['params','points','counts','kinds']:assert np.array_equal(a[key],b[key]),key
    rows=[]
    for run in ['spherical-constellation-long01','spherical-constellation-long02']:
        report=json.loads((BASE/'runs'/run/'report.json').read_text(encoding='utf-8'));ids=report['ids'];params=[];points=[];counts=[]
        for ident in ids:
            with np.load(data_path(BASE/'runs'/run/f'case-{ident}.npz')) as f:
                params.append(f['params']);p=f['points'];points.append(np.pad(p,((0,12-len(p)),(0,0))));counts.append(len(p))
        params=np.array(params);points=np.array(points);counts=np.array(counts)
        rates,states,status,minima,volume=follow(params,points,counts,100000,5000000)
        for j,ident in enumerate(ids):
            with np.load(data_path(BASE/'runs'/run/f'case-{ident}.npz')) as f:
                for key,a in [('rates',rates[j]),('states',states[j]),('status',status[j]),('minima',minima[j]),('volume_rates',volume[j])]:assert np.array_equal(a,f[key]),(ident,key)
            assert status[j,0]==0 and status[j,1]==5000000
            logdet_error=float(np.max(abs(rates[j].sum(1)-volume[j])));assert logdet_error<1e-8
            z=states[j];scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
            residuals=[float(np.max(np.linalg.norm(z[-2048+k:]-z[-2048:-k],axis=1))) for k in range(1,513)]
            periods=[k+1 for k,v in enumerate(residuals) if v<1e-8*scale]
            row=dict(id=ident,long_replay='bitwise equal',spectrum=rates[j].mean(0).tolist(),leading_windows=rates[j,:,0].tolist(),minimum_windows=rates[j].min(0).tolist(),qr_logdet_max_error=logdet_error,detected_period=min(periods) if periods else None,period_residuals=residuals)
            if ident in [1119,1790]:
                meta=json.loads((BASE/f'runs/spherical-constellation-render01/case-{ident}.json').read_text(encoding='utf-8'))
                z,branches,mins=orbit(float(params[j,0]),points[j,:counts[j]],100000,2000000)
                with np.load(data_path(BASE/f'runs/spherical-constellation-render01/case-{ident}-states.npz')) as f:assert np.array_equal(z,f['states']) and np.array_equal(branches,f['branches'])
                xy=((z-np.array(meta['center']))@np.array(meta['basis']))@np.array(meta['projection']).T
                h=np.histogram2d(xy[:,0],xy[:,1],bins=1600,range=[[-meta['halfwidth'],meta['halfwidth']]]*2)[0].T[::-1].astype(np.uint32)
                with np.load(data_path(BASE/f'runs/spherical-constellation-render01/case-{ident}-histogram.npz')) as f:assert np.array_equal(h,f['histogram'])
                assert int(h.sum())==len(z);row['dense_replay']='two million states, branch labels and histogram bitwise equal'
            rows.append(row);print(json.dumps(dict(id=ident,period=row['detected_period'],seconds=time.time()-begin)),flush=True)
    result=dict(status='passed',inventory_files=len(inventory['files']),inventory_sha256=sha(BASE/'inventory.json'),paired_input_check='params, points, counts and kinds bitwise identical to circle survey',cases=rows,seconds=time.time()-begin,scope='Numerical reproduction and data accounting. Period detection is numerical, not a certificate of exact initial-condition convergence.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
