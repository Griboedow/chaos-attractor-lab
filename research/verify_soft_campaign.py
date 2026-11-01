"""Hash accounting, long-run replay and independent QR/volume checks."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_soft_constellations import follow
from search_soft_constellations import step
from follow_spherical_constellations import qr_step
from render_soft_constellations import orbit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from data_cache import data_path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/soft-constellations'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    inventory=json.loads((BASE/'inventory.json').read_text(encoding='utf-8'))
    for row in inventory['files']:
        p=data_path(BASE/'runs'/row['run']/row['file']);assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256']
    with np.load(data_path(BASE/'runs/soft-constellation-scout01/screen.npz')) as soft,np.load(data_path(ROOT/'results/raar-constellations/runs/raar-constellation-scout01/screen.npz')) as hard:
        bases=[124,611,896,1504,1725,1873]
        for j,base in enumerate(bases):
            assert soft['params'][j*33,0]==hard['params'][base,0] and soft['params'][j*33,1]==0
            for name in ['points','counts','kinds']:assert np.array_equal(soft[name][j*33],hard[name][base])
            assert np.array_equal(soft['rates'][j*33],hard['rates'][base])
    folder=BASE/'runs/soft-constellation-long01';report=json.loads((folder/'report.json').read_text(encoding='utf-8'));ids=report['ids'];params=[];points=[];counts=[]
    for ident in ids:
        with np.load(data_path(folder/f'case-{ident}.npz')) as f:
            params.append(f['params']);p=f['points'];points.append(np.pad(p,((0,12-len(p)),(0,0))));counts.append(len(p))
    params=np.array(params);points=np.array(points);counts=np.array(counts)
    rates,states,status,minima,stats=follow(params,points,counts,100000,5000000);rows=[]
    for j,ident in enumerate(ids):
        with np.load(data_path(folder/f'case-{ident}.npz')) as f:
            for key,a in [('rates',rates[j]),('states',states[j]),('status',status[j]),('minima',minima[j]),('soft_stats',stats[j])]:assert np.array_equal(a,f[key]),(ident,key)
        assert status[j,0]==0 and status[j,1]==5000000
        row=dict(id=ident,long_replay='bitwise equal',spectrum=rates[j].mean(0).tolist(),minimum_windows=rates[j].min(0).tolist())
        # An independent LAPACK QR and determinant accumulation on one orbit.
        beta,tau=params[j];p=points[j,:counts[j]];x=states[j,0].copy();Q=np.eye(3);L=np.eye(3);s=np.zeros(3);t=np.zeros(3);volume=0.;max_step_error=0.
        for it in range(2048):
            x,J,r,peak,entropy,flag=step(x,beta,p,tau);assert not flag
            Q,a,flag=qr_step(J,Q);assert not flag
            L,R=np.linalg.qr(J@L);d=np.diag(R);L*=np.where(d<0,-1.,1.);b=np.log(abs(d))
            s+=a;t+=b;sign,logdet=np.linalg.slogdet(J);assert sign!=0;volume+=logdet;max_step_error=max(max_step_error,float(abs(a.sum()-logdet)))
        row['lapack_spectrum_max_error']=float(abs(s-t).max()/2048);row['mean_logdet_error']=float(abs(s.sum()-volume)/2048);row['maximum_step_logdet_error']=max_step_error
        assert row['lapack_spectrum_max_error']<1e-8 and row['mean_logdet_error']<1e-8
        if ident!=27:
            meta=json.loads((BASE/f'runs/soft-constellation-render01/case-{ident}.json').read_text(encoding='utf-8'))
            z,min_radius=orbit(beta,tau,p,100000,2000000)
            with np.load(data_path(BASE/f'runs/soft-constellation-render01/case-{ident}-states.npz')) as f:assert np.array_equal(z,f['states'])
            xy=(z-np.array(meta['center']))@np.array(meta['basis'])@np.array(meta['projection']).T
            h=np.histogram2d(xy[:,0],xy[:,1],bins=1600,range=[[-meta['halfwidth'],meta['halfwidth']]]*2)[0].T[::-1].astype(np.uint32)
            with np.load(data_path(BASE/f'runs/soft-constellation-render01/case-{ident}-histogram.npz')) as f:assert np.array_equal(h,f['histogram'])
            assert h.sum()==len(z);row['dense_replay']='two million consecutive states and histogram bitwise equal'
        rows.append(row);print(json.dumps(dict(id=ident,seconds=time.time()-begin)),flush=True)
    result=dict(status='passed',inventory_files=len(inventory['files']),inventory_sha256=sha(BASE/'inventory.json'),hard_control='All six exact input sets and ten leading-rate windows match the original hard-projection scout bitwise.',cases=rows,seconds=time.time()-begin,scope='Numerical reproduction and local derivative accounting, not interval certification or historical priority.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
