"""Verify retained/cache bytes, exact long-run replay and dense portrait data."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_raar_constellations import follow
from render_raar_constellations import orbit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from data_cache import data_path


ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/raar-constellations'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6);begin=time.time()
    inventory=json.loads((BASE/'inventory.json').read_text(encoding='utf-8'));checked=[]
    for row in inventory['files']:
        path=data_path(BASE/'runs'/row['run']/row['file'])
        assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256']
        checked.append(row['sha256'])
    print(json.dumps(dict(byte_checked=len(checked),seconds=time.time()-begin)),flush=True)
    report=json.loads((BASE/'runs/raar-constellation-long01/report.json').read_text(encoding='utf-8'))
    ids=report['ids'];params=[];points=[];counts=[]
    for ident in ids:
        with np.load(data_path(BASE/f'runs/raar-constellation-long01/case-{ident}.npz')) as f:
            params.append(f['params']);p=f['points'];points.append(np.pad(p,((0,12-len(p)),(0,0))));counts.append(len(p))
    params=np.array(params);points=np.array(points);counts=np.array(counts)
    rates,states,status,minima=follow(params,points,counts,100000,5000000)
    rows=[]
    for k,ident in enumerate(ids):
        with np.load(data_path(BASE/f'runs/raar-constellation-long01/case-{ident}.npz')) as f:
            for key,a in [('rates',rates[k]),('states',states[k]),('status',status[k]),('minima',minima[k])]:assert np.array_equal(a,f[key]),(ident,key)
        meta=json.loads((BASE/f'runs/raar-constellation-render01/case-{ident}.json').read_text(encoding='utf-8'))
        z,branches,mins=orbit(float(params[k][0]),points[k,:counts[k]],100000,2000000)
        with np.load(data_path(BASE/f'runs/raar-constellation-render01/case-{ident}-states.npz')) as f:
            assert np.array_equal(z,f['states']) and np.array_equal(branches,f['branches'])
        xy=((z-np.array(meta['center']))@np.array(meta['basis']))@np.array(meta['projection']).T
        h=np.histogram2d(xy[:,0],xy[:,1],bins=1600,range=[[-meta['halfwidth'],meta['halfwidth']]]*2)[0].T[::-1].astype(np.uint32)
        with np.load(data_path(BASE/f'runs/raar-constellation-render01/case-{ident}-histogram.npz')) as f:assert np.array_equal(h,f['histogram'])
        assert int(h.sum())==2000000 and np.array_equal(np.bincount(branches,minlength=counts[k]),meta['visits'])
        row=dict(id=ident,long_steps=5000000,dense_states=2000000,long_replay='bitwise equal',dense_replay='bitwise equal',histogram='bitwise equal')
        if ident in [1504,1725]:
            with np.load(data_path(BASE/f'runs/raar-constellation-branches01/case-{ident}-reduced.npz')) as f:
                row['reduced_fresh_run_bitwise_equal']=bool(np.array_equal(z[:len(f['states'])],f['states']))
                active=f['active'];subset=f['points'];assert np.array_equal(subset,points[k,active])
            inactive=np.array([j for j in range(counts[k]) if j not in active]);margin=np.inf
            for offset in range(0,len(z),100000):
                ds=np.sum((z[offset:offset+100000,None,:]-points[k,None,:counts[k],:])**2,axis=2)
                margin=min(margin,float(np.min(ds[:,inactive].min(1)-ds[:,active].min(1))))
            assert margin>0;row['observed_inactive_site_squared_distance_margin']=margin
        rows.append(row);print(json.dumps(dict(row,seconds=time.time()-begin)),flush=True)
    result=dict(status='passed',inventory_files=len(checked),inventory_sha256=sha(BASE/'inventory.json'),cases=rows,seconds=time.time()-begin,scope='Numerical replay and data accounting, not a proof of chaos, a basin or historical novelty.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
