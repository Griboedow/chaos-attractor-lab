"""Complete search replay, follow-up accounting and independent local checks."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from search_product_alphabet import step,screen,follow,checks
from follow_spherical_constellations import qr_step
from render_product_alphabet import orbit


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    sources={str(d):{p.name:sha(p) for p in sorted(d.iterdir()) if p.is_file()} for d in [args.scout,args.follow,args.long,args.dense]}
    for name in [Path(__file__).name,'search_product_alphabet.py','render_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.scout/'screen.npz') as f:
        inputs={k:f[k].copy() for k in ['params','levels','orders','points','counts','flags']}
        local=checks(inputs['params'],inputs['levels'],inputs['orders'],inputs['points'],inputs['counts'])
        values=screen(inputs['params'],inputs['levels'],inputs['orders'],5000,15000)
        for key,value in zip(['rates','tails','status','minimum_radius'],values):assert np.array_equal(value,f[key],equal_nan=True),key
        positive=(f['status'][:,0]==0)&(f['rates'].min(1)>.005);periods=np.zeros(len(positive),int)
        for i in np.flatnonzero(positive):
            z=f['tails'][i];scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
            for lag in range(1,65):
                if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<1e-10*scale:periods[i]=lag;break
        assert np.array_equal(positive,f['positive_flags']) and np.array_equal(periods,f['periods']) and np.array_equal(positive&(periods==0),inputs['flags'])
    summary=[]
    for folder in [args.follow,args.long]:
        report=json.loads((folder/'report.json').read_text(encoding='utf-8'));ids=report['ids']
        if folder==args.follow:assert ids==np.flatnonzero(inputs['flags']).tolist()
        for offset in range(0,len(ids),8):
            batch=ids[offset:offset+8];values=follow(inputs['params'][batch],inputs['levels'][batch],inputs['orders'][batch],report['warm'],report['observed'])
            for j,ident in enumerate(batch):
                row=next(r for r in report['cases'] if r['id']==ident)
                with np.load(folder/f'case-{ident}.npz') as f:
                    for key,value in zip(['rates','states','status','minimum_radius'],values):assert np.array_equal(value[j],f[key],equal_nan=True),(str(folder),ident,key)
                    assert f['status'].tolist()==row['status']==[0,report['observed']]
                    assert f['rates'].mean(0).tolist()==row['spectrum_qr_order']
                    assert float(f['rates'][:,0].min())==row['minimum_leading_window']
                    assert float(f['rates'][:,1].min())==row['minimum_second_window']
                    z=f['states'];scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
                    residuals=[float(np.max(np.linalg.norm(z[-2048+lag:]-z[-2048:-lag],axis=1))) for lag in range(1,513)]
                    assert residuals==row['return_residuals']
                    periods=[k+1 for k,v in enumerate(residuals) if v<1e-8*scale]
                    assert (min(periods) if periods else None)==row['detected_period']
                summary.append(dict(run=folder.name,id=ident,replay='All states, QR windows, status and radius bitwise equal.',period=row['detected_period'],minimum_leading=row['minimum_leading_window']))
            print(json.dumps(dict(run=folder.name,completed=min(offset+8,len(ids)),total=len(ids),seconds=time.time()-begin)),flush=True)
    dense=json.loads((args.dense/'report.json').read_text(encoding='utf-8'));rows=[]
    for ident in dense['ids']:
        beta,tau=inputs['params'][ident];levels=inputs['levels'][ident];order=int(inputs['orders'][ident])
        z,minimum=orbit(beta,tau,levels,order,dense['warm'],dense['count']);meta=json.loads((args.dense/f'case-{ident}.json').read_text(encoding='utf-8'))
        with np.load(args.dense/f'case-{ident}-states.npz') as f:assert np.array_equal(z,f['states'])
        xy=(z-np.array(meta['center']))@np.array(meta['basis'])@np.array(meta['projection']).T
        h=np.histogram2d(xy[:,0],xy[:,1],bins=1600,range=[[-meta['halfwidth'],meta['halfwidth']]]*2)[0].T[::-1].astype(np.uint32)
        with np.load(args.dense/f'case-{ident}-histogram.npz') as f:assert np.array_equal(h,f['histogram'])
        assert h.sum()==len(z) and minimum==meta['minimum_radius']
        x=z[-1].copy();Q=np.eye(3);L=np.eye(3);a=np.zeros(3);b=np.zeros(3);volume=0.
        for it in range(2048):
            x,J,r,_,_,flag=step(x,beta,tau,levels,order);assert not flag
            Q,logs,flag=qr_step(J,Q);assert not flag
            L,R=np.linalg.qr(J@L);d=np.diag(R);L*=np.where(d<0,-1.,1.);a+=logs;b+=np.log(abs(d));sign,logdet=np.linalg.slogdet(J);assert sign!=0;volume+=logdet
        row=dict(id=ident,dense_replay='All consecutive states and histogram bitwise equal.',lapack_mean_rate_error=float(np.max(abs(a-b))/2048),logdet_mean_error=float(abs(a.sum()-volume)/2048),parameter_hex=[float(v).hex() for v in [beta,tau]],level_hex=[[float(v).hex() for v in r[:order]] for r in levels])
        assert row['lapack_mean_rate_error']<1e-8 and row['logdet_mean_error']<1e-8
        rows.append(row);print(json.dumps(dict(id=ident,stage='dense and local QR',seconds=time.time()-begin)),flush=True)
    result=dict(status='passed',source_files=sources,scout_replay=f"All {len(inputs['params'])} cases, complete output arrays bitwise equal.",differential_checks=local,follow_replays=summary,dense=rows,seconds=time.time()-begin,scope='Finite numerical reproduction. Replay uses the production recurrence; independent checks concern algebraic formula, local derivatives and QR accounting. No interval or uniqueness certificate.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['scout','follow','long','dense','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
