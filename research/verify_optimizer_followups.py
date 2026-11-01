"""Check frozen Yogi attraction outcomes and the AdaDelta cycle audit."""
import argparse,hashlib,json
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def arrays(p):
    with np.load(p) as d:return {k:d[k] for k in d.files}


def read(p):return json.loads(p.read_text())


def check_manifest(root):
    a=read(root/'manifest.json');paths=[]
    for r in a['files']:
        p=root/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'];paths.append(r['path'])
    assert len(paths)==len(set(paths))==a['count']
    assert set(paths)=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for r in read(root/'copy-provenance.json')['files']:assert sha(root/r['destination'])==r['sha256']
    for r in read(root/'assets.json')['files']:assert sha(root/r['archive_path'])==r['sha256']==sha(root.parents[1]/'docs/assets'/r['path'])


def main(args):
    root=args.results/'yogi-grid-attraction';check_manifest(root)
    co=arrays(root/'cohorts/histories.npz');ca=read(root/'cohorts/analysis.json');assert ca['completed']==ca['declared']==236
    assert co['windows'].shape==(236,4,3) and np.isfinite(co['tails']).all()
    for i,r in enumerate(ca['runs']):
        w=co['windows'][i];z=co['tails'][i];tol=1e-9*(1+np.max(abs(z)))
        period=next((lag for lag in range(1,129) if np.max(abs(z[lag:]-z[:-lag]))<tol),None)
        lead=bool(w[:,0].min()>.005 and np.ptp(z[:,0])>1e-4 and period is None)
        assert period==r['period_up_to128'] and lead==r['passed_leading']
        assert (lead and w[:,1].min()>.001)==r['passed_two']
        assert np.allclose(w.mean(0),r['spectrum'],rtol=0,atol=1e-14)
    for r in ca['summary']:
        selected=[x for x in ca['runs'] if x['id']==r['id'] and x['kind']==r['kind']]
        assert len(selected)==r['total'] and sum(x['passed_leading'] for x in selected)==r['passed_leading'] and sum(x['passed_two'] for x in selected)==r['passed_two']
    vol=arrays(root/'approach/volume-checks.npz');a=read(root/'approach/analysis.json')
    assert vol['windows'].shape==(256,4,3) and np.isfinite(vol['tails']).all()
    for i,r in enumerate(a['runs']):
        ident=r['grid_id'];d=arrays(root/f'approach/case-{ident}.npz')
        assert d['snapshots'].shape==(5,32768,3) and d['tails'].shape==(32768,64,3)
        assert np.array_equal(d['snapshots'][0],d['starts']) and np.array_equal(d['snapshots'][-1],d['tails'][:,-1])
        fixed=np.max(np.linalg.norm(np.diff(d['tails'],axis=1),axis=2),axis=1)<1e-10*(1+np.max(abs(d['tails']),axis=(1,2)))
        assert np.array_equal(fixed,d['fixed_tail']) and fixed.sum()==r['detected_fixed_tails']
        assert np.sum((d['endpoint_distances']>.02)&~fixed)==r['farther_than_002_and_nonfixed']==0
        assert np.array_equal(vol['starts'][128*i:128*(i+1)],d['starts'][:128])
        w=vol['windows'][128*i:128*(i+1)]
        assert np.all(w[:,:,0]>.005,axis=1).sum()==r['all_leading_windows_above005']
        assert np.all(w[:,:,1]>.001,axis=1).sum()==r['all_second_windows_above001']
    assert np.max(abs(vol['windows'].sum(2)-np.log(vol['parameters'][:,0,None])))<1e-9
    for r in read(root/'figure02/analysis.json')['panels']:
        assert r['points']==32768 and r['included']+r['excluded']==32768
        if not r['zoom']:assert r['excluded']==0
    root=args.results/'adadelta-cycle-audit';check_manifest(root)
    d=arrays(root/'replay/histories.npz');a=read(root/'replay/analysis.json');assert a['completed']==a['declared']==5 and a['steps_per_run']==100000000
    assert np.isfinite(d['snapshots']).all() and np.isfinite(d['tails']).all()
    for i,r in enumerate(a['runs']):
        assert np.array_equal(d['snapshots'][i,-1],r['final_state'])
        assert np.array_equal(d['snapshots'][i,list(d['times']).index(2)]-d['starts'][i],r['two_step_initial_residual'])
    assert .049<abs(d['snapshots'][0,-1,0])<.051
    assert abs(d['snapshots'][3,-1,0])>1.49
    for i in [1,2]:assert np.all(abs(d['snapshots'][i,:,0])==1)
    independent=read(root/'independent/analysis.json');actual=arrays(root/'independent/actual-torch.npz')['states'];assert actual.shape==(5000,3)
    assert np.array_equal(actual[-1],independent['torch']['state_after5000'])
    with localcontext() as ctx:
        ctx.prec=95
        left,right=independent['high_precision']
        for x,y in zip(left['samples'][-1]['state'],right['samples'][-1]['state']):assert abs(Decimal(x)-Decimal(y))<Decimal('1e-60')
    assert independent['torch']['endpoint_max_difference_from90digit_exact_rationals']<5e-12
    assert a['exact_algebra_cases']==512
    print('PASS: frozen files,236cohort outcomes,65536approach histories,256volume checks,AdaDelta replays and independent arithmetic',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);main(p.parse_args())
