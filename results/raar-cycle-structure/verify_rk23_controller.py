"""Check frozen hashes and reported classifications against retained controller data."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from PIL import Image

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def read(path):
    with np.load(path) as f:return {k:f[k] for k in f.files}

def main(root,assets):
    manifest=json.loads((root/'manifest.json').read_text())
    for r in manifest['files']:
        p=root/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],p
    for r in json.loads((root/'assets.json').read_text())['files']:assert sha(assets/r['path'])==r['sha256']==sha(root/r['archive_path'])
    scout=read(root/'oscillator/scout/screen.npz');follow=read(root/'oscillator/follow/follow.npz')
    assert len(scout['params'])==1027 and scout['flags'].sum()==396
    assert np.array_equal(np.flatnonzero(scout['flags']),follow['ids'])
    flags=(follow['counters'][:,0]==0)&(follow['rates'][:,:,0].min(1)>.001)&(follow['period']==0)
    assert np.array_equal(flags,follow['flags']) and flags.sum()==393
    assert np.sum(flags&(follow['rates'][:,:,1].min(1)>.001))==3
    assert np.max(np.abs(follow['rates'].sum(2)-follow['detlogs']))<1e-10
    for name,ids in [('five-million',[381,354,883]),('hyper-five-million',[37,526,689])]:
        d=read(root/'oscillator'/name/'follow.npz');assert list(d['ids'])==ids
        assert np.all(d['counters'][:,0]==0) and np.max(np.abs(d['rates'].sum(2)-d['detlogs']))<1e-10
        assert np.all(d['rates'][:,:,0].min(1)>.001)
    d=read(root/'oscillator/cohort/cohort.npz');assert len(d['ids'])==189
    expected={381:[24,6,16],354:[24,6,32],883:[24,6,24]}
    for ident in expected:
        for group,count in zip(['start','step','neighbor'],expected[ident]):
            mask=(d['ids']==ident)&(d['groups']==group);assert d['flags'][mask].sum()==count
        mask=(d['ids']==ident)&(d['groups']=='zero_control');assert np.all(d['counters'][mask,0]==2)
    for folder,ids in [('precision',[381,354,883]),('hyper-precision',[37,526,689])]:
        for ident in ids:
            for digits in [70,90]:
                r=json.loads((root/'oscillator'/folder/f'case-{ident}-{digits}.json').read_text());rates=np.array(r['rates'])
                assert rates.shape==(10,2) and np.allclose(rates.mean(0),r['means']) and np.allclose(rates.min(0),r['minima'])
                assert rates[:,0].min()>0
    for folder in ['upstream','hyper-upstream']:
        for r in json.loads((root/'oscillator'/folder/'report.json').read_text())['runs']:
            s=read(root/'oscillator'/folder/f'case-{r["id"]}.npz')['states'];rate=read(root/'oscillator'/folder/f'case-{r["id"]}-analysis.npz')['rates']
            assert s.shape==(400001,5) and np.isfinite(s).all() and np.all(np.diff(s[:,3])>0)
            assert np.allclose(rate.mean(0),r['means']) and np.allclose(rate.min(0),r['minima'])
            assert r['max_one_step_residual']<1e-6 and r['max_window_determinant_identity_residual']<1e-10
    for ident in [381,883]:
        d=read(root/'oscillator/approach'/f'case-{ident}.npz');assert d['snapshots'].shape==(6,32768,3) and np.all(d['status']==0)
        assert np.all(np.isfinite(d['snapshots'])) and np.median(d['distances'][-1])<.003
    r=json.loads((root/'oscillator/render/analysis.json').read_text());assert r['frame_point_counts']==[[1000000,1000000]]*72
    hist=read(root/'oscillator/render/portrait-histograms.npz')['histograms'];assert np.array_equal(hist.sum((1,2)),[8000000,8000000])
    with Image.open(root/'oscillator/render/rk23-oscillator-rotation.gif') as im:assert im.n_frames==72 and im.size==(1400,660)
    pieces=[read(p) for p in sorted((root/'pi/follow').glob('follow-*.npz'))];joined={k:np.concatenate([d[k] for d in pieces]) for k in pieces[0]}
    pi=read(root/'pi/scout/screen.npz');assert len(pi['params'])==8198 and pi['flags'].sum()==1089
    assert np.array_equal(joined['ids'],np.flatnonzero(pi['flags'])) and joined['flags'].sum()==1079
    complexscout=read(root/'complex/scout/screen.npz');complexfollow=read(root/'complex/follow/follow.npz')
    assert len(complexscout['params'])==1027 and np.array_equal(complexfollow['ids'],np.flatnonzero(complexscout['flags'])) and complexfollow['flags'].sum()==35
    print(json.dumps(dict(status='passed',manifest_files=len(manifest['files']),oscillator_screen=1027,oscillator_follow=396,oscillator_retained=393,cohort=189,volume_starts=65536,scope='Integrity and finite-result consistency, not an asymptotic chaos or priority proof.')),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--assets',type=Path,required=True);a=p.parse_args();main(a.root,a.assets)
