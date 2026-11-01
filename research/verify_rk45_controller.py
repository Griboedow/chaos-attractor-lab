"""Check frozen RK45 bytes, all classifications and the displayed finite evidence."""
import argparse,json
from pathlib import Path
import numpy as np
from PIL import Image
from verify_rk23_controller import read,sha


def main(root,assets):
    js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    manifest=js(root/'manifest.json')
    for r in manifest['files']:
        p=root/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],p
    assert {r['path'] for r in manifest['files']}=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for r in js(root/'assets.json')['files']:assert sha(assets/r['path'])==r['sha256']==sha(root/r['archive_path'])
    scout=read(root/'oscillator/scout/screen.npz');follow=read(root/'oscillator/follow/follow.npz')
    assert len(scout['params'])==2051 and scout['flags'].sum()==1061 and (scout['counters'][:,0]==0).sum()==2050
    assert list(np.flatnonzero(scout['counters'][:,0]))==[802] and scout['counters'][802,0]==2
    assert np.array_equal(np.flatnonzero(scout['flags']),follow['ids'])
    flags=(follow['counters'][:,0]==0)&(follow['rates'][:,:,0].min(1)>.001)&(follow['period']==0)
    assert np.array_equal(flags,follow['flags'])
    good=follow['counters'][:,0]==0
    hyper_ids=follow['ids'][flags&(follow['rates'][:,:,1].min(1)>.001)]
    assert list(hyper_ids)==js(root/'oscillator/hyper-gallery/report.json')['ids']
    assert len(hyper_ids)==13 and good.sum()==1061 and flags.sum()==1051
    for ident in hyper_ids:
        cloud=read(root/f'oscillator/hyper-gallery/case-{ident}.npz')['states'];assert cloud.shape==(200000,3) and np.isfinite(cloud).all()
    assert np.max(np.abs(follow['rates'][good].sum(2)-follow['detlogs'][good]))<1e-10
    long=read(root/'oscillator/five-million/follow.npz')
    assert list(long['ids'])==[1890,335,448,508,2006,726] and np.all(long['flags'])
    assert np.max(np.abs(long['rates'].sum(2)-long['detlogs']))<1e-10
    assert np.all(long['rates'][list(long['ids']).index(335),:,1]>0)
    cohort=read(root/'oscillator/cohort/cohort.npz');assert len(cohort['ids'])==252
    for ident in [335,448,2006,1890]:
        for group,count in [('start',24),('step',6),('neighbor',32),('zero_control',1)]:
            mask=(cohort['ids']==ident)&(cohort['groups']==group);assert mask.sum()==count
            if group=='zero_control':assert cohort['counters'][mask,0]==2
        for digits in [70,90]:
            r=js(root/f'oscillator/precision/case-{ident}-{digits}.json');rates=np.array(r['rates'])
            assert rates.shape==(10,2) and np.allclose(rates.mean(0),r['means']) and np.allclose(rates.min(0),r['minima'])
            assert np.all(rates[:,0]>0)
            if ident==335:assert np.all(rates[:,1]>0)
    for row in js(root/'oscillator/cohort/report.json')['groups']:
        mask=(cohort['ids']==row['id'])&(cohort['groups']==row['group'])
        assert cohort['flags'][mask].sum()==row['positive'] and (cohort['period'][mask]>0).sum()==row['periodic']
        assert np.sum(cohort['flags'][mask]&(cohort['rates'][mask,:,1].min(1)>0))==row['two_positive']
    for r in js(root/'oscillator/upstream/report.json')['runs']:
        s=read(root/f'oscillator/upstream/case-{r["id"]}.npz')['states'];rates=read(root/f'oscillator/upstream/case-{r["id"]}-analysis.npz')['rates']
        assert s.shape==(400001,5) and np.isfinite(s).all() and np.all(np.diff(s[:,3])>0)
        assert np.allclose(rates.mean(0),r['means']) and np.allclose(rates.min(0),r['minima']) and np.all(rates[:,0]>0)
        if r['id']==335:assert np.all(rates[:,1]>0)
        assert r['max_one_step_residual']<1e-6 and r['max_window_determinant_identity_residual']<1e-10
    for ident in [335,448]:
        d=read(root/f'oscillator/approach/case-{ident}.npz');assert d['snapshots'].shape==(6,32768,3) and np.all(d['status']==0)
        assert np.all(np.isfinite(d['snapshots'])) and np.median(d['distances'][-1])<.003
    render=js(root/'oscillator/render/analysis.json');assert render['frame_point_counts']==[[1000000,1000000]]*72
    hist=read(root/'oscillator/render/portrait-histograms.npz')['histograms'];assert np.array_equal(hist.sum((1,2)),[8000000,8000000])
    with Image.open(root/'oscillator/render/rk45-oscillator-rotation.gif') as im:assert im.n_frames==72 and im.size==(1400,660)
    branch=js(root/'oscillator/branches/report.json');assert branch['dense_rosettes']['velocity_sign_alternates_every_step'] and branch['dense_rosettes']['points']==8000000
    for r in branch['runs']:assert r['reflection_max_residual']==0 and r['consecutive_state_max_residual']<1e-10 and sum(r['rejection_histogram'].values())==200000
    for r in js(root/'oscillator/phase/report.json')['runs']:
        d=read(root/f'oscillator/phase/case-{r["id"]}.npz');assert len(d['phase'])==300000 and r['max_complex_step_residual']<1e-12
        assert np.allclose(np.quantile(np.pi-np.abs(d['phase']),[0,.1,.5,.9,1]),r['distance_from_half_turn_quantiles'])
    for r in js(root/'oscillator/library-clouds/report.json')['runs']:
        d=read(root/f'oscillator/library-clouds/case-{r["id"]}.npz')
        for key in ['actual_to_reference','reference_to_actual','reference_halfcloud']:
            assert np.allclose(np.quantile(d[key],r['quantile_levels']),r[key+'_quantiles'])
    norm=js(root/'oscillator/error-norm-comparison/report.json')
    for r in norm['runs']:
        d=read(root/f'oscillator/error-norm-comparison/case-{r["id"]}.npz');assert np.allclose(d['rates'].mean(0),r['means']) and np.allclose(d['rates'].min(0),r['minima'])
    audit=js(root/'source/structure/analysis.json');assert audit['status']=='passed'
    assert sha(root/'source/structure/audit_rk45_structure.py')==audit['source_sha256']
    for r in audit['selected_cases']:assert r['gain_condition']==(r['id']!=1890) and r['sturm_roots']==0
    damping=read(root/'oscillator/damping/scan.npz');dr=js(root/'oscillator/damping/report.json')
    assert len(damping['ids'])==514 and np.array_equal(damping['rates'][128],damping['rates'][385])
    assert np.array_equal(damping['params'][128],scout['params'][335]) and np.array_equal(damping['params'][385],scout['params'][335])
    flags2=(damping['counters'][:,0]==0)&(damping['rates'][:,:,0].min(1)>.001)&(damping['period']==0)
    assert np.array_equal(flags2,damping['flags']) and np.array_equal(flags2&(damping['rates'][:,:,1].min(1)>0),damping['two_positive'])
    for row in dr['groups']:
        m=damping['groups']==row['group'];assert m.sum()==257
        assert damping['flags'][m].sum()==row['positive'] and damping['two_positive'][m].sum()==row['two_positive'] and (damping['period'][m]>0).sum()==row['periodic']
    good2=damping['counters'][:,0]==0;assert np.max(np.abs(damping['rates'][good2].sum(2)-damping['detlogs'][good2]))<1e-10
    result=dict(status='passed',manifest_files=len(manifest['files']),scout=2051,scout_complete=2050,follow=len(follow['ids']),follow_complete=int(good.sum()),follow_retained=int(flags.sum()),follow_two_positive=int(np.sum(flags&(follow['rates'][:,:,1].min(1)>.001))),cohort=252,volume_starts=65536,scope='Integrity and finite-result consistency, not an asymptotic chaos or priority proof.')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--assets',type=Path,required=True);a=p.parse_args();main(a.root,a.assets)
