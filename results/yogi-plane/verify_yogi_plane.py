"""Check hashes and recompute reported finite outcomes in the frozen Yogi grid."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def arrays(p):
    with np.load(p) as d:return {k:d[k] for k in d.files}


def main(args):
    root=args.root.resolve()
    def read(name):return json.loads((root/name).read_text())
    manifest=read('manifest.json');paths=[]
    for r in manifest['files']:
        p=root/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],p;paths.append(r['path'])
    assert len(paths)==len(set(paths))==manifest['count']
    assert set(paths)=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for r in read('copy-provenance.json')['files']:assert sha(root/r['destination'])==r['sha256']
    d=arrays(root/'plane/histories.npz');a=read('plane/analysis.json');protocol=read('plane/protocol.json')
    assert a['completed']==a['declared']==2048 and d['windows'].shape==(2048,4,3)
    assert d['tails'].shape==(2048,2048,3) and np.isfinite(d['tails']).all()
    assert np.array_equal(d['starts'],np.tile([.73,0,0],(2048,1)))
    assert np.max(abs(d['windows'].sum(2)-np.log(d['parameters'][:,0,None])))<1e-9
    for s,design in enumerate(protocol['slices']):
        ids=np.flatnonzero(d['slice_ids']==s)
        expected=np.array([[b,c,design['epsilon_normalized']] for c in design['beta2'] for b in design['beta1']])
        assert np.array_equal(expected,d['parameters'][ids])
    for offset in range(0,2048,128):
        block=arrays(root/f'plane/block-{offset:04d}.npz')
        for k in block:assert np.array_equal(block[k],d[k][offset:offset+128])
    periods=[]
    for z in d['tails']:
        periods.append(next((lag for lag in range(1,129) if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<1e-9*(1+np.max(abs(z)))),0))
    assert periods==a['periods']
    leading=(d['windows'][:,:,0].min(1)>.005)&(np.ptp(d['tails'][:,:,0],axis=1)>1e-4)&(np.array(periods)==0)
    two=leading&(d['windows'][:,:,1].min(1)>.001)
    assert np.array_equal(np.flatnonzero(leading),a['leading_ids']) and np.array_equal(np.flatnonzero(two),a['two_ids'])
    assert leading.sum()==356 and two.sum()==123
    assert periods.count(1)==1685 and periods.count(66)==1
    g=arrays(root/'plane-gallery/histories.npz');assert g['states'].shape==(24,100000,3) and len(set(g['ids']))==24
    assert set(g['ids']).issubset(set(a['leading_ids'])) and np.isfinite(g['states']).all()
    f=arrays(root/'plane-five-million/histories.npz');fa=read('plane-five-million/analysis.json')
    assert f['windows'].shape==(6,100,3) and fa['completed']==fa['declared']==6
    assert np.array_equal(f['parameters'],d['parameters'][f['ids']]) and np.array_equal(f['starts'],d['starts'][f['ids']])
    for i,r in enumerate(fa['runs']):
        w=f['windows'][i];assert r['id']==f['ids'][i]
        assert np.allclose(w.mean(0),r['spectrum'],rtol=0,atol=1e-14)
        assert np.allclose([w[:,1].min(),w[:,1].max()],r['second_window_range'],rtol=0,atol=1e-14)
    assert np.max(abs(f['windows'].sum(2)-np.log(f['parameters'][:,0,None])))<1e-9
    for ident in [125,798,28]:assert np.all(f['windows'][list(f['ids']).index(ident),:,:2]>0)
    assert f['windows'][list(f['ids']).index(63),:,0].min()<0
    precision=read('plane-precision/analysis.json');assert precision['completed']==precision['declared']==4
    for r in precision['runs']:
        z=arrays(root/'plane-precision'/f"case-{r['id']}-{r['digits']}digits.npz")
        assert np.array_equal(z['windows'],r['windows'])
        assert np.allclose(z['windows'].mean(0),r['spectrum_two'],atol=1e-15,rtol=0)
        assert z['windows'][:,0].min()>0
        if r['id']==125:assert z['windows'][:,1].min()>0
    u=arrays(root/'plane-upstream/histories.npz');ua=read('plane-upstream-analysis/analysis.json')
    assert u['states'].shape==(200000,6,3) and np.isfinite(u['states']).all()
    for r in ua['runs']:
        w=np.array(r['windows']);assert np.allclose(w.mean(0),r['spectrum'],atol=1e-14,rtol=0)
        assert r['autonomous_step_max_scaled_residual']<1e-12
        if r['id'] in [125,798,28]:assert w[:,:2].min()>0
        if r['id']==1506:assert r['tail_period_up_to128']==1
    render=read('plane-arc02/analysis.json');first=read('plane-arc/analysis.json')
    assert render['full_trajectory_sha256']==first['full_trajectory_sha256']
    assert arrays(root/'plane-arc02/portrait-histogram.npz')['histogram'].sum()==8000000
    animation=arrays(root/'plane-arc02/animation.npz');assert animation['states'].shape==(1000000,3)
    assert np.array_equal(animation['states'],arrays(root/'plane-arc/animation.npz')['states'])
    basis=np.array(render['basis']);assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12 and np.linalg.det(basis)>0
    assert render['animation']['frame_point_counts']==[1000000]*72
    assert min(render['late_return_residuals_1_to512'])>.09
    with Image.open(root/'plane-arc02/yogi-feathered-arc.gif') as im:
        assert im.n_frames==72 and im.size==(1200,630)
        for n in range(72):im.seek(n);assert im.info['duration']==100
    bound=read('scalar-bounds/analysis.json');assert bound['exact_controls']['cases']==1024 and bound['finite_stress']['cases']==384
    z=arrays(root/'scalar-bounds/histories.npz');assert np.isfinite(z['maxima']).all() and z['tails'].shape==(384,512,3)
    assert z['denominator_floor_errors'].max()==0
    for r in read('assets.json')['files']:
        p=root.parents[1]/'docs/assets'/r['path'];assert sha(p)==r['sha256']==sha(root/r['archive_path'])
    print('PASS: frozen bytes, all 2048 classifications, long/precision/library outcomes, figures and scalar controls',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);main(p.parse_args())
