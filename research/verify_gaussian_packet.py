"""Verify the frozen vortex and comparison packets, without certifying chaos."""
from pathlib import Path
import hashlib
import json
import numpy as np
from PIL import Image
from compare_attractor_shapes import fingerprint, distance
from compare_affine_shapes import whiten

ROOT=Path(__file__).resolve().parents[1]


def main():
    base=ROOT/'results/gaussian-vortex'
    comparison=ROOT/'results/attractor-novelty'
    manifest=json.loads((base/'manifest.json').read_text())
    listed=set()
    for entry in manifest['files']:
        path=ROOT/entry['path'];listed.add(entry['path'])
        assert path.stat().st_size==entry['bytes'],entry['path']
        with path.open('rb') as stream:
            assert hashlib.file_digest(stream,'sha256').hexdigest()==entry['sha256'],entry['path']
    for directory in [base,comparison]:
        for path in directory.rglob('*'):
            if path.is_file() and path!=base/'manifest.json':
                assert path.relative_to(ROOT).as_posix() in listed,path
    scout=np.load(base/'search-summary.npz');follow=np.load(base/'followup-summary.npz')
    assert scout['parameters'].shape==(1024,11)
    ids=np.where((scout['windows'].min(1)>.005)&(scout['speed']>1e-4))[0]
    assert len(ids)==45
    np.testing.assert_array_equal(ids,follow['ids'])
    np.testing.assert_array_equal(follow['parameters'],scout['parameters'][ids])
    assert follow['ids'][follow['windows'].min(1)<=.005].tolist()==[299,397,590,713]
    definitions=json.loads((base/'candidates.json').read_text())['candidates']
    for row in definitions:
        np.testing.assert_array_equal(row['parameters'],scout['parameters'][row['id']])
    summary=json.loads((base/'summary.json').read_text())
    long=np.load(base/'long-runs.npz')
    for filename,count in [('long-runs.npz',6),('fresh-starts.npz',16),('neighbors.npz',13)]:
        d=np.load(base/filename)
        for key in d.files:assert np.isfinite(d[key]).all(),(filename,key)
        assert d['windows'].shape==(count,int(d['observed_time']/d['window_time']),3)
        assert np.max(np.abs(d['windows'].sum(2)+3*d['parameters'][:,0,None]))<1e-7
        bounds=(np.linalg.norm(d['parameters'][:,8:11],axis=1)+np.exp(-.5)*np.sum(np.abs(d['parameters'][:,1:4])*d['parameters'][:,4:7],axis=1))/d['parameters'][:,0]
        assert np.all(d['max_radius']<bounds)
        if filename!='long-runs.npz':
            assert int((d['windows'][:,:,0].min(1)>.005).sum())==summary['cohorts'][filename[:-4]]['positive_all_windows']
    for i,row in enumerate(summary['long_runs']):
        assert row['candidate_id']==int(long['ids'][i])
        np.testing.assert_allclose(row['spectrum'],long['windows'][i].mean(0),rtol=0,atol=1e-14)
    assert long['windows'][:,:,0].min()>0
    checks=json.loads((base/'controls.json').read_text())
    assert max(checks['algebra_max_abs_residuals'].values())<1e-11
    adaptive=np.array(checks['independent_adaptive']['windows'])
    assert adaptive.shape==(12,3) and adaptive[:,0].min()>.02
    assert np.max(np.abs(adaptive.sum(1)+3*scout['parameters'][77,0]))<1e-8
    regular=json.loads((base/'regular-controls.json').read_text())['outcomes']
    assert [r['fresh_start_index'] for r in regular]==[4,6,8,12]
    assert all(r['shortest_tested_close_return']['max_return_error']<1e-7 for r in regular)
    render=json.loads((base/'render.json').read_text())
    p=np.array(render['projection_rows']);b=np.array(render['basis_columns'])
    np.testing.assert_allclose(p@p.T,np.eye(2),atol=1e-14,rtol=0)
    np.testing.assert_allclose(b.T@b,np.eye(3),atol=1e-14,rtol=0)
    assert np.linalg.det(b)>0
    h=np.load(base/'render-histogram.npz')['counts']
    assert h.shape==(2400,2400) and h.sum()==4800000 and render['clipped_points']==0
    assert np.all(h==np.floor(h)) and h.min()==0
    with Image.open(ROOT/'docs/assets/gaussian-vortex-triskelion.png') as picture:assert picture.size==(2400,2400)
    with Image.open(ROOT/'docs/assets/gaussian-vortex-triskelion.gif') as movie:assert movie.n_frames==72
    render_run=np.load(base/'render-run.npz')
    assert render_run['windows'].shape==(30,3) and render_run['windows'][:,0].min()>.019
    np.testing.assert_array_equal(render_run['windows'].mean(0),summary['render_spectrum'])
    replay=json.loads((base/'replay-checks.json').read_text())
    assert replay['render_metadata_identical'] and replay['sobol_parameter_array_identical'] and replay['portable_catalog_lorenz_identical']
    for name in ['gaussian-vortex-triskelion.png','gaussian-vortex-triskelion.gif']:
        assert hashlib.sha256((ROOT/'docs/assets'/name).read_bytes()).hexdigest()==replay[name]['identical_sha256']
    refs={path.stem:np.load(path)['points'] for path in sorted((comparison/'reference/trajectories').glob('*.npz'))}
    assert len(refs)==108
    for name,z in refs.items():
        assert z.ndim==2 and z.shape[1]==3 and len(z)>=8000 and np.isfinite(z).all(),name
    protocol=json.loads((comparison/'reference/protocol.json').read_text())
    assert len(protocol['selected'])==89 and len(protocol['excluded'])==46
    records=json.loads((comparison/'reference/nearby-runs.json').read_text())
    assert [r['name'] for r in records if r.get('success') is False]==['Halvorsen-a1.1','Halvorsen-a1.2']
    prior=json.loads((comparison/'reference/kicked-top-runs.json').read_text())
    with (ROOT/prior['source_path']).open('rb') as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==prior['source_sha256']
    old=np.load(ROOT/'results/twist-geometry/long-runs.npz')
    candidates={f'Twist-{i}':z for i,z in zip(old['ids'],old['tail'])}
    candidates.update({f'Vortex-{i}':z for i,z in zip(long['ids'][1::2],long['tail'][1::2])})
    for filename,transform in [('summary.json',lambda z:z),('affine-summary.json',whiten)]:
        report=json.loads((comparison/'comparison'/filename).read_text())
        assert report['reference_count']==108 and set(report['candidates'])==set(candidates)
        descriptors={n:fingerprint(transform(z)) for n,z in refs.items()}
        for name,z in candidates.items():
            fp=fingerprint(transform(z))
            ranks=sorted((float(distance(fp,d)),n) for n,d in descriptors.items())[:5]
            expected=report['candidates'][name]['nearest']
            assert [n for d,n in ranks]==[r['reference'] for r in expected]
            np.testing.assert_allclose([d for d,n in ranks],[r['fingerprint_distance'] for r in expected],rtol=0,atol=1e-12)
    a=refs['Lorenz'];linear=np.array([[2.,.7,-.2],[0.,.4,.1],[.3,0.,1.1]])
    assert distance(fingerprint(whiten(a)),fingerprint(whiten(a@linear+[3,-7,11])))<1e-12
    print(json.dumps({'verified_files':len(listed),'scout_rows':1024,'followup_rows':45,
                      'selected_full_spectrum_runs':6,'fresh_starts':16,'parameter_controls':13,
                      'reference_clouds':108,'comparison_candidates':6,'rendered_states':4800000,
                      'scope':'Integrity and numerical consistency, not a chaos or novelty certificate.'}))


if __name__=='__main__':
    main()
