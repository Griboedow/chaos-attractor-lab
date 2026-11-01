"""Verify frozen Yogi outcomes, structural controls, plots and byte-level provenance."""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path
import numpy as np
from PIL import Image


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(root,name):return json.loads((root/name).read_text())


def arrays(path):
    # NpzFile indexing decompresses again on every access; materialize once.
    with np.load(path) as data:return {key:data[key] for key in data.files}


def main(args):
    root=args.root.resolve();manifest=read(root,'manifest.json');paths=[]
    for row in manifest['files']:
        p=root/row['path'];assert p.is_file() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],p;paths.append(row['path'])
    assert len(paths)==len(set(paths))==manifest['count']
    assert set(paths)=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for row in read(root,'copy-provenance.json')['files']:assert sha(root/row['destination'])==row['sha256']
    d=arrays(root/'scout/scout.npz');a=read(root,'scout/analysis.json');assert len(d['parameters'])==a['completed']==a['declared']==8202
    assert np.all(d['status']==0)
    spans=np.ptp(d['tails'][:,:,0],axis=1)
    flags=np.all(d['windows']>.005,axis=1)&(spans>1e-4)&(d['periods']==0)
    assert np.array_equal(np.flatnonzero(flags),a['primary_ids']) and int(flags.sum())==268
    assert int(np.sum(np.max(abs(d['tails'][:,:,:2]),axis=(1,2))<1e-100))==a['zero_tail_count']==7901
    f=arrays(root/'follow-million/histories.npz');fa=read(root,'follow-million/analysis.json')
    assert np.array_equal(f['ids'],a['primary_ids']) and fa['completed']==fa['declared']==268
    follow=np.all(f['windows']>.005,axis=1)&(np.ptp(f['tails'][:,:,0],axis=1)>1e-4)&(f['periods']==0)
    assert np.array_equal(f['ids'][follow],fa['retained_ids']) and int(follow.sum())==266
    assert set(f['ids'][~follow])=={135,4454}
    g=arrays(root/'gallery/histories.npz');assert g['states'].shape==(24,100000,3) and len(set(g['ids']))==24 and np.isfinite(g['states']).all()
    for folder in ['spectra-million','spectra-five-million','cold-five-million']:
        data=arrays(root/f'{folder}/histories.npz');report=read(root,f'{folder}/analysis.json')
        assert report['completed']==report['declared']==len(data['ids'])
        residual=data['windows'].sum(2)-np.log(data['parameters'][:,0,None]);assert np.max(abs(residual))<1e-9
        for i,row in enumerate(report['runs']):
            assert row['id']==data['ids'][i] and np.allclose(data['windows'][i].mean(0),row['spectrum'],atol=1e-14,rtol=0)
            assert np.allclose([data['windows'][i,:,1].min(),data['windows'][i,:,1].max()],row['second_window_range'],atol=1e-14,rtol=0)
    w=arrays(root/'spectra-five-million/histories.npz');i=list(w['ids']).index(4214)
    assert w['windows'][i].shape==(100,3) and np.all(w['windows'][i,:,:2]>0)
    assert w['windows'][list(w['ids']).index(5174),:,1].min()<0
    checked=0
    for folder,n in [('precision',8),('precision-cold',2),('precision-square',2)]:
        a=read(root,f'{folder}/analysis.json');assert a['completed']==a['declared']==n
        for row in a['runs']:
            stem=f"case-{row['id']}-{row['digits']}digits";p=root/folder/(stem+'.npz');assert sha(p)==row['tail_sha256']
            data=arrays(p);assert np.array_equal(data['windows'],row['windows'])
            assert np.allclose(data['windows'].mean(0),row['spectrum_two'],atol=1e-15,rtol=0)
            if row['id']==4214:assert np.min(data['windows'][:,:2])>0
            if folder=='precision' and row['id']==7591:
                assert abs(row['spectrum_two'][0])<1e-8 and np.max(abs(data['tail'][:,:2]))<1e-50
            elif row['id'] in [44,7590,7591]:assert np.min(data['windows'][:,0])>0
            checked+=1
    co=arrays(root/'cohorts/histories.npz');ca=read(root,'cohorts/analysis.json');assert ca['completed']==ca['declared']==len(ca['runs'])==472
    for k,row in enumerate(ca['runs']):
        ww=co['windows'][k];assert np.allclose(ww.mean(0),row['spectrum'],atol=1e-14,rtol=0)
        assert np.allclose(ww.min(0),row['window_minima'],atol=1e-14,rtol=0)
        leading=bool(ww[:,0].min()>.005 and np.ptp(co['tails'][k,:,0])>1e-4 and row['period_up_to128'] is None)
        assert leading==row['passed_leading'] and (leading and ww[:,1].min()>.001)==row['passed_two']
        if row['period_up_to128'] is not None:
            assert row['period_up_to128']==1
            assert np.max(np.linalg.norm(np.diff(co['tails'][k],axis=0),axis=1))<1e-9*(1+np.max(abs(co['tails'][k])))
    for row in ca['summary']:
        items=[v for v in ca['runs'] if v['id']==row['id'] and v['kind']==row['kind']]
        assert len(items)==row['total'] and sum(v['passed_leading'] for v in items)==row['passed_leading'] and sum(v['passed_two'] for v in items)==row['passed_two']
    u=arrays(root/'upstream02/histories.npz');up=read(root,'upstream02/protocol.json');ua=read(root,'upstream-analysis02/analysis.json')
    assert u['states'].shape==(200000,9,3) and np.isfinite(u['states']).all()
    assert max(up['first20_max_scaled_state_error'])<5e-15 and up['optax_version']=='0.2.8'
    for row in ua['runs']:
        ww=np.array(row['windows']);assert ww.shape==(19,3) and np.allclose(ww.mean(0),row['spectrum'],atol=1e-14,rtol=0)
        assert row['autonomous_step_max_scaled_residual']<6e-16 and row['determinant_sum_max_residual']<1e-9
        if row['exact_fixed_line_linear_rates'] is not None:assert np.max(abs(np.array(row['spectrum'])-row['exact_fixed_line_linear_rates']))<1e-4
    actual=next(v for v in ua['runs'] if v['id']==4214);assert actual['spectrum'][1]>0 and np.min(np.array(actual['windows'])[:,1])<0
    st=read(root,'structure/analysis.json');z=d['starts'][7591];gap=Fraction.from_float(float(z[2]))-Fraction.from_float(float(z[0]))**2
    assert gap<0 and float(z[2]-z[0]*z[0])==0
    br=st['boundary_initialization'];assert gap.numerator==br['exact_margin_numerator'] and gap.denominator==br['exact_margin_denominator']
    assert st['inverse_controls']==512 and st['inverse_max_residual']<1e-12 and st['symmetry_max_residual']==0 and st['cold_bound_controls']==512
    for row in st['cold_bound_checks']:
        b,c,e=row['parameters'];delta=min(c,1-c);K=(1-b)/(np.sqrt(delta)*(1-b/np.sqrt(c)));B=max(abs(row['start'][0]),K/(1-b));V=max(row['start'][2],(2-c)*B*B)
        assert np.max(np.array(row['observed_maxima'])/np.array([B,B,V]))<=1+1e-12 and row['observed_maximum_step']<=K*(1+1e-12)
    render=read(root,'render/analysis.json');counts=np.array(render['animation_frame_point_counts']);assert counts.shape==(72,2) and np.all(counts==1000000)
    hist=arrays(root/'render/portrait-histograms.npz');assert np.all(hist['histograms'].sum((1,2))==8000000)
    im=Image.open(root/'render/yogi-rotation.gif');assert im.n_frames==72 and im.size==(1400,840)
    for k in range(im.n_frames):im.seek(k);assert im.info['duration']==100
    returns=[]
    for row in render['runs']:
        ident=row['id'];states=arrays(root/f'render/case-{ident}-animation.npz')['states'];assert states.shape==(1000000,3) and np.isfinite(states).all()
        basis=np.array(row['basis']);assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12 and np.linalg.det(basis)>0
        q=states[-4096:];res=np.array([np.max(np.linalg.norm(q[k:]-q[:-k],axis=1)) for k in range(1,513)])
        assert res.min()>1e-9*(1+np.max(abs(q)));returns.append(dict(id=ident,minimum_return_residual_up_to512=float(res.min())))
    inp=arrays(root/'render-input/endpoints.npz')
    for k,row in enumerate(read(root,'render-input/provenance.json')['records']):
        name='cold-five-million' if row['id']==7590 else 'spectra-five-million';p=root/name/'histories.npz';assert sha(p)==row['sha256'];data=arrays(p)
        for key in ['parameters','starts','ends']:assert np.array_equal(inp[key][k],data[key][row['index']])
    for row in read(root,'assets.json')['files']:
        p=root.parents[1]/'docs/assets'/row['path'];assert sha(p)==row['sha256']==sha(root/row['archive_path'])
    print(json.dumps(dict(status='PASS',manifest_files=len(paths),scouts=8202,followed=268,cohorts=472,precision_runs=checked,actual_library_cases=9,animation_frames=72,late_return_controls=returns)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);main(p.parse_args())
