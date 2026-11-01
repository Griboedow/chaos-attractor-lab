"""Verify frozen Adamax hashes and reconstruct the principal scientific summaries."""
import argparse,hashlib,json,re,subprocess
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(path.read_text(encoding='utf8'))


def normalized_step(states,p):
    x,m,u=states.T;b,c,e=p;M=b*m+(1-b)*x;U=np.maximum(c*u,np.abs(x)+e)
    return np.c_[x-M/U,M,U]


def main(args):
    root=Path(__file__).resolve().parents[1];base=root/'results/adamax-geometry';manifest=read(base/'manifest.json')
    assert manifest['count']==len(manifest['files'])
    actual={p.relative_to(base).as_posix() for p in base.rglob('*') if p.is_file() and p.name!='manifest.json'}
    assert actual=={r['path'] for r in manifest['files']}
    for row in manifest['files']:
        path=base/row['path'];assert path.stat().st_size==row['bytes'] and digest(path)==row['sha256'],path
    assets=read(base/'asset-manifest.json')['files']
    for row in assets:assert digest(root/row['path'])==row['sha256']
    if args.git_index:
        assert subprocess.check_output(['git','rev-parse','--show-object-format'],cwd=root,text=True).strip()=='sha1'
        entries=subprocess.check_output(['git','ls-files','--stage','-z','--','results/adamax-geometry','docs/assets/adamax-*'],cwd=root).split(b'\0')
        index={}
        for entry in entries:
            if entry:
                header,path=entry.split(b'\t',1);index[path.decode('utf8')]=header.split()[1].decode()
        files=[base/r['path'] for r in manifest['files']]+[base/'manifest.json']+[root/r['path'] for r in assets]
        for path in files:
            data=path.read_bytes();oid=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
            assert index[path.relative_to(root).as_posix()]==oid,path
    scout=np.load(base/'scout/scout.npz');flags=(scout['status']==0)&(scout['windows'].min(1)>.005)&(np.diff(scout['extrema'][:,:,0],axis=1)[:,0]>1e-4)
    assert len(flags)==4101 and int((scout['status']==0).sum())==4091 and int(flags.sum())==2169
    zero=np.all(scout['tails'][:,:,:2]==0,axis=(1,2));assert np.flatnonzero(flags&zero).tolist()==[3691,3979]
    assert int((flags&(np.ptp(scout['tails'][:,:,0],axis=1)<1e-8)).sum())==82
    gallery=np.load(base/'gallery/histories.npz');assert len(gallery['ids'])==len(set(gallery['ids']))==36
    long=np.load(base/'long-four-million/histories.npz');assert long['windows'].shape==(4,80,2) and np.all(long['windows'][:,:,0]>.005)
    report=read(base/'long-four-million/analysis.json')['runs']
    for i,row in enumerate(report):assert int(long['ids'][i])==row['id'] and np.allclose(long['windows'][i].mean(0),row['spectrum_two'],atol=1e-14)
    shorter=np.load(base/'long-million/histories.npz');i=list(shorter['ids']).index(1103);tail=shorter['tails'][i,-2048:]
    assert np.max(np.linalg.norm(tail[4:]-tail[:-4],axis=1))<1e-10 and shorter['windows'][i,-1,0]<0
    precision=read(base/'precision/analysis.json');assert precision['completed']==precision['declared']==8
    assert {(r['id'],r['digits']) for r in precision['runs']}=={(i,d) for i in [847,239,4099,4096] for d in [70,90]}
    assert all(len(r['leading_windows'])==10 and min(r['leading_windows'])>0 for r in precision['runs'])
    cohort=np.load(base/'cohorts/histories.npz');rows=read(base/'cohorts/analysis.json')['runs'];counts={}
    assert len(rows)==216
    for i,row in enumerate(rows):
        tail=cohort['tails'][i,-2048:];threshold=1e-10*(1+np.linalg.norm(np.ptp(tail,axis=0)))
        period=next((k for k in range(1,65) if np.max(np.linalg.norm(tail[k:]-tail[:-k],axis=1))<threshold),None)
        assert period==row['tail_period_up_to64']
        if np.all(cohort['windows'][i,:,0]>.005) and period is None:
            key=(row['id'],row['kind']);counts[key]=counts.get(key,0)+1
    for ident in [847,239,4099,4096]:
        assert counts[(ident,'cold_start')]==(8 if ident==239 else 22)
        assert counts[(ident,'parameter_neighbor')]==(14 if ident==239 else 32)
    upstream=np.load(base/'upstream-long/histories.npz');assert upstream['states'].shape==(4,200001,3)
    for z,p in zip(upstream['states'],upstream['parameters']):
        assert np.max(np.abs(normalized_step(z[:-1],p)-z[1:])/(1+np.abs(z[1:])))<1e-12
    assert np.all(np.load(base/'upstream-analysis/windows.npz')['windows'][:,:,0]>0)
    volume=np.load(base/'approach/volume-checks.npz');assert volume['windows'].shape==(256,4,2)
    assert np.all(volume['windows'][:128,:,0]>.005)
    assert int(np.all(volume['windows'][128:,:,0]>.005,axis=1).sum())==127
    for ident in [847,4099]:
        cloud=np.load(base/f'approach/case-{ident}.npz');assert cloud['snapshots'].shape==(6,32768,3) and len(cloud['starts'])==32768
        if ident==847:assert cloud['endpoint_distances'].max()<.003135
        else:
            far=np.flatnonzero(cloud['endpoint_distances']>.01);assert len(far)==436 and 87 in far
    outliers=np.load(base/'outliers/histories.npz');assert outliers['windows'].shape==(437,10,2) and np.all(outliers['windows'][:,:,0]<0)
    tails=outliers['tails'];assert np.max(np.linalg.norm(tails[:,68:]-tails[:,:-68],axis=2))<1e-10
    classified=read(base/'outlier-analysis/analysis.json');assert classified['closed_period_counts']=={'68':437}
    cycle=np.load(base/'outlier-analysis/cycle.npz');orbit=cycle['orbit'];assert orbit.shape==(68,3)
    assert np.max(np.abs(normalized_step(orbit,cycle['parameters'])-np.roll(orbit,-1,axis=0)))<1e-10
    assert min(np.linalg.norm(orbit[k]-orbit[0]) for k in range(1,68))>.00469
    assert max(abs(np.linalg.eigvals(cycle['monodromy'])))<.877
    symmetry=np.load(base/'cycle-symmetry/symmetry.npz');assert list(np.bincount(symmetry['labels']))==[221,216]
    assert np.minimum(symmetry['maximum_tail_distance_to_original'],symmetry['maximum_tail_distance_to_reflected']).max()<1e-10
    render=read(base/'render/analysis.json');hist=np.load(base/'render/portrait-histograms.npz')['histograms']
    assert np.array_equal(hist.sum(axis=(1,2)),[8000000,8000000])
    assert np.array(render['animation_frame_point_counts']).shape==(72,2) and np.all(np.array(render['animation_frame_point_counts'])==1000000)
    late=read(base/'render/late-period-audit.json')['runs']
    for row in render['runs']:
        basis=np.array(row['basis']);assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12 and np.linalg.det(basis)>0
        z=np.load(base/f"render/case-{row['id']}-animation.npz")['states'];assert z.shape==(1000000,3)
        assert np.max(np.abs(normalized_step(z[:-1],row['parameters'])-z[1:]))<1e-12
        tail=z[-4096:];returns=[float(np.max(np.linalg.norm(tail[k:]-tail[:-k],axis=1))) for k in range(1,513)]
        audit=next(r for r in late if r['id']==row['id']);assert np.allclose(returns,audit['return_residuals'],atol=1e-14) and min(returns)>.8
    with Image.open(root/'docs/assets/adamax-rotation.gif') as im:assert im.size==(1400,840) and im.n_frames==72
    for name,size in [('geometry',(2800,1680)),('approach',(2760,1270)),('sheets',(2100,980))]:
        with Image.open(root/f'docs/assets/adamax-{name}.png') as im:assert im.size==size
    h=np.load(base/'approach-final/histograms.npz');assert h['total'].shape==(12,430,430) and np.all(h['total'].sum(axis=(1,2))==32768)
    assert np.array_equal(h['regular'].sum(axis=(1,2)),[0]*6+[436]*6)
    sheets=read(base/'sheets/metadata.json');assert sum(sheets['counts'])+sheets['unknown_initial_age_count']==2000000
    assert sheets['maximum_observed_age']==4
    for path in [root/'README.md',*root.glob('docs/adamax-*.md')]:
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf8')):
            if '://' in target or target.startswith('#'):continue
            assert (path.parent/target.split('#')[0]).exists(),(path,target)
    print(json.dumps(dict(manifest_files=manifest['count'],asset_files=len(assets),git_index_checked=args.git_index,
        scout=4101,unique_gallery_cases=36,long_cases=4,precision_runs=8,cohort_runs=216,actual_library_cases=4,
        volume_checks=256,regular_outlier_continuations=437,regular_least_period=68,scientific_summary_checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--git-index',action='store_true');main(p.parse_args())
