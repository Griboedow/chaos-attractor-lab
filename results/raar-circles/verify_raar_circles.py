"""Audit completeness, numerical claims, media counts and byte provenance."""
import argparse,json
from pathlib import Path
import numpy as np
from PIL import Image
from verify_rk23_controller import read,sha


def main(args):
    root=args.root;manifest=json.loads((root/'manifest.json').read_text());seen=[]
    for f in manifest['files']:
        p=root/f['path'];assert p.stat().st_size==f['bytes'] and sha(p)==f['sha256'],p;seen.append(f['path'])
    assert len(seen)==len(set(seen))==manifest['count']
    assert sum(f['bytes'] for f in manifest['files'])==manifest['bytes']
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(seen)|{'manifest.json'}
    s=read(root/'scout/screen.npz');f=read(root/'follow/follow.npz');long=read(root/'five-million/follow.npz')
    assert len(s['params'])==2051 and np.all(s['counters'][:,0]==0) and s['flags'].sum()==53
    assert np.array_equal(f['ids'],np.flatnonzero(s['flags'])) and np.all(f['counters'][:,0]==0)
    assert f['flags'].sum()==51 and np.sum(f['period']>0)==1
    assert f['period'][list(f['ids']).index(1540)]==2
    assert set(long['ids'])=={75,820,1051,1311,1483,1679} and np.all(long['flags'])
    for d in [f,long]:
        flags=(d['counters'][:,0]==0)&(d['period']==0)&(d['rates'][:,:,0].min(1)>.001)
        assert np.array_equal(flags,d['flags'])
        assert np.max(np.abs(d['rates'].sum(2)-d['detlogs']))<1e-9
    assert long['rates'][list(long['ids']).index(820),:,1].min()>0
    galleries=[json.loads((root/x/'report.json').read_text())['ids'] for x in ['gallery','remaining-gallery']]
    assert set(galleries[0]).isdisjoint(galleries[1]) and set(sum(galleries,[]))==set(f['ids'])
    for folder,ids in zip(['gallery','remaining-gallery'],galleries):
        for i in ids:
            z=read(root/f'{folder}/case-{i}.npz')['states'];assert z.shape==(200000,3) and np.isfinite(z).all()
    for i in long['ids']:
        for digits in [70,90]:
            d=json.loads((root/f'precision/case-{i}-{digits}.json').read_text());r=np.array(d['rates']);assert r.shape==(10,2) and r[:,0].min()>0
            if i==820:assert r[:,1].min()>0
    for i in [1051,1483,820,1679]:
        d=read(root/f'upstream/case-{i}.npz');assert d['states'].shape==(400001,3) and np.isfinite(d['states']).all() and d['rates'][:,0].min()>0
        assert np.max(np.abs(d['rates'].sum(1)-d['detlogs']))<1e-9
    r=read(root/'upstream/case-820.npz')['rates'];assert r[:,1].mean()>0 and r[:,1].min()<0
    co=read(root/'cohort/cohort.npz');assert len(co['ids'])==384 and np.all(co['counters'][:,0]==0)
    assert np.all(co['flags'][co['groups']=='start'])
    volume=0
    for i in [1051,1483,820,75]:
        d=read(root/f'approach/case-{i}.npz');assert d['starts'].shape==(32768,3) and np.all(d['status']==0)
        assert np.isfinite(d['snapshots']).all();volume+=len(d['starts'])
    sl=read(root/'relaxation/scan.npz');assert len(sl['ids'])==771 and np.all(sl['counters'][:,0]==0)
    for j,counts in enumerate([(151,0,106),(218,0,39),(232,82,25)]):
        ss=slice(j*257,(j+1)*257);assert tuple(map(int,[sl['flags'][ss].sum(),sl['two_positive'][ss].sum(),(sl['period'][ss]>0).sum()]))==counts
        i=[1051,1483,820][j];assert np.array_equal(sl['params'][j*257+128],s['params'][i])
    proof=json.loads((root/'fixed-proof/report.json').read_text());assert proof['status']=='passed'
    for case in proof['cases']:
        assert case['intersection']['disjoint'];assert case['number_of_fixed_points']==(2 if case['id']==1679 else 0)
        assert case['all_fixed_points_unstable'];assert case['parameters']==s['params'][case['id']].tolist()
    algebra=json.loads((root/'algebra/report.json').read_text());assert algebra['status']=='passed'
    assert {r['id']:r['linking_number'] for r in algebra['circle_links']}=={75:0,820:0,1051:0,1311:0,1483:0,1679:-1}
    render=json.loads((root/'render/analysis.json').read_text());assert np.array(render['frame_point_counts']).shape==(72,2) and np.all(np.array(render['frame_point_counts'])==1000000)
    h=read(root/'render/portrait-histograms.npz')['histograms'];assert np.all(h.sum((1,2))==8000000)
    for row in render['runs']:
        b=np.array(row['basis']);assert np.max(np.abs(b.T@b-np.eye(3)))<1e-12 and np.linalg.det(b)>0
        assert read(root/f'render/case-{row["id"]}-animation.npz')['states'].shape==(1000000,3)
    with Image.open(root/'render/raar-circles-rotation.gif') as im:assert im.n_frames==72
    for record in json.loads((root/'assets.json').read_text())['files']:
        assert sha(args.assets/record['path'])==sha(root/record['archive_path'])==record['sha256']
    result=dict(status='passed',manifest_sha256=sha(root/'manifest.json'),files=manifest['count'],bytes=manifest['bytes'],scout=2051,flags=53,follow_retained=51,gallery=53,five_million=6,independent_precision=12,original_library=4,cohort=384,volume=volume,slices=771,scope='Recomputed archive claims and integrity; interval certificates are separately executable. Not an asymptotic chaos or priority proof.')
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--assets',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
