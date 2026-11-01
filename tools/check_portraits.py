"""Check presentation provenance, retained sample counts, dates and local links."""
from pathlib import Path
import json,re,subprocess,sys
import numpy as np
from PIL import Image
from render_portraits import prepare,digest
from animate_portraits import load_views,projected,FRAMES,DURATION,WIDTH,HEIGHT

ROOT=Path(__file__).resolve().parents[1]
catalog=json.loads((ROOT/'provenance/portrait-catalog.json').read_text(encoding='utf-8'))
render=json.loads((ROOT/'docs/assets/portraits/render.json').read_text(encoding='utf-8'))
ledger=json.loads((ROOT/'provenance/render-dates.json').read_text(encoding='utf-8'))
assert render['catalog_sha256']==digest(ROOT/'provenance/portrait-catalog.json')
assert render['renderer_sha256']==digest(ROOT/'tools/render_portraits.py')
assert [r['id'] for r in render['rows']]==[r['id'] for r in catalog['rows']]
panels=0
for study,saved in zip(catalog['rows'],render['rows']):
    path=ROOT/saved['file'];assert path.is_file() and digest(path)==saved['sha256']
    with Image.open(path) as im:
        assert im.size==(1920,1080) and im.getpixel((0,0))[:3]==(255,255,255)
    row=ledger['images'][saved['file']]
    assert row['role']=='portrait' and row['source_portrait']==study['source_portrait']
    assert row['source_portrait'] in ledger['images']
    arrays,records=prepare(study)
    for hist,actual,expected in zip(arrays,records,saved['panels']):
        actual.pop('_xy',None)
        assert int(hist.sum())==expected['rendered_count']==expected['original_count']
        for key,value in actual.items():assert value==expected[key],(study['id'],key)
        if 'camera' in actual:
            camera=np.array(actual['camera']);assert np.max(abs(camera@camera.T-np.eye(2)))<1e-14
        panels+=1
animation=json.loads((ROOT/'docs/assets/portraits/animation.json').read_text(encoding='utf-8'))
assert animation['renderer_sha256']==digest(ROOT/'tools/animate_portraits.py')
assert animation['static_renderer_sha256']==digest(ROOT/'tools/render_portraits.py')
assert animation['catalog_sha256']==digest(ROOT/'provenance/portrait-catalog.json')
featured=[r for r in catalog['rows'] if r.get('featured')]
assert [r['id'] for r in animation['rows']]==[r['id'] for r in featured]
for study,saved in zip(featured,animation['rows']):
    path=ROOT/saved['file'];assert digest(path)==saved['sha256']
    assert path.stat().st_size==saved['bytes']<50_000_000
    assert saved['frames']==FRAMES and saved['duration_ms']==[DURATION]*FRAMES
    assert saved['seconds_per_revolution']==FRAMES*DURATION/1000
    with Image.open(path) as gif:
        assert gif.size==(WIDTH,HEIGHT) and gif.n_frames==FRAMES and gif.info['loop']==0
        checkpoints=[]
        for k in range(FRAMES):
            gif.seek(k);assert gif.info['duration']==DURATION
            if k in [0,36,72,108]:
                rgb=gif.convert('RGB');assert rgb.getpixel((0,0))==(255,255,255)
                checkpoints.append(rgb.tobytes())
        assert len(set(checkpoints))==4
    date_row=ledger['images'][saved['file']]
    assert date_row['role']=='portrait' and date_row['source_portrait']==study['source_portrait']
    views,records=load_views(study)
    for actual,expected in zip(records,saved['panels']):
        for key,value in actual.items():assert value==expected[key],(study['id'],key)
        C=np.array(actual['camera']);assert np.max(abs(C@C.T-np.eye(3)))<2e-12
        if 'basis' in actual:
            B=np.array(actual['basis']);assert np.max(abs(B.T@B-np.eye(3)))<2e-12
    assert saved['frame_point_counts']==[[len(v) for v in views]]*FRAMES
    for checkpoint in saved['checkpoint_projection_bounds']:
        q=views[checkpoint['panel']];xy=projected(q,2*np.pi*checkpoint['frame']/FRAMES)
        assert np.allclose(xy.min(0),checkpoint['minimum'],rtol=0,atol=1e-14)
        assert np.allclose(xy.max(0),checkpoint['maximum'],rtol=0,atol=1e-14)
        assert abs(xy).max()<records[checkpoint['panel']]['radius']
with np.load(ROOT/'results/portrait-presentation/icosa-animation.npz') as dense,np.load(ROOT/'results/icosahedral-flow/icosa-long/sample.npz') as frozen:
    for ident in [177,183]:
        assert dense[f'points_{ident}'].shape==(320000,3)
        assert np.array_equal(dense[f'points_{ident}'][::4],frozen['tail'][list(frozen['ids']).index(ident)])
documents=['README.md','docs/attractor-gallery.md','docs/research-index.md','docs/reproduction.md','tools/README-portraits.md','tools/README-render-dates.md','results/portrait-presentation/README.md']
bad=[];links=0
for name in documents:
    file=ROOT/name
    for raw in re.findall(r'\]\(([^)]+)\)',file.read_text(encoding='utf-8')):
        if raw.startswith(('https://','http://','mailto:','#')):continue
        link=(file.parent/raw.split('#')[0]).resolve();links+=1
        if not link.exists():bad.append((name,raw))
assert not bad,bad
main=(ROOT/'README.md').read_text(encoding='utf-8')
main_images=re.findall(r'^!\[[^\]]*\]\(([^)]+)\)',main,re.M)
assert main_images==[f"docs/assets/portraits/{row['id']}.gif" for row in featured]
for row in featured:assert f"[Static portrait](docs/assets/portraits/{row['id']}.png)" in main
assert len(re.findall(r'^### ',main,re.M))==len(featured)
gallery=(ROOT/'docs/attractor-gallery.md').read_text(encoding='utf-8')
gallery_images=re.findall(r'^!\[[^\]]*\]\(([^)]+)\)',gallery,re.M)
assert len(gallery_images)==len(catalog['rows']) and all(p.startswith('assets/portraits/') for p in gallery_images)
subprocess.run([sys.executable,str(ROOT/'tools/render_dates.py'),'check'],check=True)
print(json.dumps({'families':len(catalog['rows']),'panels':panels,'main_portraits':len(main_images),'local_links':links,'input_counts_and_hashes':'passed','coordinate_checks':'passed','animation_checks':'passed','icosahedral_sample_match':'bitwise'}))
