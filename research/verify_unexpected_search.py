"""Verify frozen evidence and consequential numerical/analytical checks."""
from pathlib import Path
import ast,hashlib,json
import numpy as np
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/unexpected-search'

def main():
    manifest=json.loads((BASE/'manifest.json').read_text())
    for item in manifest['files']+json.loads((BASE/'reference-inputs.json').read_text()):
        p=ROOT/item['path'];data=p.read_bytes()
        assert len(data)==item['bytes'],p
        assert hashlib.sha256(data).hexdigest()==item['sha256'],p
    for p in (ROOT/'research/unexpected_search').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    for p in (ROOT/'research/unexpected_search').glob('*.py'):
        original=BASE/'source-original'/p.name
        if not original.exists():continue
        expected=original.read_text().replace('Path(__file__).parent','OUTPUT')
        expected=expected.replace("str(Path('D:/Codex/chaos-attractor-lab/research'))",'str(RESEARCH)')
        expected=expected.replace("'D:/Codex/chaos-attractor-lab/research'",'str(RESEARCH)')
        expected=expected.replace("Path('D:/Codex/chaos-attractor-lab')",'ROOT')
        assert p.read_text()=='from paths import OUTPUT, ROOT, RESEARCH\n'+expected,p
    pilot=np.load(BASE/'tetra-v1/data.npz');group=pilot['group']
    assert len(group)==12
    for a in group:
        assert np.allclose(a@a.T,np.eye(3)) and np.linalg.det(a)==1
        for b in group:assert any(np.array_equal(a@b,c) for c in group)
    assert int(((pilot['windows'].min(1)>.003)&(pilot['speed']>1e-4)).sum())==127
    long=np.load(BASE/'tetra-long/data.npz');passed=(long['windows'].min(1)>.003)&(long['speed']>1e-4)
    assert len(long['ids'])==36 and passed.sum()==26
    assert 896 in long['ids'][~passed]
    full=np.load(BASE/'tetra-spectra/data.npz');w=full['windows']
    assert w.shape==(12,10,3) and np.isfinite(w).all() and np.min(w[:,:,0])>0
    assert np.max(np.abs(w.sum(2)+3*pilot['gamma'][full['ids']][:,None]))<1e-8
    for row in json.loads((BASE/'tetra-dop853/summary.json').read_text()):
        assert np.min(np.array(row['windows'])[:,0])>0 and row['trace_error']<1e-10
    old=np.load(BASE/'anderson-v1/data.npz');assert len(old['status'])==4096 and np.all(old['status']==0)
    extended=np.load(BASE/'anderson-v2/data.npz');assert len(extended['positive_ids'])==132
    long=np.load(BASE/'anderson-long/data.npz');positive=(long['status']==1)&(long['windows'].min(1)>.005)
    assert positive.sum()==119 and (long['status']==0).sum()==9
    for row in json.loads((BASE/'anderson-precision/summary.json').read_text()):assert min(row['windows'])>0
    mp.iv.dps=60
    for row in json.loads((BASE/'anderson-basins/unique-root.json').read_text()):
        k=mp.iv.mpf(row['k_decimal']);c=mp.iv.mpf(row['c_decimal']);t=mp.iv.sqrt(1-1/k)
        b=k*t-mp.iv.log((1+t)/(1-t))/2
        margin=(c if float(row['c_decimal'])>0 else -c)-b
        assert margin.a>0
    for row in json.loads((BASE/'visuals/provenance.json').read_text()):
        if 'camera' not in row:continue
        camera=np.array(row['camera']);assert np.max(np.abs(camera@camera.T-np.eye(2)))<1e-12
        p=BASE/'visuals'/row['file'].replace('.png','-histogram.npz')
        h=np.load(p)['histogram'];assert int(h.sum())==row['points']==3200000
    report=json.loads((BASE/'comparison/report.json').read_text())
    assert report['known_regime_count']==108 and report['including_earlier_repo_candidates']==114
    assert json.loads((BASE/'comparison/control.json').read_text())['exact_affine_fingerprint_error']<1e-12
    print(json.dumps({'verified_files':len(manifest['files']),'reference_inputs':110,'tetra_scout_positive':127,'tetra_long_positive':26,'anderson_long_positive':119,'unique_scalar_root_cases':2,'result':'Evidence and controls pass. Historical novelty and infinite-time chaos are not certified.'}))

if __name__=='__main__':main()
