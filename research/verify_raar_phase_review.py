"""Audit the phase diagnostic and the additional literature metadata."""
import json
from pathlib import Path
import numpy as np
from verify_raar_folding import sha

def main():
    repo=Path(__file__).resolve().parents[1]
    for name in ['raar-phase-classes','raar-priority-extension']:
        root=repo/'results'/name;manifest=json.loads((root/'manifest.json').read_text())['files'];assert {r['path'] for r in manifest}=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
        for r in manifest:p=root/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==r['bytes']
    root=repo/'results/raar-phase-classes';report=json.loads((root/'final/report.json').read_text());assert sha(repo/'results/raar-circles/render/case-1051-animation.npz')==report['source_sha256'];assert report['accuracy']==.1429225
    with np.load(root/'initial/phase-statistics.npz') as a,np.load(root/'final/phase-statistics.npz') as b:
        assert a.files==b.files and all(np.array_equal(a[k],b[k]) for k in a.files);assert np.array_equal(b['actual'],np.arange(600000,1000000)%7);assert np.array_equal(b['predicted'],b['class_distances'].argmin(0));assert np.array_equal(b['nearest_distances'],b['class_distances'].min(0));assert np.mean(b['actual']==b['predicted'])==report['accuracy'];assert b['histograms'].sum()==1000000
        confusion=np.bincount(b['actual']*7+b['predicted'],minlength=49).reshape(7,7);assert np.array_equal(confusion,report['confusion']);assert confusion.sum()==400000
        with np.load(root/'intermediate/phase-statistics.npz') as c:assert c.files==b.files and all(np.array_equal(c[k],b[k]) for k in c.files)
    assert sha(root/'final/raar-phase-classes.png')==sha(repo/'docs/assets/raar-phase-classes.png')
    print('PASS: exact manifests, all three numerical versions identical, all phase labels/counts/distances and held-out accuracy')
if __name__=='__main__':main()
