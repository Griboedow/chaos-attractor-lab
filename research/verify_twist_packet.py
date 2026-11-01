"""Verify the frozen geometry packet and cross-check its numerical summaries."""
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image


def main():
    root=Path(__file__).resolve().parents[1]
    base=root/'results/twist-geometry'
    manifest=json.loads((base/'manifest.json').read_text())
    for entry in manifest['files']:
        path=root/entry['path']
        assert path.stat().st_size==entry['bytes'],entry['path']
        with path.open('rb') as stream:
            assert hashlib.file_digest(stream,'sha256').hexdigest()==entry['sha256'],entry['path']
    definitions=json.loads((base/'candidates.json').read_text())['candidates']
    search=np.load(base/'search-summary.npz')
    assert search['parameters'].shape==(4096,10)
    for row in definitions:
        np.testing.assert_array_equal(row['parameters'],search['parameters'][row['id']])
    total=0
    for filename,count in [('long-runs.npz',3),('basins.npz',96),('neighbors.npz',66)]:
        data=np.load(base/filename)
        assert len(data['ids'])==count
        assert data['windows'].shape==(count,int(data['observed'])//int(data['window']),3)
        for key in data.files:
            assert np.isfinite(data[key]).all(),(filename,key)
        assert np.max(data['radii'])<=1+1e-12
        assert np.max(np.linalg.norm(data['tail'],axis=2))<=1+1e-12
        assert np.max(np.abs(data['windows'].sum(2)-3*np.log(data['parameters'][:,0,None])))<1e-9
        total+=count
    long=np.load(base/'long-runs.npz');basins=np.load(base/'basins.npz');neighbors=np.load(base/'neighbors.npz')
    summary=json.loads((base/'summary.json').read_text())
    for row in summary['candidates']:
        i=list(long['ids']).index(row['id'])
        np.testing.assert_allclose(row['spectrum'],long['windows'][i].mean(0),rtol=0,atol=1e-14)
        w=basins['windows'][basins['ids']==row['id'],:,0]
        assert row['basin_positive_all_windows']==int((w.min(1)>0).sum())
        assert len(row['regular_basin_outcomes'])==len(w)-row['basin_positive_all_windows']
    for row in summary['neighbors']:
        w=neighbors['windows'][neighbors['ids']==row['id'],:,0]
        assert row['positive_all_windows']==int((w.min(1)>0).sum())
    checks=json.loads((base/'independent-checks.json').read_text())
    assert len(checks['precision'])==4
    assert all(min(row['windows'])>0 for row in checks['precision'])
    renders=json.loads((base/'render.json').read_text())
    for row in renders:
        projection=np.array(row['projection_after_pca']);basis=np.array(row['pca_basis_columns'])
        np.testing.assert_allclose(projection@projection.T,np.eye(2),rtol=0,atol=1e-14)
        np.testing.assert_allclose(basis.T@basis,np.eye(3),rtol=0,atol=1e-14)
        assert np.linalg.det(basis)>0
        h=np.load(base/f'render-{row["id"]}.npz')['histogram']
        assert int(h.sum())==row['observed'] and row['clipped_spatial_points']==0
        with Image.open(root/'docs/assets'/row['image']) as picture:
            assert picture.size==(row['resolution'],row['resolution'])
        with Image.open(root/'docs/assets'/f'twist-{row["label"]}.gif') as movie:
            assert movie.n_frames==row['animation']['frames']
    print(json.dumps({'verified_files':len(manifest['files']),'full_spectrum_runs':total,
                      'scout_rows':4096,'arbitrary_precision_runs':4,'rendered_points':16000000}))


if __name__=='__main__':
    main()
