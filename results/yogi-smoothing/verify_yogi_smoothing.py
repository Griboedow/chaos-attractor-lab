"""Verify the frozen smoothing packet and all56reported numerical outcomes."""
import argparse,json
from pathlib import Path
import numpy as np
from verify_optimizer_followups import check_manifest,arrays


def main(args):
    root=args.root;check_manifest(root);d=arrays(root/'control/histories.npz');a=json.loads((root/'control/analysis.json').read_text())
    assert a['completed']==a['declared']==56 and d['windows'].shape==(56,20,3)
    assert np.isfinite(d['windows']).all() and np.isfinite(d['tails']).all() and np.all(d['status']==0)
    for i,r in enumerate(a['runs']):
        z=d['tails'][i];w=d['windows'][i];tol=1e-9*(1+np.max(abs(z)))
        period=next((lag for lag in range(1,129) if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<tol),None)
        assert period==r['period_up_to128']
        assert np.allclose(w.mean(0),r['spectrum'],rtol=0,atol=1e-14) and np.allclose(w.min(0),r['window_minima'],rtol=0,atol=1e-14)
        assert np.max(abs(w.sum(1)-d['determinant_windows'][i]))<1e-7
        if r['grid_id']==125 and r['width']<=.032:assert np.all(w[:,:2]>0)
        if r['grid_id']==670 and r['width']==.01:assert period==32
        if r['grid_id']==125 and r['width']==1:assert period==4
    figures=arrays(root/'figure02/clouds.npz');assert figures['states'].shape==(8,200000,3) and np.isfinite(figures['states']).all()
    for r in json.loads((root/'figure02/render.json').read_text())['portraits']:
        b=np.array(r['basis']);assert np.max(abs(b.T@b-np.eye(3)))<1e-12 and np.linalg.det(b)>0 and r['clipped']==0
    print('PASS: frozen files,all56spectra,regular outcomes,determinant sums and8figure clouds')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);main(p.parse_args())
