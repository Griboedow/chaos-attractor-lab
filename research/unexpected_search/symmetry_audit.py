from paths import OUTPUT, ROOT, RESEARCH
"""Test a single dense trajectory against all twelve rotations without copying it."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree
from tetra_vortex import tetrahedral_group

root=OUTPUT;d=np.load(root/'tetra-v1/data.npz');group=tetrahedral_group();rng=np.random.default_rng(2026110142);rows=[]
for ident in [421,565]:
    z=np.load(root/f'tetra-render/dense-{ident}.npy',mmap_mode='r');half=len(z)//2
    a=z[rng.choice(half,30000,False)];b=z[half+rng.choice(half,6000,False)];rms=float(np.sqrt((z[::100]**2).sum(1).mean()));tree=cKDTree(a)
    errors=[float(np.sqrt(np.mean(tree.query(b@g.T)[0]**2))/rms) for g in group]
    baseline=float(np.sqrt(np.mean(tree.query(b)[0]**2))/rms)
    # Equal-volume nearest-direction cells; only a coarse occupation statistic.
    q=z[::32];sector=np.argmax(q@d['centers'][ident].T,axis=1);counts=np.bincount(sector,minlength=12)
    cov=np.cov(q.T);eig=np.linalg.eigvalsh(cov)
    rows.append({'id':ident,'group_errors':errors,'identity_error':baseline,'max_group_over_identity':max(errors)/baseline,'sector_counts':counts.tolist(),'mean_over_rms':float(np.linalg.norm(q.mean(0))/rms),'covariance_eigenvalues':eig.tolist(),'meaning':'Finite single-orbit support and occupation diagnostics, not an exact invariant-measure symmetry proof. Disjoint temporal halves; no group-replicated cloud.'})
(root/'tetra-render/symmetry.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows),flush=True)
