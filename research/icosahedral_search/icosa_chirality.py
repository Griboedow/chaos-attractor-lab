"""An orientation-sensitive degree-15 invariant, with finite-time block controls."""
from pathlib import Path
import sys
from icosa_paths import ROOT,RESEARCH,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json
import numpy as np

root=OUTPUT;d=np.load(root/'icosa-v1/input.npz');out=root/'icosa-chirality';out.mkdir(exist_ok=False)
normals=[]
for g in d['group']:
    if abs(np.trace(g)+1)<1e-10:
        ev,basis=np.linalg.eigh((g+g.T)/2);n=basis[:,-1]
        first=np.flatnonzero(np.abs(n)>1e-8)[0]
        if n[first]<0:n=-n
        normals.append(n)
normals=np.array(normals);assert len(normals)==15
def invariant(z):
    return np.prod((z/np.linalg.norm(z,axis=1)[:,None])@normals.T,axis=1)
rng=np.random.default_rng(2026110150);probes=rng.normal(size=(1000,3));p=invariant(probes)
error=max(float(np.max(np.abs(invariant(probes@g.T)-p))) for g in d['group'])
assert error<1e-16
assert np.max(np.abs(invariant(-probes)+p))<1e-20
rows=[]
for ident in [177,183,43]:
    z=np.load(root/f'icosa-long/{ident}-dense.npy',mmap_mode='r')
    values=np.concatenate([invariant(z[a:a+80000]) for a in range(0,len(z),80000)])
    blocks=np.array([a.mean() for a in np.array_split(values,16)]);rms=float(np.sqrt(np.mean(values**2)))
    row={'id':ident,'mean':float(values.mean()),'rms':rms,'normalized_mean':float(values.mean()/rms),'normalized_block_means':(blocks/rms).tolist(),'positive_fraction':float((values>0).mean()),'nonzero_points':int((np.linalg.norm(z,axis=1)>1e-12).sum())};rows.append(row)
    print(json.dumps(row),flush=True)
np.savez_compressed(out/'normals.npz',normals=normals,probes=probes)
(out/'report.json').write_text(json.dumps({'definition':'P(x) = product over the 15 unoriented twofold axes of dot(n_j, x/|x|). Axis signs are fixed by the first nonzero component being positive. It is invariant under the proper icosahedral group and odd under inversion. The sign convention is arbitrary.','rotation_error':error,'rows':rows,'limits':'Finite-time empirical moments. Blocks are not assumed independent; no calibrated p-value, asymptotic chirality theorem, or novelty probability is assigned.'},indent=2))
