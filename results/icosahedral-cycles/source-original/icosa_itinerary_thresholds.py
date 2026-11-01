"""Sensitivity of the finite saddle-neighborhood itinerary to its radius."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree

base=Path(__file__).parent;out=base/'icosa-itinerary-thresholds';out.mkdir(exist_ok=False)
foci=np.load(base/'icosa-itinerary-final/data.npz')['foci'];z=np.load(base/'icosa-long/177-dense.npy',mmap_mode='r');dist,idx=cKDTree(foci).query(z)
direction=foci/np.linalg.norm(foci,axis=1)[:,None];edge=np.abs(direction@direction.T-np.cos(np.pi/5))<1e-8;np.fill_diagonal(edge,False)
rows=[]
for radius in [.06,.08,.1,.12,.16,.2,.25,.3]:
    ids=idx[dist<radius];visits=ids[np.r_[True,ids[1:]!=ids[:-1]]];counts=np.zeros((30,30),int)
    np.add.at(counts,(visits[:-1],visits[1:]),1)
    rows.append({'radius':radius,'events':len(visits),'visited_foci':len(np.unique(visits)),'directed_edges':int((counts>0).sum()),'nongeometric_edges':int(((counts>0)&~edge).sum()),'nongeometric_transitions':int(counts[~edge].sum()),'fraction_near_focus':float((dist<radius).mean())})
(out/'report.json').write_text(json.dumps({'rows':rows,'limits':'Same finite orbit, several disjoint neighborhood radii. Missed shallow visits can create apparent long jumps. Agreement over a radius range is a robustness diagnostic, not an invariant graph proof.'},indent=2));print(json.dumps(rows),flush=True)
