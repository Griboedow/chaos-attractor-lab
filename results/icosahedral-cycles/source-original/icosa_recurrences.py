"""Short near recurrences after proper-rotation reduction; no orbit claim yet."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree

base=Path(__file__).parent;d=np.load(base/'icosa-v1/input.npz');out=base/'icosa-recurrences';out.mkdir(exist_ok=False)
group=d['group'];z=np.load(base/'icosa-long/177-dense.npy',mmap_mode='r')[::16];dt=.2
anchor=np.array([1.,np.sqrt(2)/3,np.sqrt(3)/7]);anchor/=np.linalg.norm(anchor)
labels=np.argmax(z@np.einsum('gji,j->gi',group,anchor).T,axis=1)
canonical=np.einsum('nij,nj->ni',group[labels],z);tree=cKDTree(canonical)
query_indices=np.arange(0,len(z),5);dist,neighbor=tree.query(canonical[query_indices],k=200)
candidates=[]
for i,dd,nn in zip(query_indices,dist,neighbor):
    for error,j in zip(dd,nn):
        lag=(j-i)*dt
        if 1.0<lag<30 and error<.03:
            G=group[labels[j]].T@group[labels[i]];which=int(np.argmin(np.linalg.norm(group-G,axis=(1,2))))
            discrepancy=float(np.linalg.norm(z[j]-group[which]@z[i]))
            assert abs(discrepancy-error)<1e-12
            candidates.append({'start_index':int(i),'end_index':int(j),'time_guess':float(lag),'error':discrepancy,'rotation':which,'score':float(error*np.exp(.13*lag)),'start':z[i].tolist()})
rows=[]
for row in sorted(candidates,key=lambda r:r['score']):
    if all(abs(row['start_index']-r['start_index'])>20 or abs(row['time_guess']-r['time_guess'])>.3 for r in rows):rows.append(row)
    if len(rows)==40:break
np.savez_compressed(out/'canonical.npz',points=canonical[::4],rotations=labels[::4],anchor=anchor)
(out/'candidates.json').write_text(json.dumps({'sample_dt':dt,'raw_qualifying_pairs':len(candidates),'selected':rows,'limits':'Near recurrences are starting guesses only. Canonicalization uses actual group rotations and is not applied to published orbit images.'},indent=2))
print(json.dumps({'pairs':len(candidates),'selected':len(rows),'first10':rows[:10]}),flush=True)
