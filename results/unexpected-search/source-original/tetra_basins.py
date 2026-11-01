"""Fresh volume starts, local perturbations and stable-origin controls."""
from pathlib import Path
import json,time
import numpy as np
from scipy.stats import qmc
from tetra_vortex import survey

root=Path(__file__).parent;d=np.load(root/'tetra-v1/data.npz');long=np.load(root/'tetra-spectra/data.npz');out=root/'tetra-basins';out.mkdir(exist_ok=False)
rng=np.random.default_rng(2026110137);params=[];starts=[];labels=[]
for ident in [421,565]:
    cloud=long['tail'][np.where((long['ids']==ident)&(long['h']==.0125))[0][0]]
    volume=4*qmc.Sobol(3,scramble=True,seed=2026110137+ident).random_base2(4)-2
    local=cloud[np.linspace(0,len(cloud)-1,8,dtype=int)]+1e-4*rng.normal(size=(8,3))
    origin=1e-4*rng.normal(size=(8,3))
    for label,points in [('volume',volume),('local',local),('origin',origin)]:
        for z in points:params.append(ident);starts.append(z);labels.append(label)
ids=np.array(params);starts=np.array(starts)
np.savez_compressed(out/'input.npz',ids=ids,starts=starts,labels=labels)
begin=time.time();w,z,end,speed,radius=survey(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],starts,80000,240000,.025,4000)
np.savez_compressed(out/'data.npz',ids=ids,starts=starts,labels=labels,windows=w,tail=z,end=end,speed=speed,radius=radius)
rows=[]
for i,label in zip([421,421,421,565,565,565],['volume','local','origin']*2):
    mask=(ids==i)&(np.array(labels)==label);positive=(w.min(1)>.003)&(speed>1e-4)
    rows.append({'id':i,'kind':label,'count':int(mask.sum()),'positive':int((positive&mask).sum()),'near_origin':int(((np.linalg.norm(end,axis=1)<1e-6)&mask).sum()),'other':int((mask&~positive&(np.linalg.norm(end,axis=1)>=1e-6)).sum())})
(out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2));print(json.dumps(rows),flush=True)
