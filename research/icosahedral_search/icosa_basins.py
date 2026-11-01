"""Fresh basins for the two retained geometries; all outcomes are archived."""
from pathlib import Path
import sys
from icosa_paths import ROOT,RESEARCH,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json,time
import numpy as np
from scipy.stats import qmc
from tetra_spectra import integrate

root=OUTPUT;d=np.load(root/'icosa-v1/input.npz');out=root/'icosa-basins';out.mkdir(exist_ok=False)
rng=np.random.default_rng(2026110148);ids=[];starts=[];labels=[]
for ident in [183,177]:
    cloud=np.load(root/f'icosa-long/{ident}-dense.npy',mmap_mode='r')
    volume=6*qmc.Sobol(3,scramble=True,seed=2026110148+ident).random_base2(4)-3
    local=cloud[np.linspace(0,len(cloud)-1,8,dtype=int)]+1e-4*rng.normal(size=(8,3))
    origin=1e-4*rng.normal(size=(4,3))
    for label,points in [('volume',volume),('local',local),('origin',origin)]:
        for z in points:ids.append(ident);starts.append(z);labels.append(label)
ids=np.array(ids);starts=np.array(starts);labels=np.array(labels);hs=np.full(len(ids),.0125)
design={'ids':ids.tolist(),'labels':labels.tolist(),'warm_time':1000,'observed_time':4000,'window_time':1000,'h':.0125,'positive_gate':'All four leading rates exceed 0.01; terminal speed exceeds 1e-5. Origin endpoint norm below 1e-6 is recorded separately.','seed':2026110148}
(out/'design.json').write_text(json.dumps(design,indent=2));np.savez_compressed(out/'input.npz',ids=ids,starts=starts,labels=labels)
begin=time.time();w,z,means,radius,end=integrate(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],starts,hs,1000,4000,1000,8000)
from tetra_spectra import field
speed=np.array([np.linalg.norm(field(e,d['centers'][i],d['axes'][i],d['sigma'][i],d['gamma'][i])[0]) for i,e in zip(ids,end)])
positive=(w[:,:,0].min(1)>.01)&(speed>1e-5);near=np.linalg.norm(end,axis=1)<1e-6
np.savez_compressed(out/'data.npz',ids=ids,labels=labels,starts=starts,h=hs,windows=w,tail=z,means=means,radius=radius,end=end,speed=speed)
rows=[]
for ident in [183,177]:
    for label in ['volume','local','origin']:
        mask=(ids==ident)&(labels==label)
        rows.append({'id':ident,'kind':label,'count':int(mask.sum()),'positive':int((positive&mask).sum()),'origin':int((near&mask).sum()),'other':int((mask&~positive&~near).sum())})
(out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2));print(json.dumps(rows),flush=True)
