from pathlib import Path
import sys, json, time
import numpy as np
from numba import set_num_threads
sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')
from twist_map import spectra
set_num_threads(4)
base=Path(__file__).parent
out=Path('D:/Codex/chaos-attractor-lab/results/twist-geometry')
out.mkdir(parents=True,exist_ok=True)
source=np.load(base/'pilot-v2/pilot.npz')
ids=np.array([1709,3254,307])
p=source['parameters'][ids]
definition={'parameter_order':['r','a_xy','a_yz','a_zx','k_xy','k_yz','k_zx','h_xy','h_yz','h_zx'],'candidates':[{'id':int(i),'parameters':row.tolist()} for i,row in zip(ids,p)]}
(out/'candidates.json').write_text(json.dumps(definition,indent=2)+'\n')
t=time.time()
windows,tails,ends,means,radii=spectra(p,np.tile([.1,.2,.3],(3,1)),100000,4000000,50000,200000)
np.savez_compressed(out/'long-runs.npz',ids=ids,parameters=p,starts=np.tile([.1,.2,.3],(3,1)),warmup=100000,observed=4000000,window=50000,windows=windows,tail=tails,end=ends,means=means,radii=radii)
print(json.dumps({'stage':'long','seconds':time.time()-t,'spectra':windows.mean(1).tolist(),'window_minima':windows[:,:,0].min(1).tolist(),'sum_error':(windows.sum(2)-3*np.log(p[:,0,None])).max()}),flush=True)
rng=np.random.default_rng(20261105)
starts=rng.uniform(-2,2,(96,3))
pp=np.repeat(p,32,axis=0)
w,z,e,m,r=spectra(pp,starts,100000,500000,25000,8192)
np.savez_compressed(out/'basins.npz',ids=np.repeat(ids,32),parameters=pp,starts=starts,warmup=100000,observed=500000,window=25000,windows=w,tail=z,end=e,means=m,radii=r)
print(json.dumps({'stage':'basins','seconds':time.time()-t,'positive':int((w[:,:,0].min(1)>0).sum()),'min_mean':float(w[:,:,0].mean(1).min())}),flush=True)
