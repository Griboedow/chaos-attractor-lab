from pathlib import Path
import sys,json,time
import numpy as np
from numba import set_num_threads
sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')
from gaussian_vortex import integrate
set_num_threads(4)
base=Path(__file__).parent
out=Path('D:/Codex/chaos-attractor-lab/results/gaussian-vortex');out.mkdir(exist_ok=True)
source=np.load(base/'vortex-v1/data.npz');ids=np.array([77,233,738]);p=source['parameters'][ids]
(out/'candidates.json').write_text(json.dumps({'parameter_order':['gamma','omega0','omega1','omega2','sigma0','sigma1','sigma2','center_distance','drift0','drift1','drift2'],'candidates':[{'id':int(i),'parameters':v.tolist()} for i,v in zip(ids,p)]},indent=2)+'\n')
pp=np.repeat(p,2,axis=0);starts=np.tile([.1,.2,.3],(6,1));hs=np.tile([.0125,.00625],3)
t=time.time();w,z,e,m,r=integrate(pp,starts,hs,2000.,12000.,500.,50000,4)
np.savez_compressed(out/'long-runs.npz',ids=np.repeat(ids,2),parameters=pp,starts=starts,stepsizes=hs,warm_time=2000.,observed_time=12000.,window_time=500.,tail_stride=4,windows=w,tail=z,end=e,means=m,max_radius=r)
print(json.dumps({'seconds':time.time()-t,'spectra':w.mean(1).tolist(),'min_windows':w[:,:,0].min(1).tolist(),'divergence_errors':(w.sum(2)+3*pp[:,0,None]).max(1).tolist()}),flush=True)
