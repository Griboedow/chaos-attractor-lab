"""Step refinement and full spectra for five exploratory icosahedral candidates."""
from pathlib import Path
import json,time
import numpy as np
from tetra_spectra import integrate

root=Path(__file__).parent;d=np.load(root/'icosa-v1/input.npz');out=root/'icosa-spectra';out.mkdir(exist_ok=False)
ids=np.repeat([29,85,183,43,177],3);hs=np.tile([.025,.0125,.00625],5)
(out/'design.json').write_text(json.dumps({'selected_ids':[29,85,183,43,177],'reason':'29 and 85 combine thin scout KY proxies with nearly full-group support; 183 and 43 show visible loops; 177 is a thin but symmetry-broken control. Selection precedes this full-spectrum run.','warm_time':500,'observed_time':2000,'window_time':250,'h':hs.tolist(),'keep':40000},indent=2))
begin=time.time();w,z,means,bounds,ends=integrate(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],d['starts'][ids],hs,500,2000,250,40000)
np.savez_compressed(out/'data.npz',ids=ids,h=hs,windows=w,tail=z,means=means,max_radius=bounds,end=ends)
rows=[{'id':int(i),'h':float(h),'spectrum':a.mean(0).tolist(),'min_leading':float(a[:,0].min()),'sum_error':float(np.max(np.abs(a.sum(1)+3*d['gamma'][i]))),'radius':float(r)} for i,h,a,r in zip(ids,hs,w,bounds)]
(out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2));print(json.dumps(rows),flush=True)
