"""Longer, smaller-step follow-up, selected from a frozen exploratory survey."""
from pathlib import Path
import json,time
import numpy as np
from tetra_vortex import survey

root=Path(__file__).parent
d=np.load(root/'tetra-v1/data.npz')
m=json.loads((root/'tetra-v1/metrics.json').read_text())
ids=list(dict.fromkeys([r['id'] for r in m[:24]]+[r['id'] for r in sorted(m,key=lambda r:r['symmetry_ratio'])[:16]]))
out=root/'tetra-long';out.mkdir(exist_ok=False)
(out/'selection.json').write_text(json.dumps({'ids':ids,'rule':'Union of lowest 24 scout KY proxies and lowest 16 group/self cloud discrepancy ratios. No full-spectrum interpretation of the scout proxy.','warm_time':2000,'observed_time':5000,'h':.025,'keep':50000},indent=2))
begin=time.time()
w,tail,end,speed,radius=survey(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],d['starts'][ids],80000,200000,.025,50000)
np.savez_compressed(out/'data.npz',ids=ids,windows=w,tail=tail,end=end,speed=speed,radius=radius)
print(json.dumps({'seconds':time.time()-begin,'rows':[{'id':int(i),'windows':a.tolist(),'speed':float(b)} for i,a,b in zip(ids,w,speed)]}),flush=True)
