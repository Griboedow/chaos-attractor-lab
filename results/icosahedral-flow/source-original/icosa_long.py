"""Long single-orbit checks selected after the 2,000-time-unit refinement."""
from pathlib import Path
import json,time
import numpy as np
from tetra_spectra import integrate

root=Path(__file__).parent
d=np.load(root/'icosa-v1/input.npz')
out=root/'icosa-long';out.mkdir(exist_ok=False)
ids=np.array([183,43,177]);hs=np.array([.0125,.0125,.00625])
design={'selected_ids':ids.tolist(),'reason':'183 and 43 have legible long loops; 177 retains a thinner shell. All passed three-step full-spectrum screening. This run tests long survival and coverage from the original starts.','warm_time':2000,'observed_time':16000,'window_time':2000,'h':hs.tolist(),'keep':1280000}
(out/'design.json').write_text(json.dumps(design,indent=2))
begin=time.time()
w,z,means,bounds,ends=integrate(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],d['starts'][ids],hs,2000,16000,2000,1280000)
rows=[]
for j,ident in enumerate(ids):
    np.save(out/f'{ident}-dense.npy',z[j])
    row={'id':int(ident),'h':float(hs[j]),'spectrum':w[j].mean(0).tolist(),'windows':w[j].tolist(),'min_leading':float(w[j,:,0].min()),'trace_error':float(np.max(np.abs(w[j].sum(1)+3*d['gamma'][ident]))),'radius':float(bounds[j])}
    rows.append(row)
np.savez_compressed(out/'sample.npz',ids=ids,h=hs,windows=w,tail=z[:,::16],means=means,max_radius=bounds,end=ends)
(out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2))
print(json.dumps(rows),flush=True)
