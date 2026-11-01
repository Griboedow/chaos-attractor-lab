"""Observe the first library deflation unit under its usual tolerance."""
from pathlib import Path
import inspect,json,hashlib
import numpy as np
from sklearn.decomposition import FastICA,_fastica
from pilot import run
root=Path('D:/Codex/unexpected-attractors-20261101/fastica');rows=[]
class Stop(Exception):pass
for case in [15,40,204,246,351]:
 d=json.loads((root/f'follow-01/case-{case}-input.json').read_text());X=np.array(d['X']);w0=np.array(d['start']);fun=[_fastica._logcosh,_fastica._exp,_fastica._cube][d['kind']];state={'iterations':0}
 def observe(u):
  frame=inspect.currentframe().f_back
  while frame is not None and frame.f_code.co_name!='_ica_def':frame=frame.f_back
  assert frame is not None
  if frame.f_locals['j']!=0:state['terminal']=frame.f_locals['W'][0].copy();raise Stop()
  state['iterations']+=1
  return fun(u,{'alpha':d['alpha']})
 W=np.eye(4);W[0]=w0
 try:FastICA(algorithm='deflation',whiten=False,fun=observe,max_iter=200000,tol=1e-4,w_init=W).fit(X)
 except Stop:pass
 w=state['terminal'];r=run(X,w,d['kind'],d['alpha'],1000,200000,50000)
 rows.append({'case':case,'tolerance':1e-4,'max_iter':200000,'iterations_before_first_unit_exit':state['iterations'],'terminal':w.tolist(),'continuation_warmup':1000,'continuation_observed':200000,'continuation_rates':r[0].tolist(),'period_up_to_sign':int(r[2]),'status':int(r[4])})
print(json.dumps(rows));(root/'stopping.json').write_text(json.dumps({'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Default tolerance, raised max_iter, frozen whitened inputs, explicit first-unit initialization. No claim about the rest of the demixing matrix.','rows':rows},indent=2)+'\n')
