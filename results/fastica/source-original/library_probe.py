from pathlib import Path
import inspect,sys,json,hashlib,time
import numpy as np
import sklearn
from sklearn.decomposition import FastICA
from sklearn.decomposition import _fastica
root=Path('D:/Codex/unexpected-attractors-20261101/fastica');sys.path.insert(0,str(root));from follow import map_jac
class ProbeEnd(Exception):pass
rows=[]
for case in [15,40,204,246,351]:
 data=json.loads((root/f'follow-01/case-{case}-input.json').read_text());X=np.array(data['X']);w0=np.array(data['start']);kind=data['kind'];alpha=data['alpha'];fun=[_fastica._logcosh,_fastica._exp,_fastica._cube][kind];kwargs={'alpha':alpha}
 state={'count':0,'sums':np.zeros(5),'min_update':1e9,'max_step_error':0.,'completed_first_unit':False};v=np.ones(4);v-=w0*(v@w0)/(w0@w0);v/=np.linalg.norm(v)
 def observed(u):
  global v
  frame=inspect.currentframe().f_back
  while frame is not None and frame.f_code.co_name!='_ica_def': frame=frame.f_back
  assert frame is not None
  if frame.f_locals['j']!=0:
   state['completed_first_unit']=True;state['terminal_weight']=frame.f_locals['W'][0].copy().tolist();raise ProbeEnd()
  w=frame.f_locals['w'].copy();idx=state['count'];state['count']+=1
  if idx>=52000:raise ProbeEnd()
  f,J,hn=map_jac(w,X,kind,alpha);state['min_update']=min(state['min_update'],hn)
  v=J@v;gain=np.linalg.norm(v);v/=gain
  if idx>=2000:state['sums'][(idx-2000)//10000]+=np.log(gain)/10000
  result=fun(u,kwargs)
  h=(X.T*result[0]).mean(axis=1)-np.mean(result[1])*w;ref=h/np.linalg.norm(h)
  state['max_step_error']=max(state['max_step_error'],float(abs(ref-f).max()))
  return result
 start=time.time();W=np.eye(4);W[0]=w0
 try:FastICA(algorithm='deflation',whiten=False,fun=observed,max_iter=52001,tol=1e-15,w_init=W).fit(X)
 except ProbeEnd:pass
 row={'case':case,'library':'scikit-learn','version':sklearn.__version__,'source_sha256':hashlib.sha256(Path(_fastica.__file__).read_bytes()).hexdigest(),'iterations_observed':min(state['count'],52000),'first_unit_converged':state['completed_first_unit'],'rates':state['sums'].tolist(),'minimum_update_norm':state['min_update'],'maximum_step_identity_error':state['max_step_error'],'wall_seconds':time.time()-start,'hook':'Public FastICA callable g delegates unchanged to installed nonlinearity. Observe first deflation unit; stop after its declared budget or completion.'};rows.append(row);print(json.dumps(row),flush=True)
 counts=[max(0,min(10000,min(state['count'],52000)-2000-10000*k)) for k in range(5)]
 row['window_counts']=counts;row['complete_window_rates']=[float(rate) if count==10000 else None for rate,count in zip(state['sums'],counts)]
 if 'terminal_weight' in state:row['terminal_weight']=state['terminal_weight']
 row['source_probe_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 (root/'library-controls-v2.json').write_text(json.dumps(rows,indent=2)+'\n')
