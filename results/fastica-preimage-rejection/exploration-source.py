from pathlib import Path
import importlib.util, json, time, hashlib
import numpy as np
from scipy.optimize import least_squares
ROOT=Path('D:/Codex/chaos-attractor-lab')
MODEL=ROOT/'results/fastica-tripling/code/model.py'
spec=importlib.util.spec_from_file_location('structured',MODEL); model=importlib.util.module_from_spec(spec); spec.loader.exec_module(model)
e=.04; coef=model.coefficients(e); K=.4
out=Path(__file__).resolve().parent/'pilot-01';out.mkdir(exist_ok=False)
def chart(p):
 t,u,v=p; raw=np.array([np.cos(t),np.sin(t),u,v]); n=np.sqrt(1+u*u+v*v); w=raw/n
 basis=np.column_stack([np.array([-np.sin(t),np.cos(t),0,0])/n,np.array([0,0,1,0])/n-w*u/(n*n),np.array([0,0,0,1])/n-w*v/(n*n)])
 return w,basis
w=np.array([.83,.29,-.13,.11]);w/=np.linalg.norm(w)
for _ in range(10000):w,_,_=model.update(w,coef)
targets=[]
for i in range(16):
 for _ in range(137):w,_,_=model.update(w,coef)
 y,_,_=model.update(w,coef);targets.append((w.copy(),y.copy()))
rows=[];beg=time.time()
for i,(previous,y) in enumerate(targets):
 def fun(p):return model.update(chart(p)[0],coef)[0]-y
 def jac(p):
  w,B=chart(p);return model.update(w,coef)[1]@B
 theta=np.arctan2(y[1],y[0]); seeds=[]
 for k in range(3):
  for u in [-.28,0,.28]:
   for v in [-.28,0,.28]:seeds.append([(-theta+2*np.pi*k)/3,u,v])
 seeds.append([np.arctan2(previous[1],previous[0]),previous[2]/np.linalg.norm(previous[:2]),previous[3]/np.linalg.norm(previous[:2])])
 attempts=[]; roots=[]
 for seed in seeds:
  sol=least_squares(fun,seed,jac=jac,bounds=([-np.inf,-K,-K],[np.inf,K,K]),ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=500)
  q,_=chart(sol.x); res=float(np.linalg.norm(fun(sol.x)));ratio=float(np.linalg.norm(sol.x[1:]));idx=None
  if res<1e-10 and ratio<=K:
   d=[np.linalg.norm(q-r['w']) for r in roots]
   if not d or min(d)>1e-6:
    roots.append({'w':q.tolist(),'residual':res,'ratio':ratio,'known_distance':float(np.linalg.norm(q-previous))});idx=len(roots)-1
   else:idx=int(np.argmin(d))
  attempts.append({'seed':seed,'status':int(sol.status),'nfev':sol.nfev,'residual':res,'ratio':ratio,'root':idx})
 row={'target_index':i,'known_previous':previous.tolist(),'target':y.tolist(),'roots':roots,'attempts':attempts};rows.append(row)
 print(json.dumps({'target':i,'roots':len(roots),'distances':[r['known_distance'] for r in roots],'seconds':time.time()-beg}),flush=True)
(out/'report.json').write_text(json.dumps({'epsilon':e,'tube_ratio':K,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'model_sha256':hashlib.sha256(MODEL.read_bytes()).hexdigest(),'rows':rows},indent=2)+'\n')
