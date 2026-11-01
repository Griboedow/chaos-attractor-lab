"""Five-million-step runs and prescribed fresh-start controls."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from numba import njit,prange
from model import coefficients,integrate,update

@njit(parallel=True,cache=True)
def runs(params,starts,steps,window,stride):
    rates=np.empty((len(params),steps//window,3));cloud=np.empty((len(params),steps//stride,4));minimum=np.empty((len(params),2));end=np.empty((len(params),4));status=np.empty((len(params),2),np.int64)
    for k in prange(len(params)):
        r,c,m,d,w,st,n=integrate(starts[k],params[k],10000,steps,window,stride)
        rates[k]=r;cloud[k]=c;minimum[k]=np.array([m,d]);end[k]=w;status[k]=np.array([st,n])
    return rates,cloud,minimum,end,status

root=Path(__file__).resolve().parent;out=root/'long-01';out.mkdir(exist_ok=False)
eps=np.array([.04,.1,.15,.205]);starts=np.tile([.83,.29,-.13,.11],(4,1));t=time.time()
protocol={'epsilon':eps.tolist(),'start':starts[0].tolist(),'warmup':10000,'observed':5000000,'window':500000,'cloud_stride':10,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'model_sha256':hashlib.sha256((root/'model.py').read_bytes()).hexdigest()}
(out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
r,c,m,end,st=runs(np.array([coefficients(e) for e in eps]),starts,5000000,500000,10)
np.savez_compressed(out/'runs.npz',epsilon=eps,rates=r,clouds=c,minima=m,endpoints=end,status=st)
rows=[{'epsilon':float(e),'mean':r[k].mean(0).tolist(),'rates':r[k].tolist(),'minimum_update_norm':float(m[k,0]),'minimum_qr_diagonal':float(m[k,1]),'endpoint':end[k].tolist(),'status':int(st[k,0]),'completed':int(st[k,1])} for k,e in enumerate(eps)]
(out/'report.json').write_text(json.dumps({'rows':rows,'wall_seconds':time.time()-t},indent=2)+'\n')
print(json.dumps({'stage':'long','rows':rows,'seconds':time.time()-t}),flush=True)
rng=np.random.default_rng(2026110104);starts=rng.normal(size=(128,4));starts/=np.linalg.norm(starts,axis=1)[:,None]
theta=rng.uniform(-np.pi,np.pi,64);u=rng.normal(size=(64,2));u*=rng.uniform(0,.4,64)[:,None]/np.linalg.norm(u,axis=1)[:,None]
starts[:64]=np.c_[np.cos(theta),np.sin(theta),u];starts[:64]/=np.linalg.norm(starts[:64],axis=1)[:,None]
fresh=[]
for e in [.04,.15]:
    t=time.time();coefs=np.tile(coefficients(e),(128,1));r,c,m,end,st=runs(coefs,starts,200000,50000,200000)
    fresh.append({'epsilon':e,'warmup':10000,'observed':200000,'window':50000,'seed':2026110104,'cohorts':['first 64: disk-radius uniform in [0,0.4] inside the certified tube','last 64: normalized independent normal vectors on S3'],
      'starts':starts.tolist(),'rates':r.tolist(),'minima':m.tolist(),'status':st.tolist(),'positive_cases':[int(k) for k in range(128) if st[k,0]==0 and min(r[k,:,0])>.005],'wall_seconds':time.time()-t})
    (out/'fresh-starts.json').write_text(json.dumps(fresh,indent=2)+'\n');print(json.dumps({'stage':'fresh','epsilon':e,'positive':len(fresh[-1]['positive_cases']),'failed':int(np.count_nonzero(st[:,0]))}),flush=True)
# Fixed unthinned inputs for the presentation, distinct from the long-run stride-10 preview.
for e in [.04,.15]:
    r,c,m,d,end,st,n=integrate(np.array([.83,.29,-.13,.11]),coefficients(e),10000,500000,50000,1)
    assert st==0 and n==510000
    name=f'portrait-{round(e*1000):03d}.npz';np.savez_compressed(out/name,states=c)
    (out/name.replace('.npz','.json')).write_text(json.dumps({'epsilon':e,'warmup':10000,'observed':500000,'stride':1,'initial':list([.83,.29,-.13,.11]),'sha256':hashlib.sha256((out/name).read_bytes()).hexdigest(),'rates':r.tolist(),'source_sha256':protocol['source_sha256'],'model_sha256':protocol['model_sha256']},indent=2)+'\n')
print(json.dumps({'stage':'complete','portraits':2}),flush=True)
