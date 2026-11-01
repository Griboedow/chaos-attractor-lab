"""Observe the actual public FastICA first-unit loop on the rounded data."""
from pathlib import Path
import hashlib,inspect,json,time
import numpy as np
import sklearn
from sklearn.decomposition import FastICA,_fastica
from model import dataset,coefficients,update

class ProbeEnd(Exception):pass
rows=[];root=Path(__file__).resolve().parent
for e in [.04,.15]:
    X=dataset(e);start=np.array([.83,.29,-.13,.11]);W=np.eye(4);W[0]=start;v=np.ones(4);w=start/np.linalg.norm(start);v-=w*(v@w);v/=np.linalg.norm(v)
    record={'count':0,'rates':np.zeros(5),'min_update':1e100,'max_polynomial_step_error':0.,'first_unit_finished':False,'maximum_tube_ratio':0.};coef=coefficients(e)
    def observed(u):
        global v
        frame=inspect.currentframe().f_back
        while frame is not None and frame.f_code.co_name!='_ica_def':frame=frame.f_back
        assert frame is not None
        if frame.f_locals['j']!=0:
            record['first_unit_finished']=True;record['terminal_weight']=frame.f_locals['W'][0].tolist();raise ProbeEnd()
        idx=record['count'];w=frame.f_locals['w'].copy()
        if idx>=52000:raise ProbeEnd()
        record['count']+=1;g,gp=_fastica._cube(u,{})
        # Derivative of this rounded empirical map, including its covariance residual.
        a=np.mean(gp);h=np.mean(X*g[:,None],axis=0)-a*w;H=(X.T*(3*u*u))@X/len(X)-a*np.eye(4)-np.outer(w,np.mean(X*(6*u)[:,None],axis=0))
        hn=np.linalg.norm(h);f=h/hn;z=H@v;v=(z-f*(f@z))/hn;gain=np.linalg.norm(v);v/=gain
        record['min_update']=min(record['min_update'],hn);record['maximum_tube_ratio']=max(record['maximum_tube_ratio'],float(np.linalg.norm(w[2:])/np.linalg.norm(w[:2])))
        ref,_,_=update(w,coef);record['max_polynomial_step_error']=max(record['max_polynomial_step_error'],float(abs(f-ref).max()))
        if idx>=2000:record['rates'][(idx-2000)//10000]+=np.log(gain)/10000
        return g,gp
    t=time.time()
    try:FastICA(algorithm='deflation',whiten=False,fun=observed,max_iter=52001,tol=1e-15,w_init=W).fit(X)
    except ProbeEnd:pass
    counts=[max(0,min(10000,record['count']-2000-j*10000)) for j in range(5)]
    record.update(epsilon=e,window_counts=counts,rates=[float(x) if n==10000 else None for x,n in zip(record['rates'],counts)],library='scikit-learn',version=sklearn.__version__,library_source_sha256=hashlib.sha256(Path(_fastica.__file__).read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model_sha256=hashlib.sha256((root/'model.py').read_bytes()).hexdigest(),data_sha256=hashlib.sha256(X.astype('<f8').tobytes()).hexdigest(),covariance_error=float(abs(X.T@X/len(X)-np.eye(4)).max()),wall_seconds=time.time()-t,hook='Public callable delegates unchanged to installed _cube. Observe first deflation unit without replacing its update; budget interrupt or actual library stop. Rounded literal data; whiten=False.')
    rows.append(record);(root/'library.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(record),flush=True)
