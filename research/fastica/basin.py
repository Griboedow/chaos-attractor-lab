"""Declared finite ensemble of fresh FastICA starts, including negative cases."""
from pathlib import Path
from paths import output_root
import json,hashlib,time
import numpy as np
from pilot import batch
root=output_root()
out=root/'basin-01';out.mkdir(exist_ok=False)
protocol={'cases':[15,40,204,246,351],'starts_per_case':64,'seed':2026110102,'warmup':5000,'observed':200000,'window':50000,'gate':'All four leading rates > 0.005; no recurrence at lags 1..128 up to sign; update norm >1e-8; status=0.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
rng=np.random.default_rng(protocol['seed']);X=np.zeros((320,64,4));ns=np.empty(320,np.int64);ks=np.empty_like(ns);alphas=np.empty(320);starts=rng.normal(size=(320,4));starts/=np.linalg.norm(starts,axis=1)[:,None]
for j,case in enumerate(protocol['cases']):
    d=json.loads((root/f'follow-01/case-{case}-input.json').read_text());sl=slice(64*j,64*(j+1));x=np.array(d['X']);X[sl,:len(x)]=x;ns[sl]=len(x);ks[sl]=d['kind'];alphas[sl]=d['alpha']
(out/'starts.json').write_text(json.dumps(starts.tolist())+'\n')
t=time.time();stats,ends=batch(X,ns,starts,ks,alphas,5000,200000,50000);rows=[]
for i,s in enumerate(stats):rows.append({'case':protocol['cases'][i//64],'start':i%64,'rates':s[:4].tolist(),'minimum_update_norm':s[-4],'period_up_to_sign':int(s[-3]),'recurrence_error':s[-2],'status':int(s[-1]),'retained':bool(min(s[:4])>.005 and s[-3]==0 and s[-4]>1e-8 and s[-1]==0)})
(out/'report.json').write_text(json.dumps({'wall_seconds':time.time()-t,'rows':rows},indent=2)+'\n')
for case in protocol['cases']:
    cases=[r for r in rows if r['case']==case];print(json.dumps({'case':case,'retained':sum(r['retained'] for r in cases),'periodic':sum(r['period_up_to_sign']>0 for r in cases),'rate_range':[min(min(r['rates']) for r in cases),max(max(r['rates']) for r in cases)]}),flush=True)
