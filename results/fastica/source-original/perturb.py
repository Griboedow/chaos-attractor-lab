"""Finite input-perturbation ensemble for the strongest fresh-start candidate."""
from pathlib import Path
import json,hashlib,time
import numpy as np
from pilot import batch
root=Path('D:/Codex/unexpected-attractors-20261101/fastica');out=root/'perturb-01';out.mkdir(exist_ok=False);data=json.loads((root/'follow-01/case-204-input.json').read_text());base=np.array(data['X']);start=np.array(data['start']);rng=np.random.default_rng(2026110104);inputs=[];metadata=[]
for epsilon in [1e-6,1e-3,1e-2]:
 for k in range(32):
  seed=int(rng.integers(2**32));a=base+epsilon*np.random.default_rng(seed).normal(size=base.shape);a-=a.mean(0);e,U=np.linalg.eigh(a.T@a/len(a));X=a@((U/np.sqrt(e))@U.T);inputs.append(X);metadata.append({'epsilon':epsilon,'seed':seed,'input_sha256':hashlib.sha256(X.tobytes()).hexdigest()})
protocol={'case':204,'seed':2026110104,'perturbation':'Independent standard Gaussian entries times epsilon, followed by empirical recentering and whitening. Initial weight is unchanged.','warmup':5000,'observed':200000,'window':50000,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');n=len(inputs)
stats,ends=batch(np.array(inputs),np.full(n,len(base),np.int64),np.tile(start,(n,1)),np.full(n,1,np.int64),np.full(n,data['alpha']),5000,200000,50000)
rows=[]
for m,s in zip(metadata,stats):rows.append({**m,'rates':s[:4].tolist(),'minimum_update_norm':s[-4],'period_up_to_sign':int(s[-3]),'status':int(s[-1]),'retained':bool(min(s[:4])>.005 and s[-3]==0 and s[-4]>1e-8 and s[-1]==0)})
(out/'report.json').write_text(json.dumps(rows,indent=2)+'\n')
for e in [1e-6,1e-3,1e-2]:print(json.dumps({'epsilon':e,'retained':sum(r['retained'] for r in rows if r['epsilon']==e)}),flush=True)
