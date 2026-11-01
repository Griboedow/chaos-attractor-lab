"""Consecutive states; no sign folding or temporal decimation."""
from pathlib import Path
import json,hashlib
import numpy as np
from follow import integrate
root=Path('D:/Codex/unexpected-attractors-20261101/fastica');out=root/'portraits-01';out.mkdir(exist_ok=False)
protocol={'cases':[15,204,246,351],'warmup':10000,'observed':1000000,'stride':1,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
for case in protocol['cases']:
 d=json.loads((root/f'follow-01/case-{case}-input.json').read_text());r,x,m,q,w=integrate(np.array(d['X']),np.array(d['start']),d['kind'],d['alpha'],10000,1000000,100000,1)
 np.savez_compressed(out/f'case-{case}.npz',cloud=x,rates=r,endpoint=w)
 print(json.dumps({'case':case,'min':x.min(0).tolist(),'max':x.max(0).tolist(),'successive_dot_range':[float(np.min(np.sum(x[1:]*x[:-1],axis=1))),float(np.max(np.sum(x[1:]*x[:-1],axis=1)))]}),flush=True)
