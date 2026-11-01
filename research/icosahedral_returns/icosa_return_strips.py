"""Sample two narrow candidate return strips; no covering certificate."""
from pathlib import Path
from return_paths import INPUT,STRUCTURE,CYCLES,OUTPUT
import json,time
import numpy as np
from icosa_poincare import evaluate

base=OUTPUT;out=base/'icosa-return-strips';out.mkdir(exist_ok=False)
strips=[[-.000008,.000002],[.002105,.00212]];vs=[-.0002,0.,.0002];rows=[];begin=time.time()
for strip,(lo,hi) in enumerate(strips):
    for a in np.linspace(lo,hi,9):
        for b in vs:
            try:
                q,t,J,speed=evaluate(np.array([a,b]));row={'strip':strip,'input':[float(a),b],'returned':True,'output':q.tolist(),'time':t,'derivative':J.tolist(),'speed':speed}
            except RuntimeError as error:row={'strip':strip,'input':[float(a),b],'returned':False,'error':str(error)}
            rows.append(row)
        print(json.dumps({'strip':strip,'u':float(a),'outputs':[r.get('output') for r in rows[-3:]],'seconds':time.time()-begin}),flush=True)
(out/'report.json').write_text(json.dumps({'strips':strips,'stable_levels':vs,'rows':rows,'seconds':time.time()-begin,'limits':'Finite samples of two-dimensional return maps. No interval bounds, all-point coverage, invariant set or entropy theorem.'},indent=2))
