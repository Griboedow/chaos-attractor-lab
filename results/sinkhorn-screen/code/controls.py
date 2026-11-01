"""Declared local, ordinary-Sinkhorn and nearby-start controls."""
from pathlib import Path
import json,hashlib
import numpy as np
from numba import njit
import model

@njit(cache=True)
def run(z,L,w,warm=10000,observed=200000):
    la,lb=model.targets(L);v=np.ones(len(z))/np.sqrt(len(z));rates=np.zeros(10);maximum=0.
    for k in range(warm+observed):
        z,v=model.step(z,v,L,la,lb,w);norm=np.linalg.norm(v);v/=max(norm,1e-300)
        maximum=max(maximum,np.max(np.abs(z)))
        if k>=warm:rates[(k-warm)//20000]+=np.log(max(norm,1e-300))/20000
    residual=np.max(np.abs(model.step(z,v,L,la,lb,w)[0]-z))
    return z,rates,maximum,residual

rows=[];rng=np.random.default_rng(201212565)
for n,case in [(3,43),(4,12),(4,55),(4,64)]:
    r=next(r for r in json.loads(Path(f'extension01/n{n}.json').read_text())['runs'] if r['case']==case)
    original=json.loads(Path(f'pilot02/n{n}.json').read_text())['runs'][case]
    L=np.array(r['log_kernel']);end=np.array(r['blocks'][-1]['checkpoint'])
    setups=[('ordinary',1.,np.array(original['initial_state'])),('near_origin',r['omega'],rng.normal(size=len(end))*1e-6)]
    setups.extend((f'near_endpoint_{i}',r['omega'],end+rng.normal(size=len(end))*.01) for i in range(4))
    for name,w,z in setups:
        final,rates,maximum,residual=run(z,L,w)
        rows.append(dict(shape=[n,n],case=case,control=name,omega=w,initial_state=z.tolist(),warm=10000,observed=200000,
            block_length=20000,rates=None if name=='ordinary' else rates.tolist(),
            growth_status='Not assessed for the singular ordinary update; a tangent can become exactly zero.' if name=='ordinary' else 'Finite-window estimate',
            final_state=final.tolist(),max_abs_state=float(maximum),final_step_residual=float(residual)))
Path('controls02.json').write_text(json.dumps(dict(runs=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
print(json.dumps([dict(case=r['case'],shape=r['shape'],control=r['control'],last=r['rates'][-1] if r['rates'] else None,norm=max(abs(v) for v in r['final_state'])) for r in rows]))
