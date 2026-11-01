"""Unique scalar fixed point, cold-start basins, and attraction controls."""
from pathlib import Path
import json,time
import numpy as np
import mpmath as mp
from scipy.optimize import brentq
from scipy.stats import qmc
from anderson_extended import survey,g_and_derivative

root=Path(__file__).parent;d=np.load(root/'anderson-v2/data.npz');long=np.load(root/'anderson-spectra/data.npz');out=root/'anderson-basins';out.mkdir(exist_ok=True)
mp.iv.dps=60;proofs=[];params=[];starts=[];ids=[];labels=[]
for ident in [54,449]:
    p=d['parameters'][ident];k,c,reg,_=p
    ki=mp.iv.mpf(repr(float(k)));ci=mp.iv.mpf(repr(float(c)));t=mp.iv.sqrt(1-1/ki);a=mp.iv.log((1+t)/(1-t))/2;B=ki*t-a
    margin=ci-B if c>0 else -ci-B
    assert margin.a>0
    fixed=brentq(lambda x:c+k*np.tanh(x)-x,c-k-1,c+k+1,xtol=1e-14)
    proofs.append({'id':ident,'k_decimal':repr(float(k)),'c_decimal':repr(float(c)),'unique_root_margin_interval':str(margin),'root_approx':fixed,'argument':'For r(x)=c+k*tanh(x)-x, stationary points are +/-atanh(sqrt(1-1/k)); their values are c+/-B. The positive interval |c|-B excludes two extra roots. Bounded increasing g implies every plain Picard trajectory converges, hence to this unique root. Interval arithmetic uses 60 decimal digits.'})
    xs=np.linspace(c-k-2,c+k+2,257)
    cold=np.column_stack([xs,c+k*np.tanh(xs),c+k*np.tanh(c+k*np.tanh(xs))])
    arbitrary=(qmc.Sobol(3,scramble=True,seed=2026110140+ident).random_base2(5)*2-1)*(k+abs(c)+2)
    cloud=long['tail'][np.where(long['ids']==ident)[0][0]]
    rng=np.random.default_rng(2026110141+ident);local=cloud[np.linspace(0,len(cloud)-1,16,dtype=int)]+1e-5*rng.normal(size=(16,3))
    near=np.ones((8,3))*fixed+1e-5*rng.normal(size=(8,3))
    for label,points in [('cold_scalar',cold),('arbitrary_history',arbitrary),('local',local),('near_root',near)]:
        for z in points:params.append(p);starts.append(z);ids.append(ident);labels.append(label)
ids=np.array(ids);params=np.array(params);starts=np.array(starts);labels=np.array(labels)
(out/'unique-root.json').write_text(json.dumps(proofs,indent=2))
begin=time.time();w,z,end,status,steps,bounds=survey(params,starts,50000,500000,2000)
np.savez_compressed(out/'data.npz',ids=ids,parameters=params,starts=starts,labels=labels,windows=w,tail=z,end=end,status=status,steps=steps,max_abs=bounds)
rows=[]
for ident in [54,449]:
    for label in ['cold_scalar','arbitrary_history','local','near_root']:
        mask=(ids==ident)&(labels==label);positive=(status==1)&(w.min(1)>.005)
        rows.append({'id':ident,'kind':label,'count':int(mask.sum()),'positive':int((positive&mask).sum()),'fixed':int(((status==0)&mask).sum()),'other':int((mask&~positive&(status!=0)).sum())})
(out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2));print(json.dumps(rows),flush=True)
