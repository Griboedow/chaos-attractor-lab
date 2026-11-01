"""Arbitrary precision check of the actual nonlinear map and tangent vector."""
from pathlib import Path
import json,time
import numpy as np
import mpmath as mp

root=Path(__file__).parent;d=np.load(root/'anderson-v2/data.npz');out=root/'anderson-precision';out.mkdir(exist_ok=False)
rows=[]
for ident,dps in [(54,80),(54,110),(449,80),(2268,80)]:
    begin=time.time();mp.mp.dps=dps
    k,c,reg,kind=[mp.mpf(repr(float(a))) for a in d['parameters'][ident]]
    z=[mp.mpf(repr(float(a))) for a in d['starts'][ident]];v=[mp.mpf(1),mp.mpf(0),mp.mpf(0)];acc=mp.mpf(0);windows=[];tail=[]
    for t in range(60000):
        if int(kind)==4:g=[c+k*x for x in z];gp=[k]*3
        else:
            tan=[mp.tanh(x) for x in z];g=[c+k*a for a in tan];gp=[k*(1-a*a) for a in tan]
        r=[a-b for a,b in zip(g,z)];mu=sum(r)/3;dev=[a-mu for a in r];den=reg+sum(a*a for a in dev);gm=sum(g)/3;cross=sum(a*b for a,b in zip(dev,g))
        grad=[]
        for j in range(3):
            rp=gp[j]-1;grad.append(gp[j]/3-(rp*cross/3+mu*((g[j]-gm)*rp+dev[j]*gp[j]))/den+mu*cross*2*dev[j]*rp/(den*den))
        z=[z[1],z[2],gm-mu*cross/den];v=[v[1],v[2],sum(a*b for a,b in zip(grad,v))];norm=mp.sqrt(sum(a*a for a in v));v=[a/norm for a in v]
        if t>=10000:
            acc+=mp.log(norm)
            if (t-10000+1)%10000==0:windows.append(float(acc/10000));acc=mp.mpf(0)
            if t>=59000:tail.append([float(a) for a in z])
    np.savez_compressed(out/f'{ident}-{dps}.npz',tail=np.array(tail))
    row={'id':ident,'dps':dps,'warm':10000,'observed':50000,'windows':windows,'mean':float(np.mean(windows)),'seconds':time.time()-begin};rows.append(row)
    (out/'summary.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
