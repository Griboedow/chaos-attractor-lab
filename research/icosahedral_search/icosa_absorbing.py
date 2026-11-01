"""Outward interval radial bounds for the actual stored Gaussian coefficients."""
from pathlib import Path
import sys
from icosa_paths import ROOT,RESEARCH,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json
import numpy as np
import mpmath as mp

root=OUTPUT;d=np.load(root/'icosa-v1/input.npz');out=root/'icosa-absorbing';out.mkdir(exist_ok=False)
mp.iv.dps=70;iv=mp.iv
def exact(a):
    n,q=float(a).as_integer_ratio();return iv.mpf(n)/q
C=iv.mpf('1.0000000001');rows=[]
for ident in [177,183,43,29,85]:
    c=[[exact(v) for v in a] for a in d['centers'][ident]];u=[[exact(v) for v in a] for a in d['axes'][ident]]
    norms=[iv.sqrt(sum(v*v for v in a)) for a in c];assert all(v.b<C.a for v in norms)
    U=sum(iv.sqrt(sum(v*v for v in a)) for a in u);s=exact(d['sigma'][ident]);g=exact(d['gamma'][ident])
    for n in range(12,61):
        R=iv.mpf(n)/10;y=R-C
        if y.a<s.b:continue
        vortex=U*y*iv.exp(-y*y/(2*s*s));damping=g*R
        if vortex.b<(iv.mpf('.9')*damping).a:break
    else:raise ValueError('No declared-radius certificate')
    rows.append({'id':ident,'radius':f'{n/10:.1f}','center_norm_bound':'1.0000000001','sum_axis_norms':str(U),'vortex_bound_at_R':str(vortex),'damping_at_R':str(damping),'strict_margin':str(damping-vortex)})
report={'precision':70,'method':'All binary64 coefficients are converted exactly by integer ratios. For r >= R >= C+sigma, each localized term has magnitude <= |u_j| (r-C) exp(-(r-C)^2/(2 sigma^2)). The ratio of this envelope to gamma*r decreases strictly for r >= C+sigma. The strict inequality at R therefore holds for every r >= R. Thus the sphere is inward pointing, the ball is positively invariant, and every bounded set enters it in finite time.','scope':'Global boundedness and an absorbing ball for the actual stored coefficients. This does not prove chaotic attraction or exhaustive basin coverage.','cases':rows}
(out/'certificate.json').write_text(json.dumps(report,indent=2));print(json.dumps(rows),flush=True)
