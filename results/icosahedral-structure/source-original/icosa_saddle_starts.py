"""Probe the basin near two representative saddle types of case 177."""
from pathlib import Path
import json,time
import numpy as np
from tetra_spectra import field,integrate

base=Path(__file__).parent;d=np.load(base/'icosa-v1/input.npz');inv=json.loads((base/'icosa-equilibria/inventory.json').read_text());out=base/'icosa-saddle-starts';out.mkdir(exist_ok=False)
case=next(c for c in inv if c['id']==177);ident=177;c=d['centers'][ident];u=d['axes'][ident];s=d['sigma'][ident];g=d['gamma'][ident]
starts=[];labels=[];origins=[]
for radius in [1.0408942,1.2374565]:
    e=np.array(next(r['point'] for r in case['roots'] if abs(r['radius']-radius)<1e-6));j=field(e,c,u,s,g)[1];ev,V=np.linalg.eig(j)
    for eps in [1e-4,1e-6]:
        if radius<1.1:
            k=np.where((ev.real>0)&(ev.imag>0))[0][0];a=V[:,k].real;a/=np.linalg.norm(a);b=V[:,k].imag;b-=a*(a@b);b/=np.linalg.norm(b)
            directions=[a,-a,b,-b]
        else:
            k=np.argmax(ev.real);a=V[:,k].real;a/=np.linalg.norm(a);directions=[a,-a]
        for k,v in enumerate(directions):starts.append(e+eps*v);origins.append(e);labels.append({'equilibrium_radius':radius,'epsilon':eps,'direction_index':k})
starts=np.array(starts);ids=np.full(len(starts),ident);hs=np.full(len(starts),.00625)
(out/'design.json').write_text(json.dumps({'labels':labels,'warm_time':500,'observed_time':2000,'window_time':500,'h':.00625,'limits':'Finite perturbations near numerical saddle points; this cannot prove a homoclinic connection or an arbitrary-neighborhood basin property.'},indent=2))
np.savez_compressed(out/'input.npz',starts=starts,equilibria=np.array(origins));begin=time.time()
w,z,means,radius,end=integrate(d['centers'][ids],d['axes'][ids],d['sigma'][ids],d['gamma'][ids],starts,hs,500,2000,500,16000)
np.savez_compressed(out/'data.npz',windows=w,tail=z,means=means,radius=radius,end=end)
rows=[dict(label,mean_spectrum=a.mean(0).tolist(),min_leading=float(a[:,0].min()),end_radius=float(np.linalg.norm(e))) for label,a,e in zip(labels,w,end)]
(out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2));print(json.dumps(rows),flush=True)
