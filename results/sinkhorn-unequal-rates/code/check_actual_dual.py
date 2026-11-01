"""Check retained states against independent probability arithmetic and inversion."""
from pathlib import Path
import hashlib,json
import mpmath as mp
import numpy as np
from scipy.special import logsumexp
import model,model_dual,precision_dual

mp.mp.dps=100;rows=[];rng=np.random.default_rng(201212903)
inputs=[('balanced-long01/3x2.json',[150,987,732,717,703]),('dual-long01/2x2.json',[156,193,239,709]),('dual-long01/2x3.json',[513,627]),('dual-long01/3x2.json',[104,127,209,428,637]),('rectangular-long01/2x4.json',[326,279])]
for filename,cases in inputs:
    data=json.loads(Path(filename).read_text())
    for case in cases:
        r=next(r for r in data['runs'] if r['case']==case);L=np.array(r['log_kernel']);m,n=L.shape;p=m-1;la,lb=model.targets(L);wx=r.get('wx',r.get('omega'));wy=r.get('wy',r.get('omega'));z=np.array(r['blocks'][-1]['checkpoint']);v=rng.normal(size=len(z));v/=np.linalg.norm(v)
        zz,dz=model_dual.step(z,v,L,la,lb,wx,wy);K,a,b,mwx,mwy=precision_dual.prepare(r)
        mz,md=precision_dual.step([precision_dual.exact(x) for x in z],[precision_dual.exact(x) for x in v],K,a,b,mwx,mwy)
        state_error=float(max(abs(precision_dual.exact(x)-y) for x,y in zip(zz,mz)));derivative_error=float(max(abs(precision_dual.exact(x)-y) for x,y in zip(dz,md)))
        h=1e-5;fd=(model_dual.step(z+h*v,v,L,la,lb,wx,wy)[0]-model_dual.step(z-h*v,v,L,la,lb,wx,wy)[0])/(2*h);fd_error=float(np.max(abs(fd-dz)))
        xn=np.r_[0.,zz[:p]];g=lb-logsumexp(L+xn[:,None],axis=0);old_y=(zz[p:]-wy*(g[1:]-g[0]))/(1-wy);f=la-logsumexp(L+np.r_[0.,old_y][None,:],axis=1);old_x=(zz[:p]-wx*(f[1:]-f[0]))/(1-wx);inverse_error=float(np.max(abs(z-np.r_[old_x,old_y])))
        J=np.column_stack([model_dual.step(z,e,L,la,lb,wx,wy)[1] for e in np.eye(len(z))]);det_error=float(abs(np.linalg.det(J)-(1-wx)**p*(1-wy)**(n-1)))
        assert state_error<5e-12 and derivative_error<5e-12 and fd_error<1e-7 and inverse_error<5e-11 and det_error<5e-12
        rows.append(dict(parent=filename,parent_sha256=hashlib.sha256(Path(filename).read_bytes()).hexdigest(),case=case,shape=[m,n],state=z.tolist(),direction=v.tolist(),state_error=state_error,derivative_error=derivative_error,finite_difference_error=fd_error,inverse_error=inverse_error,determinant_error=det_error))
result=dict(checks=rows,precision_decimal_digits=100,sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path('model.py'),Path('model_dual.py'),Path('precision_dual.py'))})
Path('actual-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(checks=len(rows),max_errors={k:max(r[k] for r in rows) for k in ('state_error','derivative_error','finite_difference_error','inverse_error','determinant_error')})))
