"""Independent high-precision one-step and inverse checks on retained states."""
from pathlib import Path
import hashlib,json
import numpy as np
import mpmath as mp
import model,precision

def inverse(z,L,la,lb,w):
    m,n=L.shape;p=m-1;c=1-w;x=np.r_[0.,z[:p]]
    g=np.array([lb[j]-model.lse(L[:,j]+x) for j in range(n)])
    y=np.r_[0.,(z[p:]-w*(g[1:]-g[0]))/c]
    f=np.array([la[i]-model.lse(L[i]+y) for i in range(m)])
    xx=(z[:p]-w*(f[1:]-f[0]))/c
    return np.r_[xx,y[1:]]

mp.mp.dps=100;rng=np.random.default_rng(201212566);rows=[]
for n,case in [(3,43),(4,12),(4,55),(4,64)]:
    parent=Path(f'extension01/n{n}.json');r=next(r for r in json.loads(parent.read_text())['runs'] if r['case']==case)
    L=np.array(r['log_kernel']);la,lb=model.targets(L);w=r['omega'];K,a,b,ww=precision.prepare(r)
    data=np.load(f'extension01/n{n}-final-25000-steps.npz');states=data['states'][list(data['case']).index(case)]
    for index in (0,1234,9999,20000,24999):
        z=states[index];v=rng.normal(size=len(z));v/=np.linalg.norm(v)
        actual,dv=model.step(z,v,L,la,lb,w)
        high,hdv=precision.step([precision.exact(t) for t in z],[precision.exact(t) for t in v],K,a,b,ww)
        err=float(max(abs(precision.exact(x)-y) for x,y in zip(actual,high)))
        derr=float(max(abs(precision.exact(x)-y) for x,y in zip(dv,hdv)))
        inv=float(np.max(abs(inverse(actual,L,la,lb,w)-z)))
        h=1e-4;fd=(model.step(z+h*v,v,L,la,lb,w)[0]-model.step(z-h*v,v,L,la,lb,w)[0])/(2*h)
        ferr=float(np.max(abs(fd-dv)))
        assert err<2e-10 and derr<2e-11 and inv<2e-10 and ferr<2e-7
        rows.append(dict(shape=[n,n],case=case,final_window_index=index,mp100_update_error=err,
            mp100_direction_error=derr,inverse_error=inv,finite_difference_error=ferr))
Path('actual-check.json').write_text(json.dumps(dict(checks=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
print(json.dumps({key:max(r[key] for r in rows) for key in ('mp100_update_error','mp100_direction_error','inverse_error','finite_difference_error')}))
