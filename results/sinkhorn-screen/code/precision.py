"""Independent arbitrary-precision probability-domain Sinkhorn checks."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,hashlib,json,time
import mpmath as mp
import numpy as np

def exact(value):
    p,q=float(value).as_integer_ratio()
    return mp.mpf(p)/q

def prepare(row):
    L=[[exact(v) for v in line] for line in row['log_kernel']]
    K=[[mp.exp(v) for v in line] for line in L]
    a=[mp.fsum(line) for line in K]
    b=[mp.fsum(line[j] for line in K) for j in range(len(K[0]))]
    return K,a,b,exact(row['omega'])

def step(z,v,K,a,b,w):
    m,n=len(a),len(b);c=1-w
    x=[mp.mpf(0)]+z[:m-1];y=[mp.mpf(0)]+z[m-1:]
    vx=[mp.mpf(0)]+v[:m-1];vy=[mp.mpf(0)]+v[m-1:]
    ey=[mp.exp(t) for t in y];s=[mp.fsum(k*t for k,t in zip(row,ey)) for row in K]
    f=[mp.log(aa/ss) for aa,ss in zip(a,s)]
    df=[-mp.fsum(k*t*vv for k,t,vv in zip(row,ey,vy))/ss for row,ss in zip(K,s)]
    xn=[mp.mpf(0)]+[c*x[i]+w*(f[i]-f[0]) for i in range(1,m)]
    vn=[mp.mpf(0)]+[c*vx[i]+w*(df[i]-df[0]) for i in range(1,m)]
    ex=[mp.exp(t) for t in xn];s=[mp.fsum(K[i][j]*ex[i] for i in range(m)) for j in range(n)]
    g=[mp.log(bb/ss) for bb,ss in zip(b,s)]
    dg=[-mp.fsum(K[i][j]*ex[i]*vn[i] for i in range(m))/s[j] for j in range(n)]
    yn=[c*y[j]+w*(g[j]-g[0]) for j in range(1,n)]
    vd=[c*vy[j]+w*(dg[j]-dg[0]) for j in range(1,n)]
    return xn[1:]+yn,vn[1:]+vd

def run(job):
    row,dps,warm,observed,folder=job;mp.mp.dps=dps;t=time.perf_counter()
    K,a,b,w=prepare(row);z=[exact(v) for v in row['blocks'][-1]['checkpoint']]
    d=len(z);v=[1/mp.sqrt(d)]*d;rates=[];states=np.empty((observed,d));total=mp.mpf(0)
    for k in range(warm+observed):
        z,v=step(z,v,K,a,b,w);norm=mp.sqrt(mp.fsum(vv*vv for vv in v));v=[vv/norm for vv in v]
        if k>=warm:
            total+=mp.log(norm);states[k-warm]=[float(zz) for zz in z]
            if (k-warm+1)%2000==0:rates.append(float(total/2000));total=mp.mpf(0)
    key=f"n{len(K)}-case{row['case']}-{dps}digits"
    np.save(Path(folder)/(key+'-all-observed-states.npy'),states)
    result=dict(case=row['case'],shape=[len(K),len(K[0])],precision_decimal_digits=dps,
        warm=warm,observed=observed,checkpoint_step=row['blocks'][-1]['stop_step'],
        initial_state=row['blocks'][-1]['checkpoint'],input_note='Binary64 kernel entries, omega and checkpoint converted exactly via integer ratios; targets recomputed at working precision.',
        block_length=2000,rates=rates,mean=states.mean(0).tolist(),std=states.std(0).tolist(),
        min_state=states.min(0).tolist(),max_state=states.max(0).tolist(),
        final_state=[mp.nstr(zz,dps) for zz in z],elapsed_seconds=time.perf_counter()-t)
    (Path(folder)/(key+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--case',type=int,default=55)
    ap.add_argument('--warm',type=int,default=2000);ap.add_argument('--observed',type=int,default=20000)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    row=next(r for r in json.loads(args.input.read_text())['runs'] if r['case']==args.case)
    with ProcessPoolExecutor(max_workers=2) as pool:
        jobs=[pool.submit(run,(row,dps,args.warm,args.observed,str(args.output))) for dps in (70,100)]
        for f in as_completed(jobs):
            r=f.result();print(json.dumps(dict(case=r['case'],digits=r['precision_decimal_digits'],rates=r['rates'],elapsed=r['elapsed_seconds'])),flush=True)
    meta=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parent_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),mpmath=mp.__version__)
    (args.output/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
