"""Independent probability-domain state/tangent iteration at 70 and 100 digits."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,hashlib,json,time
import mpmath as mp
import numpy as np

def exact(value):
    p,q=float(value).as_integer_ratio();return mp.mpf(p)/q

def prepare(row):
    K=[[mp.exp(exact(v)) for v in line] for line in row['log_kernel']]
    a=[mp.fsum(line) for line in K];b=[mp.fsum(line[j] for line in K) for j in range(len(K[0]))]
    return K,a,b,exact(row.get('wx',row.get('omega'))),exact(row.get('wy',row.get('omega')))

def step(z,v,K,a,b,wx,wy):
    m,n=len(a),len(b);x=[mp.mpf(0)]+z[:m-1];y=[mp.mpf(0)]+z[m-1:];vx=[mp.mpf(0)]+v[:m-1];vy=[mp.mpf(0)]+v[m-1:]
    ey=[mp.exp(t) for t in y];s=[mp.fsum(k*t for k,t in zip(row,ey)) for row in K]
    f=[mp.log(aa/ss) for aa,ss in zip(a,s)];df=[-mp.fsum(k*t*vv for k,t,vv in zip(row,ey,vy))/ss for row,ss in zip(K,s)]
    xn=[mp.mpf(0)]+[(1-wx)*x[i]+wx*(f[i]-f[0]) for i in range(1,m)];vn=[mp.mpf(0)]+[(1-wx)*vx[i]+wx*(df[i]-df[0]) for i in range(1,m)]
    ex=[mp.exp(t) for t in xn];s=[mp.fsum(K[i][j]*ex[i] for i in range(m)) for j in range(n)]
    g=[mp.log(bb/ss) for bb,ss in zip(b,s)];dg=[-mp.fsum(K[i][j]*ex[i]*vn[i] for i in range(m))/s[j] for j in range(n)]
    yn=[(1-wy)*y[j]+wy*(g[j]-g[0]) for j in range(1,n)];vd=[(1-wy)*vy[j]+wy*(dg[j]-dg[0]) for j in range(1,n)]
    return xn[1:]+yn,vn[1:]+vd

def run(job):
    row,dps,warm,observed,folder=job;mp.mp.dps=dps;t=time.perf_counter();K,a,b,wx,wy=prepare(row)
    start=row['blocks'][-1]['checkpoint'];z=[exact(v) for v in start];d=len(z);v=[1/mp.sqrt(d)]*d;rates=[];states=[];total=mp.mpf(0);complete=True
    for k in range(warm+observed):
        z,v=step(z,v,K,a,b,wx,wy)
        if not all(mp.isfinite(zz) for zz in z) or max(abs(zz) for zz in z)>100000:complete=False;break
        norm=mp.sqrt(mp.fsum(vv*vv for vv in v));v=[vv/norm for vv in v]
        if k>=warm:
            total+=mp.log(norm);states.append([float(zz) for zz in z])
            if (k-warm+1)%2000==0:rates.append(float(total/2000));total=mp.mpf(0)
    states=np.asarray(states).reshape((-1,d));key=f"{len(a)}x{len(b)}-case{row['case']}-{dps}digits"
    np.save(Path(folder)/(key+'-all-observed-states.npy'),states)
    result=dict(case=row['case'],shape=[len(a),len(b)],precision_decimal_digits=dps,checkpoint_step=row['blocks'][-1]['stop_step'],initial_state=start,warm=warm,requested_observed=observed,completed=complete,accepted_steps=k+1 if complete else k,observed=len(states),input_note='Exact integer-ratio conversion of binary64 kernel, rates and checkpoint; targets recomputed at working precision.',block_length=2000,rates=rates,mean=states.mean(0).tolist() if len(states) else None,std=states.std(0).tolist() if len(states) else None,final_state=[mp.nstr(zz,dps) for zz in z],elapsed_seconds=time.perf_counter()-t)
    (Path(folder)/(key+'.json')).write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--cases',type=int,nargs='+',required=True);ap.add_argument('--warm',type=int,default=2000);ap.add_argument('--observed',type=int,default=20000)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);rows=[r for r in json.loads(args.input.read_text())['runs'] if r['case'] in args.cases]
    assert len(rows)==len(args.cases)
    with ProcessPoolExecutor(max_workers=2) as pool:
        jobs=[pool.submit(run,(r,dps,args.warm,args.observed,str(args.output))) for r in rows for dps in (70,100)]
        for f in as_completed(jobs):
            r=f.result();print(json.dumps(dict(case=r['case'],digits=r['precision_decimal_digits'],completed=r['completed'],rates=r['rates'],elapsed=r['elapsed_seconds'])),flush=True)
    (args.output/'metadata.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parent_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),mpmath=mp.__version__),indent=2)+'\n')
