"""Arbitrary-precision squared-distance posterior and nested-reflection control."""
import argparse,hashlib,itertools,json,sys,time
from pathlib import Path
import numpy as np


def iterate(x,beta,tau,levels,mp):
    mean=[]
    for j,grid in enumerate(levels):
        logs=[-(x[j]-b)**2/tau for b in grid];top=max(logs)
        w=[mp.exp(v-top) for v in logs];total=sum(w)
        mean.append(sum(a*b for a,b in zip(w,grid))/total)
    rb=[2*b-a for a,b in zip(x,mean)];r=mp.sqrt(sum(v*v for v in rb))
    assert r>0
    ra=[2*v/r-v for v in rb]
    return [beta*(a+b)/2+(1-beta)*c for a,b,c in zip(ra,x,mean)],r


def main(args):
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    args.out.mkdir(parents=True,exist_ok=False);mp.mp.dps=args.digits
    def save(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    protocol=dict(ids=args.ids,digits=args.digits,source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),warm=10000,observed=50000,windows=10,initial=[.317,-.211,.419],method='Exact binary64 input lifts; scalar squared-distance Gaussian weights; nested reflections; finite separation without the production Jacobian.',mpmath=mp.__version__)
    save(args.out/'protocol.json',protocol);(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    with np.load(args.source) as f:raw={i:(f['params'][i].copy(),f['levels'][i,:].copy(),int(f['orders'][i])) for i in args.ids}
    rows=[]
    for ident in args.ids:
        begin=time.time();p,grid,order=raw[ident];beta,tau=[mp.mpf(float(a)) for a in p];levels=[[mp.mpf(float(v)) for v in row[:order]] for row in grid]
        eps=mp.mpf(10)**(-args.digits//2+5);x=[mp.mpf(v) for v in [.317,-.211,.419]];v=[mp.mpf(v) for v in [1.,.371,-.217]];length=mp.sqrt(sum(a*a for a in v));y=[a+eps*b/length for a,b in zip(x,v)]
        rates=[mp.mpf(0) for _ in range(10)];minimum=mp.inf;sample=[]
        for it in range(60000):
            x,r=iterate(x,beta,tau,levels,mp);y,_=iterate(y,beta,tau,levels,mp);minimum=min(minimum,r)
            d=[b-a for a,b in zip(x,y)];length=mp.sqrt(sum(v*v for v in d));assert length>0
            if it>=10000:rates[(it-10000)//5000]+=mp.log(length/eps)
            y=[a+eps*b/length for a,b in zip(x,d)]
            if it>=50000:sample.append([float(v) for v in x])
        estimates=[v/5000 for v in rates];ev=np.linalg.eigvalsh(np.cov(np.array(sample).T))
        row=dict(id=ident,digits=args.digits,epsilon=str(eps),mean=float(sum(estimates)/10),minimum_window=float(min(estimates)),windows=[str(v) for v in estimates],minimum_radius=str(minimum),variance_fractions=(ev/ev.sum()).tolist(),seconds=time.time()-begin)
        rows.append(row);np.savez_compressed(args.out/f'case-{ident}.npz',states=np.array(sample),rates=np.array([float(v) for v in estimates]));save(args.out/f'case-{ident}.json',row);print(json.dumps(row),flush=True)
    save(args.out/'report.json',dict(protocol,cases=rows,scope='Finite precision sensitivity control; not a certified attractor or priority claim.'))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--mpmath-path',type=Path);p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--digits',type=int,required=True);main(p.parse_args())
