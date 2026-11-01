"""Independent posterior logits, nested reflections and finite separations."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np


def iterate(x,beta,tau,points,mp):
    logits=[(2*sum(a*b for a,b in zip(x,p))-sum(b*b for b in p))/tau for p in points]
    offset=max(logits);weights=[mp.exp(a-offset) for a in logits];total=sum(weights)
    mean=[sum(w*p[j] for w,p in zip(weights,points))/total for j in range(3)]
    rb=[2*b-a for a,b in zip(x,mean)];radius=mp.sqrt(rb[0]**2+rb[1]**2)
    if not radius>0:raise ValueError('Nonunique circle projection')
    a=[rb[0]/radius,rb[1]/radius,mp.mpf(0)];ra=[2*c-d for c,d in zip(a,rb)]
    return [beta*(c+d)/2+(1-beta)*e for c,d,e in zip(ra,x,mean)],radius


def main(args):
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    args.out.mkdir(parents=True,exist_ok=False);rows=[]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=args.ids,digits=args.digits,warm=10000,observed=50000,windows=10,initial=[.317,-.211,.419],method='Exact binary64 input lifts. Posterior computed from inner-product logits; nested reflections; renormalized finite separation with no production Jacobian.',mpmath=mp.__version__)
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    with np.load(args.source) as f:raw={i:(f['params'][i].copy(),f['points'][i,:f['counts'][i]].copy()) for i in args.ids}
    for digits in args.digits:
        mp.mp.dps=digits;eps=mp.mpf(10)**(-30 if digits==70 else -40)
        for ident in args.ids:
            begin=time.time();beta,tau=[mp.mpf(float(v)) for v in raw[ident][0]];points=[[mp.mpf(float(v)) for v in p] for p in raw[ident][1]]
            x=[mp.mpf(v) for v in [.317,-.211,.419]];v=[mp.mpf(v) for v in [1.,.371,-.217]];length=mp.sqrt(sum(a*a for a in v));y=[a+eps*b/length for a,b in zip(x,v)]
            rates=[mp.mpf(0) for _ in range(10)];min_radius=mp.inf;sample=[]
            for it in range(60000):
                x,radius=iterate(x,beta,tau,points,mp);y,_=iterate(y,beta,tau,points,mp);min_radius=min(min_radius,radius)
                d=[b-a for a,b in zip(x,y)];length=mp.sqrt(sum(a*a for a in d))
                if not length>0:raise ValueError('Finite separation collapsed')
                if it>=10000:rates[(it-10000)//5000]+=mp.log(length/eps)
                y=[a+eps*b/length for a,b in zip(x,d)]
                if it>=50000:sample.append([float(v) for v in x])
            estimates=[v/5000 for v in rates];row=dict(id=ident,digits=digits,epsilon=str(eps),mean=float(sum(estimates)/10),minimum_window=float(min(estimates)),windows=[str(v) for v in estimates],minimum_radius=str(min_radius),seconds=time.time()-begin)
            rows.append(row);np.savez_compressed(args.out/f'case-{ident}-digits-{digits}.npz',states=np.array(sample),rates=np.array([float(v) for v in estimates]))
            (args.out/f'case-{ident}-digits-{digits}.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(row),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,scope='Finite-time sensitivity control. No proof of an asymptotic chaotic attractor or historical uniqueness.'),indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--digits',type=int,nargs='+',default=[70,100]);p.add_argument('--mpmath-path',type=Path);main(p.parse_args())
