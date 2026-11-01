"""Independent arbitrary-precision nested-reflector and finite-separation control."""
import argparse,hashlib,json,time
from pathlib import Path
import mpmath as mp
import numpy as np


def iterate(x,beta,points):
    distances=[sum((a-b)**2 for a,b in zip(x,p)) for p in points];order=sorted(range(len(points)),key=lambda i:distances[i]);i=order[0];gap=distances[order[1]]-distances[i]
    if not gap>0:raise ValueError('Unresolved exact nearest-point tie')
    b=points[i];rb=[2*c-a for a,c in zip(x,b)];radius=mp.sqrt(rb[0]**2+rb[1]**2)
    if not radius>0:raise ValueError('Circle projection is not unique')
    a=[rb[0]/radius,rb[1]/radius,mp.mpf(0)];ra=[2*c-d for c,d in zip(a,rb)]
    return [beta*(c+d)/2+(1-beta)*e for c,d,e in zip(ra,x,b)],i,gap,radius


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=args.ids,digits=args.digits,warm=10000,observed=50000,windows=10,initial=[.317,-.211,.419],method='Nested reflections with scalar arbitrary-precision arithmetic. Renormalized finite separation, no production Jacobian. Exact lifts of binary64 beta and point coordinates.',mpmath=mp.__version__)
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    with np.load(args.source) as f:raw={i:(float(f['params'][i,0]),f['points'][i,:f['counts'][i]].copy()) for i in args.ids}
    for digits in args.digits:
        mp.mp.dps=digits;eps=mp.mpf(10)**(-30 if digits==70 else -40)
        for i in args.ids:
            begin=time.time();beta=mp.mpf(raw[i][0]);points=[[mp.mpf(float(v)) for v in p] for p in raw[i][1]];x=[mp.mpf(v) for v in [.317,-.211,.419]];v=[mp.mpf(v) for v in [1.,.371,-.217]];length=mp.sqrt(sum(a*a for a in v));y=[a+eps*b/length for a,b in zip(x,v)];rates=[mp.mpf(0) for _ in range(10)];mismatch=0;mingap=mp.inf;minaxis=mp.inf;visits=np.zeros(len(points),int);sample=[]
            for it in range(60000):
                x,k,gap,axis=iterate(x,beta,points);y,j,_,_=iterate(y,beta,points);mismatch+=int(k!=j);mingap=min(mingap,gap);minaxis=min(minaxis,axis);d=[b-a for a,b in zip(x,y)];length=mp.sqrt(sum(a*a for a in d))
                if not length>0:raise ValueError('Finite separation collapsed')
                if it>=10000:rates[(it-10000)//5000]+=mp.log(length/eps);visits[k]+=1
                y=[a+eps*b/length for a,b in zip(x,d)]
                if it>=50000:sample.append([float(v) for v in x])
            estimates=[v/5000 for v in rates];row=dict(id=i,digits=digits,epsilon=str(eps),mean=float(sum(estimates)/10),minimum_window=float(min(estimates)),windows=[str(v) for v in estimates],branch_disagreements=mismatch,minimum_gap=str(mingap),minimum_axis=str(minaxis),visits=visits.tolist(),seconds=time.time()-begin)
            rows.append(row);np.savez_compressed(args.out/f'case-{i}-digits-{digits}.npz',states=np.array(sample),rates=np.array([float(v) for v in estimates]),visits=visits);(args.out/f'case-{i}-digits-{digits}.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8');print(json.dumps(row),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,scope='Finite-time sensitivity control. Does not prove historical uniqueness or a chaotic invariant attractor of a discontinuous map.'),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--digits',type=int,nargs='+',default=[70,100]);main(p.parse_args())
