"""Independent full QR growth from high-precision finite separations.

Uses scalar nested reflections and no production Jacobian. This is a numerical
precision control, not interval validation or an asymptotic Lyapunov theorem.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from check_soft_sphere_precision import iterate


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    args.out.mkdir(parents=True,exist_ok=False)
    for name in [Path(__file__).name,'check_soft_sphere_precision.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:
        raw={i:(f['params'][i].copy(),f['points'][i,:f['counts'][i]].copy()) for i in args.ids}
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=args.ids,digits=args.digits,warm=10000,observed=100000,windows=20,initial=[.317,-.211,.419],method='Three renormalized finite separations. Modified Gram-Schmidt. Independent scalar nested-reflection formula; no production Jacobian.',input='Exact lifts of the stored binary64 inputs',mpmath=mp.__version__)
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    rows=[]
    for digits in args.digits:
        mp.mp.dps=digits;eps=mp.mpf(10)**(-30 if digits==70 else -40)
        for ident in args.ids:
            begin=time.time();beta,tau=[mp.mpf(float(v)) for v in raw[ident][0]];points=[[mp.mpf(float(v)) for v in p] for p in raw[ident][1]]
            x=[mp.mpf(v) for v in [.317,-.211,.419]]
            Q=[[mp.mpf(i==j) for i in range(3)] for j in range(3)]
            sums=[[mp.mpf(0) for _ in range(3)] for _ in range(20)]
            min_radius=mp.inf;orth_error=mp.mpf(0)
            for it in range(110000):
                y,r=iterate(x,beta,tau,points,mp);min_radius=min(min_radius,r)
                vectors=[]
                for q in Q:
                    yp,_=iterate([a+eps*b for a,b in zip(x,q)],beta,tau,points,mp)
                    vectors.append([(a-b)/eps for a,b in zip(yp,y)])
                Q=[]
                for j,v in enumerate(vectors):
                    for q in Q:
                        dot=sum(a*b for a,b in zip(q,v));v=[a-dot*b for a,b in zip(v,q)]
                    length=mp.sqrt(sum(a*a for a in v));assert length>0
                    Q.append([a/length for a in v])
                    if it>=10000:sums[(it-10000)//5000][j]+=mp.log(length)
                if it%5000==0:
                    orth_error=max(orth_error,max(abs(sum(a*b for a,b in zip(Q[i],Q[j]))-mp.mpf(i==j)) for i in range(3) for j in range(3)))
                x=y
            rates=[[v/5000 for v in row] for row in sums]
            record=dict(id=ident,digits=digits,epsilon=str(eps),spectrum=[float(sum(row[j] for row in rates)/20) for j in range(3)],minimum_windows=[float(min(row[j] for row in rates)) for j in range(3)],windows=[[str(v) for v in row] for row in rates],minimum_radius=str(min_radius),maximum_sampled_orthogonality_error=str(orth_error),seconds=time.time()-begin)
            rows.append(record);(args.out/f'case-{ident}-{digits}.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({k:v for k,v in record.items() if k!='windows'}),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,scope='Finite-time full-spectrum precision control. Two positive estimates do not prove asymptotic hyperchaos or novelty.'),indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--digits',type=int,nargs='+',default=[70,100]);p.add_argument('--mpmath-path',type=Path,required=True);main(p.parse_args())
