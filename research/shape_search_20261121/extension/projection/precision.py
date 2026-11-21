"""Independent finite-displacement spectrum in 256-bit MPFR arithmetic."""
import argparse
import json
from pathlib import Path
import sys
import gmpy2 as g

BASE=Path(__file__).resolve().parents[2]/'projection'
sys.path.insert(0,str(BASE))
from mpfr_check import add,norm


def metric_step(params):
    beta=g.mpfr(params['beta']); inv=g.mpfr(params['inv_tau']);tube=g.mpfr(params['tube_radius'])
    sc=[g.mpfr(v) for v in params['scales']];warp=g.mpfr(params.get('warp',0))
    sites=[list(map(g.mpfr,row)) for row in params['sites']]
    def step(x):
        logits=[-inv*sum((v-w)**2 for v,w in zip(x,b)) for b in sites];peak=max(logits)
        weights=[g.exp(v-peak) for v in logits];den=sum(weights)
        m=[sum(t*b[i] for t,b in zip(weights,sites))/den for i in range(3)]
        v=[2*a-b for a,b in zip(m,x)];h=[u/s for u,s in zip(v,sc)]
        h[2]-=warp*(h[0]**2-h[1]**2)
        rho=g.sqrt(h[0]**2+h[1]**2);center=[h[0]/rho,h[1]/rho,0]
        diff=[a-b for a,b in zip(h,center)];l=norm(diff);p=[a+tube*b/l for a,b in zip(center,diff)]
        p[2]+=warp*(p[0]**2-p[1]**2);out=[a*b for a,b in zip(sc,p)]
        return tuple(beta*(2*a-b+c)/2+(1-beta)*d for a,b,c,d in zip(out,v,x,m))
    return step


def check(record,bits,burn,window):
    with g.context(precision=bits):
        fn=metric_step(record['parameters']);eps=g.mpfr(2)**-(bits//3)
        x=tuple(map(g.mpfr,record['initial']));q=[tuple(g.mpfr(i==j) for i in range(3))for j in range(3)]
        rates=[[g.mpfr(0)]*3 for _ in range(4)]
        for t in range(0,burn+4*window,4):
            ys=[add(x,[eps*v for v in vec]) for vec in q]
            for _ in range(4):x=fn(x);ys=[fn(y) for y in ys]
            q=[]
            for k,y in enumerate(ys):
                v=[a-b for a,b in zip(y,x)]
                for _ in range(2):
                    for direction in q:
                        dot=sum(a*b for a,b in zip(v,direction));v=[a-dot*b for a,b in zip(v,direction)]
                size=norm(v);q.append([a/size for a in v])
                if t>=burn:rates[(t-burn)//window][k]+=g.log(size/eps)/window
        return dict(case=record['case'],bits=bits,burn=burn,window=window,windows=[[float(v) for v in row] for row in rates],shift=f'2^-{bits//3}',stride=4)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--records',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--cases',type=int,nargs='+',required=True);p.add_argument('--bits',type=int,default=256);p.add_argument('--burn',type=int,default=8192);p.add_argument('--window',type=int,default=8192)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    for k in a.cases:
        dest=a.out/f'mpfr-{k}.json'
        if not dest.exists():dest.write_text(json.dumps(check(json.loads((a.records/f'case-{k}.json').read_text()),a.bits,a.burn,a.window),indent=2)+'\n',encoding='utf8')
