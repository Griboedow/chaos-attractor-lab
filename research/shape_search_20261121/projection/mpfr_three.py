"""Independent full spectrum from one MPFR orbit and three finite shadows."""
import argparse
import json
import os
from pathlib import Path
import gmpy2 as g

from mpfr_check import build_map,norm,add


def replay(record,bits,burn,window,windows,stride):
    if burn%stride or window%stride:raise ValueError('Window and burn must be multiples of stride')
    with g.context(precision=bits):
        step=build_map(record['parameters'])
        scale=g.mpfr(2)**-(bits//3)
        x=tuple(g.mpfr(v) for v in record['initial'])
        unit=[tuple(g.mpfr(j==i) for j in range(3)) for i in range(3)]
        sums=[[g.mpfr(0)]*3 for _ in range(windows)]
        for n in range(0,burn+window*windows,stride):
            shadows=[add(x,[scale*v for v in u]) for u in unit]
            for _ in range(stride):
                x=step(x)
                shadows=[step(y) for y in shadows]
            unit=[]
            for j,y in enumerate(shadows):
                vec=[v-w for v,w in zip(y,x)]
                for _ in range(2):
                    for u in unit:
                        c=sum(a*b for a,b in zip(u,vec))
                        vec=[v-c*w for v,w in zip(vec,u)]
                mag=norm(vec)
                unit.append([v/mag for v in vec])
                if n>=burn:sums[(n-burn)//window][j]+=g.log(mag/scale)/window
        return dict(case=record['case'],bits=bits,mpfr=g.mpfr_version(),gmpy2=g.version(),method='Four independent finite MPFR trajectories. Twice-reorthogonalized finite displacements and no map Jacobian.',shift_power_of_two=-(bits//3),burn=burn,window_length=window,qr_stride=stride,windows=[[float(v) for v in row] for row in sums],final_state=[str(v) for v in x])


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    paths=[args.source/f'case-{i}.json' for i in args.cases] if args.cases else sorted(args.source.glob('case-*.json'))
    for path in paths:
        rec=json.loads(path.read_text())
        out=args.out/f'mpfr-three-{rec["case"]}-{args.bits}.json'
        if out.exists():continue
        doc=replay(rec,args.bits,args.burn,args.window,args.windows,args.stride)
        tmp=out.with_suffix('.tmp.json')
        tmp.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
        os.replace(tmp,out)
        print(rec['case'],min(v[2] for v in doc['windows']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cases',type=int,nargs='+')
    p.add_argument('--bits',type=int,default=192)
    p.add_argument('--burn',type=int,default=8192)
    p.add_argument('--window',type=int,default=8192)
    p.add_argument('--windows',type=int,default=4)
    p.add_argument('--stride',type=int,default=4)
    main(p.parse_args())
