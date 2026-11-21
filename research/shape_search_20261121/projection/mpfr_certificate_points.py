"""Independent MPFR checks of certificate centers and stored inverse Jacobians."""
import argparse
import json
from pathlib import Path
import gmpy2 as g
from mpfr_check import build_map


def vec(v):return [g.mpfr(x) for x in v]


def apply(step,x,p):
    for _ in range(p):x=step(x)
    return x


def mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def inverse(a):
    x=[list(v)+[g.mpfr(i==j) for j in range(3)] for i,v in enumerate(a)]
    for col in range(3):
        pivot=max(range(col,3),key=lambda j:abs(x[j][col]))
        x[col],x[pivot]=x[pivot],x[col]
        val=x[col][col];x[col]=[v/val for v in x[col]]
        for j in range(3):
            if j==col:continue
            t=x[j][col];x[j]=[v-t*w for v,w in zip(x[j],x[col])]
    return [r[3:] for r in x]


def residual(step,doc,target,identity,basis=None):
    x=vec(doc['center']);p=doc['period'];h=g.mpfr(2)**-160
    y=apply(step,x,p)
    jac=[[g.mpfr(0)]*3 for _ in range(3)]
    for j in range(3):
        lo=x.copy();hi=x.copy();lo[j]-=h;hi[j]+=h
        a=apply(step,lo,p);b=apply(step,hi,p)
        for i in range(3):jac[i][j]=(b[i]-a[i])/(2*h)-int(identity and i==j)
    if basis is not None:
        chart,inv=basis;jac=mm(inv,mm(jac,chart))
    prod=mm([vec(row) for row in doc['preconditioner']],jac)
    return max(abs(v-w) for v,w in zip(y,target)),max(sum(abs(z-int(i==j)) for j,z in enumerate(row)) for i,row in enumerate(prod))


def main(args):
    all_rows=[]
    with g.context(precision=512):
        for prefix in (Path(),Path('three')):
            base=args.source/prefix
            for group in ('cycles','sinks','snapback'):
                for path in sorted((base/group).glob('*.json')):
                    doc=json.loads(path.read_text())
                    if group=='snapback' and 'fixed' not in doc:continue
                    if group!='snapback' and not doc.get('cycles'):continue
                    case=doc['case'];rec=json.loads((base/f'case-{case}.json').read_text());step=build_map(rec['parameters'])
                    tests=[]
                    if group=='snapback':
                        root=vec(doc['fixed']['center'])
                        tests.append(residual(step,doc['fixed'],root,True))
                        b=[vec(row) for row in doc['basis']];bi=inverse(b)
                        target=root
                        for box in doc['blocks']:
                            tests.append(residual(step,box,target,False,(b,bi)))
                            target=vec(box['center'])
                    else:
                        for cert in doc['cycles']:
                            tests.append(residual(step,cert,vec(cert['center']),True))
                    all_rows.append(dict(group=(prefix/group).as_posix(),case=case,centers=len(tests),max_state_residual=float(max(t[0] for t in tests)),max_inverse_jacobian_residual=float(max(t[1] for t in tests))))
    report=dict(bits=512,mpfr=g.mpfr_version(),gmpy2=g.version(),difference_step='2^-160',method='Independent numerical control, using centered finite differences of direct MPFR orbits and an independent matrix inversion. The interval proofs are rechecked separately.',rows=all_rows)
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(len(all_rows),sum(r['centers'] for r in all_rows),max(r['max_state_residual'] for r in all_rows),max(r['max_inverse_jacobian_residual'] for r in all_rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
