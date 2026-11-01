"""Propose real-root intervals for a complete inverse-image classification."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as P


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);d=json.loads(args.source.read_text());rows=[]
    for case in d['targets']:
        beta,r,t,cx,cy,cz=case['params'];yx,yy,yz=case['target'];ct=np.cos(t);st=np.sin(t);D=P([1,0,1]);S=P([0,2]);C=P([1,0,-1]);X=(cx-yx)*D+r*C;Y=(cy-yy)*D+r*ct*S
        f=-(yx-(1-beta)*cx)*S+(ct*(yy-(1-beta)*cy)+st*(yz-(1-beta)*cz))*C
        e=-(cx-yx)*S*D+ct*(cy-yy)*C*D-r*st*st*S*C
        polynomial=f*f*(X*X+Y*Y)-beta*beta*e*e;roots=polynomial.roots();real=sorted(z.real for z in roots if abs(z.imag)<1e-7);intervals=[[float(z-1e-7*(1+abs(z))),float(z+1e-7*(1+abs(z)))] for z in real]
        assert all(b<c for (a,b),(c,d) in zip(intervals,intervals[1:]));rows.append(dict(label=f'target-{case["index"]}',params=case['params'],target=case['target'],root_intervals=intervals,proposal_coefficients=polynomial.coef.tolist(),proposal_complex_roots=[[float(z.real),float(z.imag)] for z in roots]))
    report=dict(cases=rows,source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),scope='Binary64 polynomial roots propose intervals only; exact counts and radial signs must be independently validated.')
    (args.out/'input.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print([len(r['root_intervals']) for r in rows])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
