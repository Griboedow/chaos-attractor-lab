"""Retained numerical screen for thinner stable widths at the same saddle."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from numba import njit
from search_raar_circles import step


@njit(cache=True)
def samples(c,B,r,p,q):
    Y=np.empty_like(q);axes=np.empty(len(q));S=B*r
    inverse=np.linalg.inv(S)
    for k in range(len(q)):
        x=c+S@q[k];axis=1e300
        for _ in range(7):
            x,J,gap,flag=step(x,p);assert flag==0;axis=min(axis,gap)
        Y[k]=inverse@(x-c);axes[k]=axis
    return Y,axes


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    raw=json.loads(args.input.read_text());c=np.array(raw['centers'][0]);B=np.array(raw['bases'][0]);r0=np.array(raw['radii'][0]);p=np.array(raw['params'])
    grid=np.linspace(-1,1,17);q=np.array([[u,v,z] for u in grid for v in grid for z in grid]);rows=[];arrays=[]
    for v in [.005,.002,.001,.0005,.0002,.0001,5e-5,2e-5,1e-5]:
        for z in [r0[2],1e-5,3e-6,1e-6,3e-7]:
            r=np.array([r0[0],v,z]);Y,axes=samples(c,B,r,p,q)
            stable=np.max(abs(Y[:,1:]),axis=1);required=abs(Y[:,0])<=1
            bound=float(stable[required].max()) if required.any() else 0.
            left=float(Y[q[:,0]==-1,0].min());right=float(-Y[q[:,0]==1,0].max())
            rows.append({'radii':r.tolist(),'maximum_required_stable':bound,'left_exit':left,'right_exit':right,'minimum_axis':float(axes.min()),'candidate':bound<.95 and min(left,right)>1.05})
            arrays.append(Y)
    candidates=[i for i,v in enumerate(rows) if v['candidate']]
    assert candidates
    chosen=min(candidates,key=lambda i:np.prod(rows[i]['radii']))
    refined=json.loads(json.dumps(raw));refined['radii'][0]=rows[chosen]['radii']
    (args.out/'input.json').write_text(json.dumps(refined,indent=2)+'\n')
    np.savez_compressed(args.out/'samples.npz',normalized=q,images=np.array(arrays))
    report={'original_input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'chosen':chosen,'trials':rows,'scope':'Numerical proposals only. The chosen thinner window still requires full interval self-covering and critical-face proofs.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    for name in [Path(__file__).name,'search_raar_circles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps({'chosen':chosen,'trial':rows[chosen],'candidate_count':len(candidates),'trials':len(rows)}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path('results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'));p.add_argument('--out',type=Path,required=True);main(p.parse_args())
