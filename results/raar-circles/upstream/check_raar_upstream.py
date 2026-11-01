"""Run the pinned, unchanged ProxToolbox evaluator with circle adapters.

This tests the published evaluator, not a full image-reconstruction experiment.
The original Algorithm constructor and relaxation calculation are retained.
"""
import argparse,concurrent.futures,hashlib,json,shutil,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from numba import njit
from search_raar_circles import step
from follow_raar_circles import orthogonalize


class CircleA:
    def __init__(self,experiment):pass
    def eval(self,x):
        length=np.hypot(x[0],x[1])
        if length<1e-12:raise ValueError('Projection axis A')
        return np.array([x[0]/length,x[1]/length,0.])


class CircleB:
    def __init__(self,experiment):
        p=experiment.parameters;self.radius=p[1];self.center=p[3:]
        self.cos=np.cos(p[2]);self.sin=np.sin(p[2])
    def eval(self,x):
        # Independent orthonormal-plane basis, instead of I-n n^T.
        d=x-self.center;a=d[0];b=self.cos*d[1]+self.sin*d[2]
        length=np.hypot(a,b)
        if length<1e-12:raise ValueError('Projection axis B')
        return self.center+self.radius/length*np.array([a,self.cos*b,self.sin*b])


@njit(cache=True)
def analyze(states,params):
    q,_,_=orthogonalize(np.array([[1.,.2,-.3],[.4,1.,.15],[-.25,.35,1.]]))
    rates=np.zeros((29,3));detlogs=np.zeros(29);residual=0.;minaxis=np.inf
    for n in range(100000,400000):
        expected,j,axis,status=step(states[n],params)
        if status:raise ValueError('Transcription failed')
        minaxis=min(minaxis,axis);residual=max(residual,np.max(np.abs(expected-states[n+1])))
        q,r,status=orthogonalize(j@q)
        if status:raise ValueError('QR failed')
        if n>=110000:
            k=(n-110000)//10000;rates[k]+=r;detlogs[k]+=np.log(abs(np.linalg.det(j)))
    return rates/10000,detlogs/10000,residual,minaxis


def run(task):
    ident,params,outdir=task;out=Path(outdir);sys.path.insert(0,str(out/'upstream'))
    from proxtoolbox.algorithms.RAAR import RAAR
    exp=SimpleNamespace(proxOperators=[CircleA,CircleB],parameters=params,
        norm_data=1,Nx=3,Ny=1,Nz=1,product_space_dimension=1,truth=None,
        truth_dim=0,norm_truth=0,diagnostic=False,MAXIT=400000,TOL=0,
        lambda_0=params[0],lambda_max=params[0],lambda_switch=1,debug=False)
    solver=RAAR(exp,None);s=np.array([.317,-.211,.419]);states=np.empty((400001,3));states[0]=s
    start=time.time();relax_error=0.
    for n in range(400000):
        solver.iter=n;relax_error=max(relax_error,abs(solver.computeRelaxation()-params[0]))
        s=solver.evaluate(s);states[n+1]=s
    assert np.isfinite(states).all()
    rates,detlogs,residual,minaxis=analyze(states,params);assert residual<1e-10
    np.savez_compressed(out/f'case-{ident}.npz',states=states,params=params,rates=rates,detlogs=detlogs)
    row=dict(id=ident,params=params.tolist(),steps=400000,warm=100000,alignment=10000,observed=290000,windows=29,
        means=rates.mean(0).tolist(),minima=rates.min(0).tolist(),max_step_residual=residual,min_projection_axis_distance=minaxis,
        max_determinant_identity=float(np.max(abs(rates.sum(1)-detlogs))),max_relaxation_roundoff=relax_error,seconds=time.time()-start)
    (out/f'case-{ident}.json').write_text(json.dumps(row,indent=2),encoding='utf-8');return row


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);package=args.out/'upstream/proxtoolbox'
    for sub in ['', 'algorithms','utils']:
        dest=package/sub;dest.mkdir(parents=True,exist_ok=True);(dest/'__init__.py').write_text('')
    for a,b in [('upstream-RAAR.py','algorithms/RAAR.py'),('upstream-algorithm.py','algorithms/algorithm.py'),('upstream-cell.py','utils/cell.py')]:
        shutil.copy2(args.upstream/a,package/b)
    shutil.copy2(args.upstream/'upstream-LICENSE',args.out/'upstream/LICENSE')
    shutil.copy2(args.upstream/'repository-source.json',args.out/'upstream/source.json')
    (args.out/'upstream/ATTRIBUTION.txt').write_text('ProxToolbox, D. Russell Luke and contributors, University of Goettingen.\nhttps://gitlab.gwdg.de/nam/ProxPython\nPinned commit 702154fc33a0b05747ef7f2980c88613236e3626.\nThe three upstream Python files and LICENSE are copied without byte changes. Empty package initializers are supplied by this benchmark. Circle adapters and the experiment driver are original benchmark code. The repository LICENSE at this commit is Creative Commons Attribution 4.0 International; the older documentation says BSD. The pinned LICENSE is retained verbatim. No endorsement is implied.\n',encoding='utf-8')
    with np.load(args.source) as f:params=f['params']
    protocol=dict(ids=args.ids,source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        method='Unmodified upstream RAAR.evaluate and Algorithm constructor/computeRelaxation, using independent circle projector adapters. Full imaging package and monitor not invoked.',
        spectrum='Analytic transcription Jacobian evaluated along actual upstream states, not an independent tangent implementation.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    for name in ['check_raar_upstream.py','search_raar_circles.py','follow_raar_circles.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run,(i,params[i],str(args.out))) for i in args.ids]):
            row=future.result();rows.append(row);print(json.dumps(row),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,runs=rows),indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',default=[1051,1483,820,1679]);main(p.parse_args())
