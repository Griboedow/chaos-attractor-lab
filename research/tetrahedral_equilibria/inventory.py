"""Numerical equilibrium inventory; explicitly not an exhaustive root certificate."""
from pathlib import Path
import argparse,json,sys
import numpy as np
from scipy.optimize import root,brentq
from scipy.stats import qmc
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'unexpected_search'))
from tetra_spectra import field

parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
base=Path(__file__).resolve().parents[2]/'results/unexpected-search';d=np.load(base/'tetra-v1/data.npz');out=args.output;out.mkdir(parents=True,exist_ok=False);rows=[]
for ident in [421,565]:
    c=d['centers'][ident];u=d['axes'][ident];s=d['sigma'][ident];g=d['gamma'][ident]
    def fun(z):return field(z,c,u,s,g)[0]
    def jac(z):return field(z,c,u,s,g)[1]
    seeds=list(6*qmc.Sobol(3,scramble=True,seed=2026110144+ident).random_base2(10)-3);axisroots=[]
    for axis in [np.array([1.,0.,0.]),np.ones(3)/np.sqrt(3)]:
        xs=np.linspace(-3,3,2001);values=np.array([np.dot(axis,fun(x*axis)) for x in xs])
        for j in range(len(xs)-1):
            if values[j]*values[j+1]<0:
                x=brentq(lambda a:np.dot(axis,fun(a*axis)),xs[j],xs[j+1],xtol=1e-14);z=x*axis;axisroots.append(z);seeds.append(z)
    roots=[np.zeros(3)];residuals=[];solved=0
    for seed in seeds:
        result=root(fun,seed,jac=jac,tol=1e-11)
        if np.linalg.norm(fun(result.x))<1e-10 and np.linalg.norm(result.x)<3:
            solved+=1
            for matrix in d['group']:
                z=matrix@result.x
                if min(np.linalg.norm(z-a) for a in roots)>1e-6:roots.append(z)
    points=np.load(base/f'tetra-render/sample-{ident}.npz')['points'];tree=cKDTree(points)
    info=[]
    for z in sorted(roots,key=lambda z:(round(np.linalg.norm(z),7),*np.round(z,7))):
        eigen=np.linalg.eigvals(jac(z));info.append({'point':z.tolist(),'radius':float(np.linalg.norm(z)),'residual':float(np.linalg.norm(fun(z))),'eigenvalues':[[float(v.real),float(v.imag)] for v in eigen],'positive_real_count':int((eigen.real>1e-8).sum()),'nearest_cloud_distance':float(tree.query(z)[0])})
    row={'id':ident,'seeds':len(seeds),'root_solver_accepted':solved,'roots':info,'axis_root_count':len(axisroots),'not_exhaustive':True};rows.append(row)
    print(json.dumps({'id':ident,'roots':len(roots),'radii':sorted(set(round(np.linalg.norm(z),7) for z in roots)),'positive_counts':{str(i):sum(r['positive_real_count']==i for r in info) for i in range(4)}}),flush=True)
(out/'inventory.json').write_text(json.dumps(rows,indent=2))
