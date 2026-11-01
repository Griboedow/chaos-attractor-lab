"""Nonexhaustive equilibrium inventory for the two icosahedral geometries."""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import root
from scipy.stats import qmc
from scipy.spatial import cKDTree
from tetra_spectra import field

base=Path(__file__).parent;d=np.load(base/'icosa-v1/input.npz');out=base/'icosa-equilibria';out.mkdir(exist_ok=False);rows=[]
for ident,radius in [(177,1.8),(183,3.3)]:
    c=d['centers'][ident];u=d['axes'][ident];s=d['sigma'][ident];g=d['gamma'][ident]
    def fun(z):return field(z,c,u,s,g)[0]
    def jac(z):return field(z,c,u,s,g)[1]
    seeds=2*radius*(qmc.Sobol(3,scramble=True,seed=2026110152+ident).random_base2(10)-.5)
    cloud=np.load(base/f'icosa-long/{ident}-dense.npy',mmap_mode='r')[::16];seeds=np.vstack([seeds,cloud[::160]])
    roots=[];solved=0;begin=time.time()
    for seed in seeds:
        result=root(fun,seed,jac=jac,tol=1e-11)
        if np.linalg.norm(fun(result.x))<1e-10 and np.linalg.norm(result.x)<radius:
            solved+=1
            if not roots or min(np.linalg.norm(result.x-a) for a in roots)>1e-6:
                for matrix in d['group']:
                    z=matrix@result.x
                    if not roots or min(np.linalg.norm(z-a) for a in roots)>1e-6:roots.append(z)
    tree=cKDTree(cloud);info=[]
    for z in sorted(roots,key=lambda z:(round(np.linalg.norm(z),7),*np.round(z,7))):
        ev=np.linalg.eigvals(jac(z));info.append({'point':z.tolist(),'radius':float(np.linalg.norm(z)),'residual':float(np.linalg.norm(fun(z))),'eigenvalues':[[float(v.real),float(v.imag)] for v in ev],'positive_real_count':int((ev.real>1e-8).sum()),'nearest_cloud_distance':float(tree.query(z)[0])})
    row={'id':ident,'seeds':len(seeds),'root_solver_accepted':solved,'roots':info,'not_exhaustive':True,'seconds':time.time()-begin};rows.append(row)
    print(json.dumps({'id':ident,'roots':len(roots),'seconds':row['seconds'],'radii':sorted(set(round(np.linalg.norm(z),7) for z in roots)),'positive_counts':{str(i):sum(r['positive_real_count']==i for r in info) for i in range(4)}}),flush=True)
(out/'inventory.json').write_text(json.dumps(rows,indent=2))
