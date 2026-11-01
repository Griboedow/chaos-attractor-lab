"""Dense numerical portraits using the declared factorized posterior map."""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from numba import njit
from search_product_alphabet import step
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from render_portraits import render
from render_twist import camera


@njit(cache=True)
def orbit(beta,tau,levels,order,warm,count):
    x=np.array([.317,-.211,.419]);states=np.empty((count,3));minimum=np.inf
    for it in range(warm+count):
        x,J,r,peak,entropy,flag=step(x,beta,tau,levels,order)
        if flag or not np.isfinite(x).all():raise ValueError('Singular or nonfinite orbit')
        minimum=min(minimum,r)
        if it>=warm:states[it-warm]=x
    return states,minimum


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);args.images.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1];begin=time.time();rows=[]
    protocol=dict(ids=args.ids,warm=100000,count=args.count,initial=[.317,-.211,.419],camera='Right-handed PCA basis and one fixed orthogonal camera; equal units.',density='Every consecutive state; no interpolation or smoothing.')
    def save(path,data):path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    save(args.out/'protocol.json',protocol)
    for name in [Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for ident in args.ids:
        source=args.source/f'case-{ident}.npz'
        with np.load(source) as f:params=f['params'];levels=f['levels'];order=int(f['order']);points=f['points']
        z,minimum=orbit(float(params[0]),float(params[1]),levels,order,100000,args.count)
        center=z.mean(0);_,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
        if np.linalg.det(basis)<0:basis[:,-1]*=-1
        assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12
        projection=camera(.6,.4);xy=(z-center)@basis@projection.T;half=float(abs(xy).max()*1.06)
        h=np.histogram2d(xy[:,0],xy[:,1],bins=1600,range=[[-half,half]]*2)[0].T[::-1].astype(np.uint32)
        assert h.sum()==len(z)
        states=args.out/f'case-{ident}-states.npz';hist=args.out/f'case-{ident}-histogram.npz'
        np.savez_compressed(states,states=z,params=params,levels=levels,order=order,points=points)
        np.savez_compressed(hist,histogram=h)
        portrait=render(dict(id=f'product-alphabet-{ident}',panels=[dict(label=str(ident),mode='histogram',source=str(hist.resolve()),key='histogram')]),args.images)
        subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',portrait['file']],check=True)
        row=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),states_sha256=hashlib.sha256(states.read_bytes()).hexdigest(),histogram_sha256=hashlib.sha256(hist.read_bytes()).hexdigest(),center=center.tolist(),basis=basis.tolist(),projection=projection.tolist(),halfwidth=half,minimum_radius=minimum,portrait=portrait)
        rows.append(row);save(args.out/f'case-{ident}.json',row);print(json.dumps(dict(id=ident,seconds=time.time()-begin)),flush=True)
    save(args.out/'report.json',dict(protocol,cases=rows,seconds=time.time()-begin))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out','images']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--count',type=int,default=2000000);main(p.parse_args())
