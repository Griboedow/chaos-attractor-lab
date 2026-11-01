"""Dense portraits of actual Gaussian-posterior iterates, without interpolation."""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from numba import njit
from search_soft_constellations import step
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from render_portraits import render
from render_twist import camera


@njit(cache=True)
def orbit(beta,tau,points,warm,count):
    x=np.array([.317,-.211,.419]);states=np.empty((count,3));min_radius=np.inf
    for it in range(warm+count):
        x,j,r,peak,entropy,flag=step(x,beta,points,tau)
        if flag or not np.isfinite(x).all():raise ValueError('Singular or nonfinite orbit')
        min_radius=min(min_radius,r)
        if it>=warm:states[it-warm]=x
    return states,min_radius


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);args.images.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1];rows=[];begin=time.time()
    protocol=dict(ids=args.ids,warm=100000,count=args.count,initial=[.317,-.211,.419],camera='Right-handed PCA basis followed by the same orthogonal camera for every case; equal spatial units.',density='All consecutive states; no smoothing, interpolation or synthetic samples.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'search_soft_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for ident in args.ids:
        source=args.source/f'case-{ident}.npz'
        with np.load(source) as f:params=f['params'];points=f['points'];kind=f['kind'];base_id=f['base_id']
        z,min_radius=orbit(float(params[0]),float(params[1]),points,100000,args.count)
        center=z.mean(0);_,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
        if np.linalg.det(basis)<0:basis[:,-1]*=-1
        assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12
        projection=camera(.6,.4);xy=(z-center)@basis@projection.T;half=float(abs(xy).max()*1.06)
        h=np.histogram2d(xy[:,0],xy[:,1],bins=1600,range=[[-half,half]]*2)[0].T[::-1].astype(np.uint32);assert h.sum()==len(z)
        states_path=args.out/f'case-{ident}-states.npz';hist_path=args.out/f'case-{ident}-histogram.npz'
        np.savez_compressed(states_path,states=z,params=params,points=points,kind=kind,base_id=base_id)
        np.savez_compressed(hist_path,histogram=h)
        portrait=render(dict(id=f'soft-constellation-{ident}',panels=[dict(label=str(ident),mode='histogram',source=str(hist_path.resolve()),key='histogram')]),args.images)
        subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',portrait['file']],check=True)
        record=dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),states_sha256=hashlib.sha256(states_path.read_bytes()).hexdigest(),histogram_sha256=hashlib.sha256(hist_path.read_bytes()).hexdigest(),center=center.tolist(),basis=basis.tolist(),projection=projection.tolist(),halfwidth=half,minimum_radius=min_radius,portrait=portrait)
        rows.append(record);(args.out/f'case-{ident}.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(json.dumps(dict(id=ident,seconds=time.time()-begin)),flush=True)
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,seconds=time.time()-begin),indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--images',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--count',type=int,default=2000000);main(p.parse_args())
