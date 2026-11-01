"""Follow linear unstable seeds near the closest numerically found fixed points.

These finite clouds are not validated invariant manifolds or homoclinic proofs.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.spatial import cKDTree
from twist_map import step_tangent,spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    skeleton=json.loads(args.skeleton.read_text());long=np.load(args.reference)
    protocol=dict(skeleton=str(args.skeleton),skeleton_sha256=hashlib.sha256(args.skeleton.read_bytes()).hexdigest(),
        reference=str(args.reference),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        selection='For each1709/3254, choose the numerical fixed point closest to the archived reference cloud among all found roots.',
        radii=[1e-6,1e-8],starts_per_radius=128,warmup=4096,observed=50000,window=10000,tail=256,
        one_dimensional='64 equally spaced positive magnitudes in[epsilon/4,epsilon]and their negatives.',
        two_dimensional='128 equally spaced angles on a circle in an orthonormal real basis of the complex unstable eigenspace.',
        scope='Linear unstable seeds, not exact manifold points. Finite approach does not prove a homoclinic connection or attractor classification.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_twist_unstable_seeds.py','twist_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[]
    for ident in [1709,3254]:
        row=next(x for x in skeleton['runs'] if x['id']==ident);p=np.array(row['parameters'])
        fixed=min((x for x in row['periodic_orbits'] if x['least_period']==1),key=lambda x:min(x['nearest_reference_distance']))
        z=np.array(fixed['states'][0]);_,j=step_tangent(z,np.eye(3),p);values,vectors=np.linalg.eig(j)
        unstable=np.flatnonzero(abs(values)>1+1e-7)
        if len(unstable)==1:
            basis=vectors[:,unstable[0]].real[:,None];basis/=np.linalg.norm(basis)
            mag=np.linspace(.25,1.,64);unit=np.r_[mag,-mag][:,None]@basis.T
        else:
            assert len(unstable)==2 and abs(values[unstable[0]].imag)>.01
            v=vectors[:,unstable[0]];basis=np.linalg.qr(np.c_[v.real,v.imag])[0]
            theta=2*np.pi*np.arange(128)/128;unit=np.c_[np.cos(theta),np.sin(theta)]@basis.T
        starts=np.concatenate([z+eps*unit for eps in protocol['radii']]);pars=np.tile(p,(len(starts),1))
        windows,tail,ends,means,radii=spectra(pars,starts,4096,50000,10000,256)
        ref=long['tail'][list(long['ids']).index(ident)];dist=cKDTree(ref).query(tail.reshape(-1,3),workers=4)[0].reshape(len(starts),-1)
        residual=np.max(np.linalg.norm(tail[:,3:]-tail[:,:-3],axis=2),axis=1)
        file=out/f'case-{ident}.npz';np.savez_compressed(file,parameters=p,fixed=z,basis=basis,starts=starts,windows=windows,tail=tail,ends=ends,means=means,radii=radii,distance=dist,period3_residual=residual)
        result=dict(id=ident,sha256=hashlib.sha256(file.read_bytes()).hexdigest(),fixed_point=fixed,groups=[])
        for k,eps in enumerate(protocol['radii']):
            sl=slice(k*128,(k+1)*128);positive=windows[sl,:,0].min(1)>0;periodic=residual[sl]<1e-8
            result['groups'].append(dict(radius=eps,positive_all_windows=int(positive.sum()),period3=int(periodic.sum()),
                unclassified=int((~positive&~periodic).sum()),minimum_leading_window=float(windows[sl,:,0].min()),
                positive_reference_distance_quantiles=np.quantile(dist[sl][positive],[.5,.95,.99,1]).tolist() if positive.any() else []))
        rows.append(result);(out/'analysis.json').write_text(json.dumps(dict(runs=rows),indent=2)+'\n');print(ident,result['groups'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--skeleton',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
