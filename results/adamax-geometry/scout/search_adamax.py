"""Declared shape-first survey of scalar Adamax with its maximum memory switch."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from numba import set_num_threads
from adamax_map import scout,step_tangent


def controls():
    rng=np.random.default_rng(2026110123);errors=[];branches=[]
    for _ in range(400):
        p=np.r_[rng.uniform(.02,.99,2),5e-6];z=np.r_[rng.uniform(-2,2,2),rng.uniform(.02,3)]
        if abs(z[0])<1e-4 or abs(abs(z[0])+p[2]-p[1]*z[2])<1e-4:continue
        reset=abs(z[0])+p[2]>p[1]*z[2];sign=np.sign(z[0]);b,c,e=p
        def scalar(w):
            x,m,u=w;mn=b*m+(1-b)*x;un=sign*x+e if reset else c*u
            return np.array([x-mn/un,mn,un])
        val,j,actualbranch=step_tangent(z,np.eye(3),p)
        cs=np.column_stack([np.imag(scalar(z.astype(complex)+1e-25j*np.eye(3)[k]))/1e-25 for k in range(3)])
        error=float(np.max(abs(j-cs))/(1+np.max(abs(j))));assert error<1e-12
        assert np.max(abs(val-scalar(z)))<1e-12 and actualbranch==reset
        if reset:assert np.array_equal(j[:,2],np.zeros(3))
        else:assert abs(np.linalg.det(j)-b*c)<1e-10
        errors.append(error);branches.append(bool(reset))
    return dict(samples=len(errors),reset_branch_samples=sum(branches),maximum_scaled_jacobian_error=max(errors),
                method='Independent complex-step evaluation on the selected smooth branch; skip declared kink neighborhoods. Rank loss and memory-branch determinant checked.')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    seed=2026110122;u=qmc.Sobol(2,scramble=True,seed=seed).random_base2(12);p=np.c_[.02+.97*u[:,0],.02+.97*u[:,1],np.full(len(u),5e-6)]
    p[2048:,1]=1-10**(-4+3.95*u[2048:,1]);p=np.r_[p,[[.9,.999,5e-6],[.9,.99,5e-6],[.5,.9,5e-6],[.9,.5,5e-6],[.1,.1,5e-6]]]
    starts=np.tile([.73,0.,0.],(len(p),1));protocol=dict(seed=seed,count=len(p),alpha=.002,epsilon=1e-8,normalized_epsilon=5e-6,
        objective='theta^2/2',coordinates='(theta,M,U)/alpha',warmup=10000,observed=20000,window=5000,tail=256,
        design='4096Sobol rows: beta1 uniform[.02,.99); first2048beta2uniform[.02,.99),last2048beta2=1-10^(-4+3.95u); five explicit controls.',
        initialization='theta=.73*alpha,zero optimizer memories. Bias-corrected warmup, then autonomous observation.',
        primary_rule='All four leading windows>.005,completed,observed x range>1e-4.',
        source='Kingma&Ba2015Algorithm2; PyTorch2.15 Adamax epsilon inside gradient absolute value before maximum.',
        scope='Exploratory numerical shape search. Optimizer chaos is known; no new model or priority claim.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(out/'controls.json').write_text(json.dumps(controls(),indent=2)+'\n')
    for name in ['search_adamax.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    before=time.perf_counter();windows,tails,ends,status,extrema,resets,margins=scout(p,starts,10000,20000,5000,256)
    file=out/'scout.npz';np.savez_compressed(file,parameters=p,starts=starts,windows=windows,tails=tails,ends=ends,status=status,extrema=extrema,resets=resets,minimum_margins=margins)
    keep=(status==0)&(windows.min(1)>.005)&(extrema[:,1,0]-extrema[:,0,0]>1e-4)
    result=dict(completed=int((status==0).sum()),declared=len(p),primary_count=int(keep.sum()),primary_ids=np.flatnonzero(keep).tolist(),
        controls=[dict(id=int(i),parameters=p[i].tolist(),windows=windows[i].tolist(),extrema=extrema[i].tolist(),reset_fraction=float(resets[i]/20000)) for i in range(4096,len(p))],
        seconds=time.perf_counter()-before,sha256=hashlib.sha256(file.read_bytes()).hexdigest())
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['completed','declared','primary_count','controls','seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
