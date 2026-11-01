"""Recheck the entire declared polymer scout in positive logarithmic coordinates."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from search_rho_adsi import GPU,write
from supramolecular import REFERENCE,equilibrium
from supramolecular_log import field_jac,independent_field,from_raw,to_raw


def controls(gpu,out):
    records=[]
    for p in [REFERENCE,np.array([1.,200.,8.,4.6]),np.array([.03,5000.,3.,4.6])]:
        raw=equilibrium(p)['state']*np.array([1.001,.999,1.001]);q=from_raw(raw,p)
        f,j=field_jac(q,p);jj=np.column_stack([independent_field(q+1e-25j*np.eye(3)[i],p).imag/1e-25 for i in range(3)])
        algebra=float(max(np.max(abs(f-independent_field(q,p))/(1+abs(f))),np.max(abs(j-jj)/(1+abs(jj)))))
        assert algebra<1e-10
        v=np.array([1.,2.,3.])/np.sqrt(14.)
        def augmented(t,z):
            f,j=field_jac(z[:3],p);return np.r_[f,j@z[3:]]
        cpu=solve_ivp(augmented,(0,2),np.r_[q,v],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
        assert cpu.success
        result=gpu.scan(p[None],q[None],0,2,rtol=1e-11,maxstep=.005)
        norm=np.linalg.norm(cpu.y[3:,-1]);errors=dict(algebra=algebra,state=float(np.max(abs(result['states'][0,:3]-cpu.y[:3,-1]))),
            tangent=float(np.max(abs(result['states'][0,3:6]-cpu.y[3:,-1]/norm))),rate=float(abs(result['rate'][0]-np.log(norm)/2)))
        assert result['status'][0]==0 and max(errors.values())<1e-7
        records.append(dict(parameters=p.tolist(),errors=errors))
    write(out/'independent-controls.json',records)


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source/'declared-inputs.npz';data=np.load(source);params=data['parameters'];raw=data['initial']
    starts=np.array([from_raw(x,p) for x,p in zip(raw,params)])
    write(out/'protocol.json',dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        design='Identical parameter tuples and physical starts as the raw-coordinate scout; no new coverage.',lanes=len(params),
        state='x=a1/kd2,u=log(m0/kd3),v=log(L-2), L=polymer mass/polymer number; no concentration clipping.',
        warm=args.warm,observation=args.duration,rtol=args.rtol,atol=args.rtol*.001,maxstep=.02,
        method='Adaptive FP64 Dormand-Prince5(4), exact transformed tangent, full tangent-norm error scale.',
        retention='All4 log-coordinate leading windows>.005; exploratory>.001. Minimum observed L>10 is a separate approximation screen.',
        limitations='Approximate polymer moment closure. Finite-time exponents change under nonlinear coordinates; no chaos or priority claim. Very low polymer populations can invalidate the physical closure even when numerically resolved.'))
    np.savez_compressed(out/'declared-inputs.npz',parameters=params,initial=starts,physical_initial=raw)
    for name in ['recheck_supramolecular_log.py','supramolecular_log.py','supramolecular_log.cl','supramolecular.py','search_rho_adsi.py','rho_adsi.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=GPU(kernel='supramolecular_log.cl',field=field_jac);controls(gpu,out)
    result=gpu.scan(params,starts,args.warm,args.duration,rtol=args.rtol,maxstep=.02)
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    primary=(result['status']==0)&(result['windows'].min(axis=1)>.005);valid=primary&(result['stats'][:,7]>10)
    weak=(result['status']==0)&(result['windows'].min(axis=1)>.001)
    write(out/'analysis.json',dict(total=len(params),statuses={str(i):int((result['status']==i).sum()) for i in range(6)},
        primary_ids=np.flatnonzero(primary).tolist(),long_polymer_primary_ids=np.flatnonzero(valid).tolist(),
        exploratory_ids=np.flatnonzero(weak).tolist(),wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],
        sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest()))
    print('complete',int((result['status']==0).sum()),'primary',int(primary.sum()),'long-polymer primary',int(valid.sum()),'weak',int(weak.sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--warm',type=float,default=100.);p.add_argument('--duration',type=float,default=800.);p.add_argument('--rtol',type=float,default=1e-9)
    main(p.parse_args())
