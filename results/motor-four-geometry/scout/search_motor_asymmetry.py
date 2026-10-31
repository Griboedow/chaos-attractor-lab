"""Declared N-row survey of the published motor transition-rate asymmetry."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from motor_rows import field_jac,source_force_constraint_field
from motor_rows_gpu import MotorRowsGPU,save,default_tangent


def write(p,x):
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    n=args.rows;d=3*n-1
    if n<3:raise ValueError('This survey requires at least three rows')
    rho=np.geomspace(.3,20,12)
    gain=np.geomspace(1.03,3.,24)
    asymmetry=np.array([0.,.05,.15,.35,.7,1.5,3.])
    tuples=np.array([[a,q*(a+1),b] for b in asymmetry for a in rho for q in gain])
    params=np.repeat(tuples,2,axis=0)
    initial=np.zeros((len(params),d));initial[:,:n]=1.;initial[:,n:2*n]=params[:,2,None]
    starts=np.random.default_rng(2026103142).uniform(-.05,.05,(2,n))
    starts-=starts.mean(axis=1)[:,None]
    initial[:,2*n:]=np.tile(starts[:,:n-1],(len(tuples),1))
    np.savez_compressed(out/'declared-inputs.npz',parameters=params,initial=initial,
                        tuple_id=np.repeat(np.arange(len(tuples)),2),start_id=np.tile([0,1],len(tuples)))
    write(out/'protocol.json',dict(links=n,dimension=d,parameters='rho,r,b',rho=rho.tolist(),
        relative_gain=gain.tolist(),asymmetry=asymmetry.tolist(),count=len(params),seed=2026103142,
        initial='c_j=1, s_j=b, two zero-sum position vectors drawn uniformly then centered; no historical endpoint.',
        source='Alouges et al., arXiv:2409.03506v2, Eq.15; exact first-harmonic closure.',
        scope='Dimensionless family with varied prescribed transition-rate asymmetry; not a calibrated experimental activity path.',
        invariant='c_j^2+s_j^2 <= 1+b^2; choose A=-.1/sqrt(1+b^2) for nonnegative source rates.',
        observation=args.duration,warm=args.warm,rtol=1e-9,max_wall=args.wall,
        geometry_goal='Inspect literal state coordinates and declared linear physical observables. No inferred topological novelty.'))
    for name in ['search_motor_asymmetry.py','motor_rows.py','motor_rows_gpu.py','motor_rows_gpu.cl']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=MotorRowsGPU(n)
    controls=[]
    for b in [0.,.35,3.]:
        p=np.array([3.7,6.,b]);x=np.r_[np.ones(n),np.full(n,b),starts[0,:n-1]];v=default_tangent(n)
        f,j=field_jac(x,p,n)
        field_error=float(max(abs(f-source_force_constraint_field(x,p,n))))
        def augmented(t,z):
            f,j=field_jac(z[:d],p,n)
            return np.r_[f,j@z[d:]]
        sol=solve_ivp(augmented,(0,2),np.r_[x,v],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.005)
        assert sol.success
        pilot=gpu.scan(p[None],x[None],transient=0,duration=2,rtol=1e-11,step_factor=.5,
                       time_mode='original',max_wall=30,tangents=v[None])
        norm=np.linalg.norm(sol.y[d:,-1])
        error=float(max(abs(sol.y[:d,-1]-pilot['final_states'][0])))
        tangent_error=float(max(abs(sol.y[d:,-1]/norm-pilot['final_tangents'][0])))
        rate_error=float(abs(np.log(norm)/2-pilot['lle'][0]))
        assert pilot['status'][0]==0 and max(error,tangent_error,field_error)<1e-8 and rate_error<1e-7
        controls.append(dict(parameters=p.tolist(),field_error=field_error,state_error=error,
                             tangent_error=tangent_error,rate_error=rate_error))
    write(out/'independent-controls.json',controls)
    result=gpu.scan(params,initial,transient=args.warm,duration=args.duration,rtol=1e-9,
                    time_mode='original',step_factor=.5,max_steps=15000000,max_wall=args.wall)
    result.update(tuple_id=np.repeat(np.arange(len(tuples)),2),start_id=np.tile([0,1],len(tuples)),
                  protocol_sha256=hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest())
    save(out/'complete-or-budgeted.npz',result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rows',type=int,default=3)
    p.add_argument('--duration',type=float,default=1200.)
    p.add_argument('--warm',type=float,default=300.)
    p.add_argument('--wall',type=float,default=1200.)
    main(p.parse_args())
