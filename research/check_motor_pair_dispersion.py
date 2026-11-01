"""Finite close-pair dispersion: compare the motor braid with its larger regime.

Inspired by Wernecke et al.2017. This is a declared finite ensemble, not their
10,000-pair test and not a proof of partially predictable chaos or novelty.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from motor_rows import rhs
from motor_rows_gpu import MotorRowsGPU
from audit_motor_permutation_geometry import full


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);began=time.perf_counter()
    rng=np.random.default_rng(2026110106);jobs=[];starts=[];parameters=[];source_records=[]
    for source in args.sources:
        data=np.load(source);states=data['states'];p=data['parameters']
        weights=np.r_[np.ones(8),np.full(4,p[0]/p[1])]/np.sqrt(1+p[2]**2)
        physical=full(states,4)*weights;center=physical.mean(0);extent=float(np.sqrt(np.mean(np.sum((physical-center)**2,axis=1))))
        ids=np.rint(np.linspace(.1,.9,32)*(len(states)-1)).astype(int)
        # The same orbit phases and random directions are used at every delta.
        directions=rng.normal(size=(len(ids),11));directions/=np.linalg.norm(full(directions,4)*weights,axis=1)[:,None]
        source_records.append(dict(source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),parameters=p.tolist(),
            center=center.tolist(),extent=extent,weights=weights.tolist(),sample_indices=ids.tolist()))
        for delta in [1e-8,1e-6,1e-4]:
            for k,index in enumerate(ids):
                base=states[index];partner=base+delta*extent*directions[k]
                starts.extend([base,partner]);parameters.extend([p,p]);jobs.append(dict(source_index=len(source_records)-1,phase_index=k,delta=delta))
    starts=np.array(starts);parameters=np.array(parameters);state=starts.copy();gpu=MotorRowsGPU(4);tangent=None
    times=np.r_[0.,1.,2.,5.,10.,np.arange(20,1001,10),np.arange(1050,10001,50)]
    protocol=dict(sources=source_records,jobs=jobs,seed=2026110106,pairs=len(jobs),times=times.tolist(),
        metric='Full12D physical(c,s,X) metric from previous symmetry audit; perturbation norm=delta*source RMS extent.',
        selection='32 evenly spaced indices between10%and90%of each source history; shared phases/directions across three deltas.',
        method='GPU exact constrained11D field,rtol2e-11,step_factor.25; endpoints retained at every declared time. Tangents carried but are not pair-distance estimates.',
        correlation='Mean centered inner product divided by source extent squared; also normalized mean squared separation and mean distance. Source-time mean is an empirical finite baseline.',
        prior='https://doi.org/10.1038/s41598-017-01083-x',
        limitation='96 pairs per model,32 per delta; finite sampling and integration error. No asymptotic scaling, natural-measure, PPC or novelty proof.')
    write(out/'protocol.json',protocol);np.savez_compressed(out/'declared-inputs.npz',initial=starts,parameters=parameters)
    for name in ['check_motor_pair_dispersion.py','motor_rows_gpu.py','motor_rows_gpu.cl','motor_rows.py','audit_motor_permutation_geometry.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    # One representative pair per model at an early time, independently on CPU.
    control=gpu.scan(parameters[::192],starts[::192],transient=0,duration=2,rtol=2e-11,step_factor=.25,time_mode='original',max_wall=60)
    errors=[]
    for k,(p,q) in enumerate(zip(parameters[::192],starts[::192])):
        cpu=solve_ivp(lambda t,x:rhs(t,x,p,4),(0,2),q,method='Radau',rtol=2e-12,atol=2e-14,max_step=.005)
        assert cpu.success and control['status'][k]==0
        error=float(np.max(abs(cpu.y[:,-1]-control['final_states'][k])));assert error<1e-8;errors.append(error)
    write(out/'independent-controls.json',dict(two_tau_Radau_state_errors=errors))
    history=[state.copy()];metadata=[]
    for k in range(1,len(times)):
        result=gpu.scan(parameters,state,transient=0,duration=float(times[k]-times[k-1]),rtol=2e-11,step_factor=.25,
            time_mode='original',tangents=tangent,max_wall=120,max_steps=2000000)
        if np.any(result['status']!=0):
            write(out/'failure.json',dict(time=float(times[k]),statuses=result['status'].tolist()));return
        state=result['final_states'];tangent=result['final_tangents'];history.append(state.copy())
        metadata.append(dict(end=float(times[k]),accepted=int(result['step_counts'][:,0].sum()),rejected=int(result['step_counts'][:,1].sum())))
        if k%25==0 or k==len(times)-1:
            np.savez_compressed(out/'histories.npz',time=times[:k+1],states=np.array(history),parameters=parameters)
            write(out/'progress.json',dict(completed_time=float(times[k]),declared_end=float(times[-1]),blocks=metadata,wall_seconds=time.perf_counter()-began))
            print('pair observation',times[k],flush=True)
    histories=np.array(history);records=[]
    for si,s in enumerate(source_records):
        weights=np.array(s['weights']);center=np.array(s['center']);extent=s['extent']
        for delta in [1e-8,1e-6,1e-4]:
            ids=np.array([i for i,j in enumerate(jobs) if j['source_index']==si and j['delta']==delta])
            a=full(histories[:,2*ids],4)*weights;b=full(histories[:,2*ids+1],4)*weights
            distance=np.linalg.norm(a-b,axis=2)/extent
            correlation=np.sum((a-center)*(b-center),axis=2)/extent**2
            records.append(dict(source_index=si,delta=delta,mean_distance=distance.mean(1).tolist(),
                median_distance=np.median(distance,axis=1).tolist(),mean_square_distance=np.mean(distance**2,axis=1).tolist(),
                cross_correlation=correlation.mean(1).tolist(),correlation_standard_error=(correlation.std(1,ddof=1)/np.sqrt(len(ids))).tolist()))
    write(out/'analysis.json',dict(times=times.tolist(),records=records,wall_seconds=time.perf_counter()-began,
        history_sha256=hashlib.sha256((out/'histories.npz').read_bytes()).hexdigest(),
        uncertainty='Reported standard errors summarize these32 phases, not independent natural-measure draws or a confidence guarantee.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
