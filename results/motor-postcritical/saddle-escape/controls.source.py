"""Follow both outgoing directions of a period-four motor saddle.

The pair (rows 2,3 in zero-based notation) is retained as an exact invariant
cluster while the unstable pair (6,7) splits. This is a numerical connection
experiment, not a rigorous proof of a heteroclinic orbit or network.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import linear_sum_assignment, minimize_scalar

from motor_postcritical_cycles import field_jac, lift, restrict, audit, ROOT
from motor_rows import rhs as full_rhs

GROUPS=[[0],[1],[2,3],[4],[5],[6],[7]]
WEIGHTS=np.array([len(q) for q in GROUPS])/8
M=len(GROUPS)
D=3*M-1


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def physical(full):
    full=np.asarray(full)
    x=np.concatenate([full[...,16:],-np.sum(full[...,16:],axis=-1,keepdims=True)],axis=-1)
    return np.stack([full[...,:8],full[...,8:16],x],axis=-1)


def batch_rhs(t,flat,p,count):
    state=flat.reshape(count,D)
    c,s=state[:,:M],state[:,M:2*M]
    x=np.c_[state[:,2*M:],-state[:,2*M:]@WEIGHTS[:-1]/WEIGHTS[-1]]
    v=-p[0]*x+p[1]*(s-(s@WEIGHTS)[:,None])
    return np.c_[1-c-v*s,p[2]-s+v*c,v[:,:-1]].ravel()


def run(source,output,periods=400,max_step=.05):
    source,output=source.resolve(),output.resolve()
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    data=np.load(source/'cycle.npz');p=data['parameters'];period=float(data['time'][-1]);base=data['full_states'][0]
    eig,vec=np.linalg.eig(data['full_monodromy']);index=np.argmax(abs(eig));vector=vec[:,index].real
    q=physical(vector);contrast=(q[6]-q[7])/2
    if contrast[0]<0:contrast*=-1
    delta=np.zeros((8,3));delta[6]=contrast;delta[7]=-contrast
    reduced=np.r_[delta[:,0],delta[:,1],delta[:-1,2]];reduced/=np.linalg.norm(reduced)
    starts=np.array([base+sign*eps*reduced for eps in (1e-7,1e-9) for sign in (1,-1)])
    quotient=np.array([restrict(z,GROUPS) for z in starts]);count=len(starts)
    independent=np.array([field_jac(z,p,WEIGHTS)[0] for z in quotient])
    batch_error=float(np.max(abs(independent-batch_rhs(0,quotient.ravel(),p,count).reshape(count,D))))
    assert batch_error<1e-12
    def one(t,z):return field_jac(z,p,WEIGHTS)[0]
    short=[]
    for method in ['DOP853','Radau']:
        s=solve_ivp(one,(0,period),quotient[0],method=method,rtol=2e-12,atol=2e-14,max_step=.0125)
        if not s.success:raise RuntimeError(s.message)
        short.append(s.y[:,-1])
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),parameters=p.tolist(),
        source=source.relative_to(ROOT).as_posix(),source_sha256=hashlib.sha256((source/'cycle.npz').read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),period=period,
        unstable_multiplier=float(eig[index].real),unstable_exponent=float(np.log(abs(eig[index]))/period),
        unstable_pair=[6,7],preserved_pair=[2,3],invariant_groups=GROUPS,
        initial_direction_full_reduced=reduced.tolist(),
        cases=[dict(index=i,epsilon=eps,sign=sign) for i,(eps,sign) in enumerate([(e,s) for e in (1e-7,1e-9) for s in (1,-1)])],
        periods_requested=periods,dt_sample=period/40,max_step=max_step,rtol=1e-11,atol=1e-13,
        algebra=audit(GROUPS,p),batch_field_error=batch_error,
        independent_DOP853_Radau_state_error=float(np.max(abs(short[0]-short[1]))),
        status='running',blocks=[],caveat='Exact invariant cluster restriction, approximate orbit and outgoing tangent. Numerical connection evidence only.')
    (output/'source.py').write_bytes(Path(__file__).read_bytes())
    began=time.perf_counter();elapsed=0.;alltimes=[0.];allstates=[starts]
    for k in range(0,periods,25):
        length=min(25,periods-k)*period
        times=np.linspace(0,length,min(25,periods-k)*40+1)
        sol=solve_ivp(lambda t,z:batch_rhs(t,z,p,count),(0,length),quotient.ravel(),method='DOP853',
                      rtol=1e-11,atol=1e-13,max_step=max_step,t_eval=times)
        if not sol.success:
            record.update(status='failed',error=sol.message);break
        block=sol.y.T.reshape(-1,count,D)
        full=np.array([[lift(z,GROUPS) for z in row] for row in block])
        alltimes.extend((elapsed+sol.t[1:]).tolist());allstates.extend(full[1:])
        quotient=block[-1];elapsed+=length
        record['blocks'].append(dict(end_tau=elapsed,completed_periods=k+min(25,periods-k),nfev=sol.nfev))
        record.update(elapsed_tau=elapsed,wall_seconds=time.perf_counter()-began)
        np.savez_compressed(output/'trajectory.npz',parameters=p,time=np.array(alltimes),
                            states=np.array(allstates),final_quotient_state=quotient,initial_full_state=starts)
        write(output/'summary.json',record)
        print('escape periods',k+min(25,periods-k),'wall',round(time.perf_counter()-began,1),flush=True)
    if record['status']=='running':record['status']='complete'
    write(output/'summary.json',record)


def analyze(source,folder):
    output=folder/'diagnostics.json'
    if output.exists():raise FileExistsError(output)
    a=np.load(folder/'trajectory.npz');meta=json.loads((folder/'summary.json').read_text());ref=np.load(source/'cycle.npz')
    qref=physical(ref['full_states']);qref[-1]=qref[0]
    period=float(ref['time'][-1]);spline=CubicSpline(ref['time'],qref,bc_type='periodic')
    coarse=np.arange(80)*period/80;coarse_q=spline(coarse)
    full=a['states'];times=a['time'];q=physical(full);records=[]
    for case in range(full.shape[1]):
        distance=np.linalg.norm(q[:,case,6]-q[:,case,7],axis=1)
        strobes=np.arange(0,len(times),40);sd=distance[strobes]
        mask=(sd>1e-7)&(sd<1e-3)&(times[strobes]<2000)
        slope=np.polyfit(times[strobes][mask],np.log(sd[mask]),1) if mask.sum()>5 else None
        thresholds={}
        for threshold in [1e-5,1e-4,1e-3,.01,.1]:
            hits=np.flatnonzero(sd>threshold)
            thresholds[str(threshold)]=float(times[strobes[hits[0]]]) if len(hits) else None
        matches=[]
        for index in strobes[::5]:
            current=q[index,case]
            def score(phase,details=False):
                target=spline(phase%period)
                cost=np.sum((current[:,None,:]-target[None,:,:])**2,axis=2)
                rows,cols=linear_sum_assignment(cost)
                return (float(np.sqrt(cost[rows,cols].sum())),cols) if details else float(cost[rows,cols].sum())
            costs=[]
            for target in coarse_q:
                cost=np.sum((current[:,None,:]-target[None,:,:])**2,axis=2)
                rows,cols=linear_sum_assignment(cost);costs.append(cost[rows,cols].sum())
            guess=coarse[int(np.argmin(costs))]
            fitted=minimize_scalar(score,bounds=(guess-period/80,guess+period/80),method='bounded',options={'xatol':1e-12})
            dist,perm=score(fitted.x,True)
            pairs=[]
            for targetpair in [[2,3],[6,7]]:
                pairs.append([int(j) for j in np.flatnonzero(np.isin(perm,targetpair))])
            # Which current pair maps to the reference's expanding pair?
            matches.append(dict(tau=float(times[index]),distance_to_permuted_saddle=dist,
                phase=float(fitted.x%period),permutation_current_from_reference=perm.tolist(),
                contracting_pair=pairs[0],expanding_pair=pairs[1]))
        pairs=[]
        for i in range(8):
            for j in range(i+1,8):
                last=np.linalg.norm(q[-4000:,case,i]-q[-4000:,case,j],axis=1)
                pairs.append(dict(rows=[i,j],last100periods_max_distance=float(max(last)),last100periods_rms=float(np.sqrt(np.mean(last**2)))))
        records.append(dict(case=case,initial=meta['cases'][case],linear_growth_fit=None if slope is None else float(slope[0]),
            linear_growth_fit_samples=int(mask.sum()),threshold_first_stroboscopic_tau=thresholds,
            orbit_matches=matches,last_pair_statistics=pairs))
        print('analyzed escape',case,'growth',None if slope is None else slope[0],'last match',matches[-1],flush=True)
    write(output,dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases=records,
          expected_local_growth=meta['unstable_exponent'],phase_fit='Periodic cubic interpolation of 1,001 archived orbit samples; all 8! row permutations minimized by Hungarian assignment at each fitted phase.',
          caveat='Matching and escape are floating-point evidence, not interval bounds or a proof of an infinite-time connection.',
          data_sha256=hashlib.sha256((folder/'trajectory.npz').read_bytes()).hexdigest(),
          code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))


def controls(source,folder,half_folder=None):
    output=folder/'growth-controls.json'
    if output.exists():raise FileExistsError(output)
    a=np.load(folder/'trajectory.npz');meta=json.loads((folder/'summary.json').read_text())
    diagnostic=json.loads((folder/'diagnostics.json').read_text());ref=np.load(source/'cycle.npz')
    base=ref['full_states'][0];p=ref['parameters'];q=physical(a['states']);times=a['time'];rows=[]
    for case in range(a['states'].shape[1]):
        strobes=np.arange(0,len(times),40)
        amplitude=np.linalg.norm(q[strobes,case,6]-q[strobes,case,7],axis=1)
        mask=(amplitude>1e-7)&(amplitude<1e-3)&(times[strobes]<2000)
        fit=np.polyfit(times[strobes][mask],np.log(amplitude[mask]),1)
        target=q[-1,case];match=diagnostic['cases'][case]['orbit_matches'][-1]
        phase=match['phase'];permutation=np.array(match['permutation_current_from_reference'])
        for _ in range(3):
            sol=solve_ivp(lambda t,z:full_rhs(t,z,p,8),(0,phase),base,method='DOP853',
                          rtol=1e-12,atol=1e-14,max_step=.0125)
            if not sol.success:raise RuntimeError(sol.message)
            mapped=physical(sol.y[:,-1])[permutation]
            direction=physical(full_rhs(0,sol.y[:,-1],p,8))[permutation]
            phase+=float(np.sum((target-mapped)*direction)/np.sum(direction**2))
        rows.append(dict(case=case,fitted_growth_rate=float(fit[0]),fitted_log_intercept=float(fit[1]),
            extrapolated_tau_to_pair_distance_1e_minus_4=float((np.log(1e-4)-fit[1])/fit[0]),
            endpoint_phase=phase,endpoint_permutation=permutation.tolist(),
            direct_ODE_endpoint_match_max_abs=float(np.max(abs(target-mapped))),
            direct_ODE_endpoint_match_l2=float(np.linalg.norm(target-mapped))))
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases=rows,
        floquet_growth=meta['unstable_exponent'],expected_epsilon_100_delay=float(np.log(100)/meta['unstable_exponent']),
        measured_epsilon_100_delays=[rows[j+2]['extrapolated_tau_to_pair_distance_1e_minus_4']-rows[j]['extrapolated_tau_to_pair_distance_1e_minus_4'] for j in [0,1]],
        invariant_pair='Rows2,3 zero-based are retained exactly by the seven-cluster coordinates',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        caveat='Numerical finite-time connection and local-growth controls, not a rigorous heteroclinic theorem.')
    if half_folder:
        b=np.load(half_folder/'trajectory.npz')
        assert np.array_equal(a['parameters'],b['parameters']) and np.array_equal(a['initial_full_state'],b['initial_full_state'])
        assert np.max(abs(a['time']-b['time']))<1e-9
        record['half_step']=dict(folder=half_folder.relative_to(ROOT).as_posix(),
            initial_states_exactly_equal=True,maximum_full_state_difference=float(np.max(abs(a['states']-b['states']))),
            final_full_state_difference=float(np.max(abs(a['states'][-1]-b['states'][-1]))),
            data_sha256=hashlib.sha256((half_folder/'trajectory.npz').read_bytes()).hexdigest())
    write(output,record);(folder/'controls.source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['run','analyze','controls'])
    parser.add_argument('--source',type=Path,default=ROOT/'results/motor-postcritical/six-cluster-cycle-rho-3p268')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--periods',type=int,default=400)
    parser.add_argument('--max-step',type=float,default=.05)
    parser.add_argument('--half-folder',type=Path)
    a=parser.parse_args()
    if a.command=='run':run(a.source,a.output,a.periods,a.max_step)
    elif a.command=='analyze':analyze(a.source.resolve(),a.output.resolve())
    else:controls(a.source.resolve(),a.output.resolve(),a.half_folder.resolve() if a.half_folder else None)
