"""Separate within-cluster chaos estimates from conditional transverse growth."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
import motor_rows
import motor_postcritical_cycles as cycles


def groups_from_state(state):
    rows=np.stack([state[:8],state[8:16],np.r_[state[16:],-sum(state[16:])]],axis=1)
    todo=set(range(8));groups=[]
    while todo:
        first=min(todo);group=[k for k in sorted(todo) if np.linalg.norm(rows[k]-rows[first])<1e-10]
        groups.append(group);todo-=set(group)
    return groups


def contrast_matrix(state,p,w,j):
    c,s,x=cycles.unpack(state,w);v=-p[0]*x+p[1]*(s-np.dot(w,s));rho,r,_=p
    return np.array([[-1,-v[j]-r*s[j],rho*s[j]],[v[j],-1+r*c[j],-rho*c[j]],[0,r,-rho]])


def run_lane(lane,source,out,wall,observed_tau=2000.,max_step=.03,rtol=1e-11,atol=1e-13):
    old=np.load(source);i=int(np.flatnonzero(old['source_lane_indices']==lane)[0])
    full=old['final_states'][i];p=old['parameters'][i];groups=groups_from_state(full)
    w=np.array(list(map(len,groups)))/8.;x=cycles.restrict(full,groups);d=len(x)
    paired=[j for j,g in enumerate(groups) if len(g)>1]
    if not paired:raise ValueError('A clustered endpoint is required')
    out.mkdir(parents=True)
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_lane=lane,
        source_path=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        groups=groups,weights=w.tolist(),parameters=p.tolist(),quotient_dimension=d,
        initial_projection_error=float(np.max(abs(cycles.lift(x,groups)-full))),
        restriction_controls=cycles.audit(groups,p),windows=[],status='running',
        settings=dict(warm_tau=200.,observed_tau_requested=observed_tau,qr_interval=1.,max_step=max_step,rtol=rtol,atol=atol),
        interpretation='A declared invariant cluster restriction, initialized by a measured tiny projection of an archived endpoint. Quotient and omitted contrast growth are computed separately.',
        limitations=['Finite-time conditional exponents, not transverse stability proof or asymptotic chaos.',
            'Positive transverse growth would rule out interpreting this restricted numerical trajectory alone as evidence of full-space attraction.',
            'No change to the original cohort trajectories or counts.'])
    max_error=0.
    for trial in (x,x+.013*np.sin(np.arange(d))):
        lifted=cycles.lift(trial,groups);_,jf=motor_rows.field_jac(lifted,p,8)
        for j in paired:
            a,b=groups[j][0],groups[j][-1];embedding=np.zeros((23,3))
            embedding[a,0]=1;embedding[b,0]=-1;embedding[8+a,1]=1;embedding[8+b,1]=-1
            if a<7:embedding[16+a,2]=1
            if b<7:embedding[16+b,2]=-1
            max_error=max(max_error,float(np.max(abs(jf@embedding-embedding@contrast_matrix(trial,p,w,j)))))
    record['transverse_full_jacobian_error']=max_error
    if max_error>1e-10:raise ArithmeticError('Contrast decomposition control failed')
    sources=[Path(__file__),Path(cycles.__file__),Path(motor_rows.__file__)]
    record['source_snapshot_sha256']={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in sources}
    for q in sources:(out/q.name).write_bytes(q.read_bytes())
    v=np.sin(np.arange(d)+1.);v/=np.linalg.norm(v);bases=np.repeat(np.eye(3)[None],len(paired),axis=0)
    logs=np.zeros(1+3*len(paired));observed=0.;clock=0.;times=[];states=[];begin=time.perf_counter()
    def save():
        record.update(elapsed_tau=clock,observed_tau=observed,wall_seconds=time.perf_counter()-begin)
        (out/'analysis.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        np.savez_compressed(out/'history.npz',time=np.array(times),states=np.array(states),parameters=p,
            final_state=x,final_quotient_tangent=v,final_contrast_bases=bases,log_sums=logs)
    def rhs(t,z):
        y=z[:d];f,j=cycles.field_jac(y,p,w)
        parts=[f,j@z[d:2*d]]
        for k,a in enumerate(paired):parts.append((contrast_matrix(y,p,w,a)@z[2*d+9*k:2*d+9*k+9].reshape(3,3)).ravel())
        return np.concatenate(parts)
    for block,length in enumerate([200.]+[500.]*(int(observed_tau)//500)):
        block_logs=np.zeros_like(logs)
        for step in range(int(length)):
            if time.perf_counter()-begin>wall:
                record['status']='wall budget incomplete';save();return
            z=np.r_[x,v,bases.ravel()]
            sol=solve_ivp(rhs,(0,1.),z,method='DOP853',rtol=rtol,atol=atol,max_step=max_step)
            if not sol.success:raise RuntimeError(sol.message)
            end=sol.y[:,-1];x=end[:d];v=end[d:2*d];norm=np.linalg.norm(v);v/=norm
            block_logs[0]+=np.log(norm)
            for k,a in enumerate(paired):
                mat=end[2*d+9*k:2*d+9*k+9].reshape(3,3);q,r=np.linalg.qr(mat)
                block_logs[1+3*k:4+3*k]+=np.log(abs(np.diag(r)));bases[k]=q
            clock+=1.;times.append(clock);states.append(x.copy())
        if block:
            observed+=length;logs+=block_logs
            record['windows'].append(dict(duration=length,end_tau=clock,quotient_leading_ftle=float(block_logs[0]/length),
                contrast_ftle=(block_logs[1:].reshape(-1,3)/length).tolist()))
            record['quotient_leading_ftle']=float(logs[0]/observed)
            record['contrast_ftle']=(logs[1:].reshape(-1,3)/observed).tolist()
        save();print('cluster',lane,'tau',clock,'quotient',record.get('quotient_leading_ftle'),'contrasts',record.get('contrast_ftle'),flush=True)
    record['status']='complete';save()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--lanes',type=int,nargs='+',default=[36,42,46])
    parser.add_argument('--wall-per-lane',type=float,default=1200.)
    parser.add_argument('--observed',type=float,default=2000.)
    parser.add_argument('--max-step',type=float,default=.03)
    parser.add_argument('--rtol',type=float,default=1e-11)
    parser.add_argument('--atol',type=float,default=1e-13)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output')
    if args.observed<=0 or args.observed%500 or min(args.max_step,args.rtol,args.atol)<=0:raise ValueError('Positive settings and a positive multiple of500 observed time required')
    source=Path(__file__).resolve().parents[1]/'results/motor-rows/activity-path/N8-lifetime-history.npz'
    for lane in args.lanes:run_lane(lane,source,args.output/f'lane-{lane}',args.wall_per_lane,args.observed,args.max_step,args.rtol,args.atol)


if __name__=='__main__':main()
