"""Replay the stated AdaDelta cycle convention, retaining epsilon exactly."""
import argparse,hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads


@njit(parallel=True)
def replay(parameters,starts,times):
    snapshots=np.empty((len(parameters),len(times),3));tails=np.empty((len(parameters),2048,3))
    for i in prange(len(parameters)):
        alpha,beta,curvature,epsilon2=parameters[i];x,d,v=starts[i];cursor=0
        for t in range(times[-1]+1):
            if t==times[cursor]:snapshots[i,cursor]=(x,d,v);cursor=min(cursor+1,len(times)-1)
            if t>times[-1]-2048:tails[i,t-(times[-1]-2048)-1]=(x,d,v)
            if t<times[-1]:
                g=curvature*x;vn=beta*v+(1-beta)*g*g
                delta=g*np.sqrt((d+epsilon2)/(vn+epsilon2))
                dn=beta*d+(1-beta)*delta*delta
                x=x-alpha*delta;d=dn;v=vn
    return snapshots,tails


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    parameters=np.array([[.1,.9,1.,1e-8],[.1,.9,1.,0.],[2.,.9,1.,1e-8],[3.,.9,1.,1e-8],[.1,.9,1.,1e-6]])
    starts=np.array([[1.,400.,1.],[1.,400.,1.],[1.,1.,1.],[1.,4/9,1.],[1.,400.,1.]])
    times=np.unique(np.r_[0,1,2,5000,np.geomspace(1e4,1e8,161).astype(np.int64)])
    protocol=dict(parameter_columns=['alpha','beta','curvature','epsilon_squared'],parameters=parameters.tolist(),starts_x_delta_memory_gradient_memory=starts.tolist(),
        times=times.tolist(),update='V=beta*v+(1-beta)*(curvature*x)^2; delta=curvature*x*sqrt((d+epsilon_squared)/(V+epsilon_squared)); D=beta*d+(1-beta)*delta^2; X=x-alpha*delta.',
        source='Bock dissertation DOI10.5283/epub.58326, equation5.2, section5.2.3, theorem5.4 and figure52. Typeset physical pages96,106,107 inspected.',
        design='First row reproduces stated figure52 parameters and initialization. Remaining rows separately declare zero-epsilon, critical-product, above-critical and larger-epsilon controls.',
        scope='Finite replay and exact algebra control. No global convergence or divergence theorem is inferred from a finite trace; absence of noncritical period2 follows from the separate written proof.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    rng=np.random.default_rng(2026110152);controls=0
    for _ in range(512):
        beta=Q(int(rng.integers(0,1000)),1000);curv=Q(int(rng.integers(1,1000)),71)
        a=Q(int(rng.integers(1,1000)),137);b=Q(int(rng.integers(1,1000)),113);E=Q(int(rng.integers(1,1000)),10**8)
        va=curv**2*(a*a+beta*b*b)/(1+beta);vb=curv**2*(beta*a*a+b*b)/(1+beta)
        assert (va+E)*b*b-(vb+E)*a*a == (b*b-a*a)*(curv**2*beta*(a*a+b*b)/(1+beta)+E)
        assert curv**2*beta*(a*a+b*b)/(1+beta)+E>0
        controls+=1
    exact_delta_memory_change=(Q(1,10)*(Q(400)+Q(1,10**8))/(Q(1)+Q(1,10**8))+Q(9,10)*400)-400
    assert exact_delta_memory_change<0
    start=time.perf_counter();snapshots,tails=replay(parameters,starts,times)
    assert np.isfinite(snapshots).all() and np.isfinite(tails).all()
    np.savez_compressed(out/'histories.npz',parameters=parameters,starts=starts,times=times,snapshots=snapshots,tails=tails)
    rows=[]
    for i in range(len(parameters)):
        row=dict(index=i,parameters=parameters[i].tolist(),state_after_two=snapshots[i,list(times).index(2)].tolist(),
            two_step_initial_residual=(snapshots[i,list(times).index(2)]-starts[i]).tolist(),state_after5000=snapshots[i,list(times).index(5000)].tolist(),
            final_state=snapshots[i,-1].tolist(),final_tail_x_span=float(np.ptp(tails[i,:,0])),final_tail_max_abs_x=float(np.max(abs(tails[i,:,0]))))
        rows.append(row);print(row,flush=True)
    result=dict(runs=rows,completed=len(rows),declared=len(parameters),steps_per_run=int(times[-1]),wall_seconds=time.perf_counter()-start,
        exact_algebra_cases=controls,source_start_one_step_delta_memory_change=str(exact_delta_memory_change),
        identity='(Va+E)b^2-(Vb+E)a^2=(b^2-a^2)[curvature^2*beta*(a^2+b^2)/(1+beta)+E]',
        scope=protocol['scope'])
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
