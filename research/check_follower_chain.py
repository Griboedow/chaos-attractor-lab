"""Independent adaptive full-spectrum and contact audit for a saved chain scout."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from search_follower_chain import field


def coefficients(p):
    n = (len(p)-1)//4+1
    m = n-1
    stiffness = np.r_[1., p[1:1+m]]
    drag = np.r_[1., p[1+m:1+2*m]]
    length = np.r_[1., p[1+2*m:1+3*m]]
    rest = np.r_[0., np.cumsum(p[1+3*m:])]
    difference = np.eye(n)-np.eye(n, k=-1)
    elastic = difference.T@np.diag(stiffness)@difference
    tail = np.cumsum(drag[::-1])[::-1]
    base = tail[np.maximum.outer(np.arange(n), np.arange(n))]*np.outer(length, length)
    return p[0], elastic, base, length, rest


@njit(cache=True)
def field_jac(x, force, elastic, base, length, rest):
    n = len(x)
    resistance = np.empty((n, n))
    sine = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            resistance[i,j] = base[i,j]*np.cos(x[i]-x[j])
            sine[i,j] = base[i,j]*np.sin(x[i]-x[j])
    f = np.linalg.solve(resistance, force*length*np.sin(x-x[-1])-elastic@(x-rest))
    rhs = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            rhs[i,j] = force*length[i]*np.cos(x[i]-x[-1])*((i==j)-(n-1==j))-elastic[i,j]
            for k in range(n):
                rhs[i,j] += sine[i,k]*((i==j)-(k==j))*f[k]
    return f, np.linalg.solve(resistance, rhs)


@njit(cache=True)
def augmented(z, force, elastic, base, length, rest):
    n = len(length)
    f, j = field_jac(z[:n], force, elastic, base, length, rest)
    result = np.empty(n+n*n+1)
    result[:n] = f
    result[n:-1] = (j@z[n:-1].reshape(n,n)).ravel()
    result[-1] = np.trace(j)
    return result


def contact_audit(states, length):
    directions = np.stack((np.cos(states), np.sin(states)), axis=2)*length[None,:,None]
    positions = np.concatenate((np.zeros((len(states),1,2)), np.cumsum(directions,axis=1)),axis=1)
    def cross(a,b):
        return a[:,0]*b[:,1]-a[:,1]*b[:,0]
    any_crossing = np.zeros(len(states),bool)
    pairs = []
    for i in range(len(length)):
        for j in range(i+2,len(length)):
            a,b,c,d = positions[:,i],positions[:,i+1],positions[:,j],positions[:,j+1]
            intersect = (cross(b-a,c-a)*cross(b-a,d-a)<0)&(cross(d-c,a-c)*cross(d-c,b-c)<0)
            any_crossing |= intersect
            pairs.append(dict(segments=[i+1,j+1],sampled_crossings=int(intersect.sum())))
    return dict(samples=len(states),sampled_strict_intersections=int(any_crossing.sum()),
                sampled_intersection_fraction=float(any_crossing.mean()),pairs=pairs,
                warning='Zero-thickness nonadjacent segment crossings at saved samples only; no continuous clearance certificate.')


def write(path, record):
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    source = np.load(args.source)
    p = source['parameters'][args.id]
    n = (len(p)-1)//4+1
    x = source['states'][args.id,:n].copy()
    initial = x.copy()
    coeff = coefficients(p)
    rng = np.random.default_rng(20261031824)
    error = 0.
    for y in [x,*rng.uniform(-4,4,(16,n))]:
        f,j = field_jac(y,*coeff)
        reference = np.column_stack([field(y+1e-25j*np.eye(n)[i],p).imag/1e-25 for i in range(n)])
        error = max(error,float(np.max(abs(j-reference))/(1+np.max(abs(reference)))),
                    float(np.max(abs(f-field(y,p)))/(1+np.max(abs(f)))))
    assert error < 1e-10
    short = []
    for method in ['DOP853','Radau']:
        sol = solve_ivp(lambda t,y:field(y,p),(0,2),x,method=method,rtol=2e-11,atol=2e-13,max_step=.005)
        assert sol.success
        short.append(sol.y[:,-1])
    short_error = float(np.max(abs(short[0]-short[1])))
    assert short_error < 1e-7
    write(out/'protocol.json',dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        id=args.id,links=n,parameters=p.tolist(),initial=initial.tolist(),warm=args.warm,duration=args.duration,
        method='DOP853 full variational QR; independent general-matrix complex-step and Radau controls',
        rtol=1e-10,atol=1e-12,max_step=args.step,qr=.1,sampling=.01,
        field_jac_relative_error=error,short_2_tau_DOP853_Radau_error=short_error,
        selection='Retained GPU scout; no new-attractor claim.'))
    for name in ['check_follower_chain.py','search_follower_chain.py','follower_link.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    basis = np.eye(n)
    logs = np.zeros(n)
    trace = 0.
    times,states,windows = [],[],[]
    previous = np.zeros(n)
    prevtrace = 0.
    section_times,section_states = [],[]
    began = time.perf_counter()
    def section(t,z):
        return z[0]
    section.direction = 1
    total = int(round((args.warm+args.duration)/.1))
    warm_steps = int(round(args.warm/.1))
    window_steps = int(round(args.duration/.4))
    assert total == warm_steps+4*window_steps
    for iteration in range(total):
        start,end = iteration*.1,(iteration+1)*.1
        grid = np.linspace(start,end,11)
        sol = solve_ivp(lambda t,z:augmented(z,*coeff),(start,end),np.r_[x,basis.ravel(),0.],
                        method='DOP853',rtol=1e-10,atol=1e-12,max_step=args.step,t_eval=grid,events=section)
        if not sol.success:
            raise RuntimeError(sol.message)
        x = sol.y[:n,-1]
        basis,r = np.linalg.qr(sol.y[n:-1,-1].reshape(n,n))
        if iteration >= warm_steps:
            logs += np.log(abs(np.diag(r)))
            trace += sol.y[-1,-1]
            times.extend(grid[1:])
            states.extend(sol.y[:n,1:].T)
            section_times.extend(sol.t_events[0])
            section_states.extend(np.asarray(sol.y_events[0]).reshape(-1,n+n*n+1)[:,:n])
        if iteration+1 > warm_steps and (iteration+1-warm_steps)%window_steps == 0:
            window = dict(end=end,spectrum=((logs-previous)/(args.duration/4)).tolist(),
                          mean_divergence=(trace-prevtrace)/(args.duration/4))
            windows.append(window)
            previous,prevtrace = logs.copy(),trace
            np.savez_compressed(out/f'checkpoint-{iteration+1:06}.npz',time=times,states=states,parameters=p,
                                final_state=x,final_basis=basis,log_sums=logs,integrated_trace=trace,
                                section_times=section_times,section_states=np.array(section_states).reshape(-1,n))
            write(out/'progress.json',dict(windows=windows,elapsed=end,wall_seconds=time.perf_counter()-began))
            print(args.id,window,flush=True)
    states = np.array(states)
    np.savez_compressed(out/'trace.npz',time=times,states=states,parameters=p,final_state=x,final_basis=basis,
                        log_sums=logs,integrated_trace=trace,section_times=section_times,
                        section_states=np.array(section_states).reshape(-1,n))
    result = dict(spectrum=(logs/args.duration).tolist(),windows=windows,
                  trace_residual_per_tau=float(abs(logs.sum()-trace)/args.duration),
                  min_state=states.min(axis=0).tolist(),max_state=states.max(axis=0).tolist(),sections=len(section_times),
                  crossing_check=contact_audit(states,coeff[3]),wall_seconds=time.perf_counter()-began,
                  status='Finite-time numerical evidence only; no historical or topological novelty claim.')
    write(out/'analysis.json',result)
    assert result['trace_residual_per_tau'] < 1e-6


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--id',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--warm',type=float,default=100.)
    p.add_argument('--duration',type=float,default=2000.)
    p.add_argument('--step',type=float,default=.01)
    main(p.parse_args())
