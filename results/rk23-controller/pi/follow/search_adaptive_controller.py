"""Screen explicitly defined Heun/Euler PI-controller maps on y'=-y.

This is a transcription experiment, not an execution of SUNDIALS. The raw
proposal rule is a subset of its documented Soderlind controller. Policy 1
adds our declared error rejection and factor limiter, not ARKODE's policy.
States are logarithms of (proposed h, |y|/atol, previous accepted error).
All outcomes, including escapes from the declared numerical box, are saved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from numba import njit, prange, set_num_threads
from scipy.stats import qmc


@njit(cache=True)
def advance(z, p, policy):
    a, b, safety = p
    x, y, previous = z
    dx = np.array([1.0, 0.0, 0.0])
    dy = np.array([0.0, 1.0, 0.0])
    dz = np.array([0.0, 0.0, 1.0])
    rejected = 0
    if policy == 1 and 2*x+y-np.log(2.0) > 0.0:
        # Recompute the same step with an I-controller shrink. For this
        # particular embedded pair this makes its error exactly safety**2.
        x = 0.5*(np.log(2.0)-y)+np.log(safety)
        dx = -0.5*dy
        rejected = 1
    h = np.exp(x)
    stability = 1.0-h+0.5*h*h
    error = 2*x+y-np.log(2.0)
    de = 2*dx+dy
    logfactor = np.log(safety)-a*error+b*previous
    dfactor = -a*de+b*dz
    limited = 0
    if policy == 1:
        if logfactor < np.log(0.1):
            logfactor = np.log(0.1)
            dfactor *= 0.0
            limited = 1
        elif logfactor > np.log(5.0):
            logfactor = np.log(5.0)
            dfactor *= 0.0
            limited = 1
    out = np.array([x+logfactor, y+np.log(stability), error])
    jac = np.empty((3,3))
    jac[0] = dx+dfactor
    jac[1] = dy+h*(h-1.0)/stability*dx
    jac[2] = de
    return out, jac, rejected, limited


def start(safety):
    # h0=.1, u0=1, and an I-controller first step before PI history exists.
    h0 = 0.1
    error = h0*h0/2
    return np.log([h0*safety/np.sqrt(error), 1-h0+h0*h0/2, error])


@njit(cache=True)
def run_one(p, policy, initial, warm, steps, windows):
    z = initial.copy()
    tangent = np.array([1.0, 0.371, -0.217])
    tangent /= np.sqrt(np.sum(tangent*tangent))
    rates = np.zeros(windows)
    low = np.full(3, np.inf)
    high = np.full(3, -np.inf)
    rejects = 0
    limits = 0
    status = 0
    completed = 0
    for i in range(warm+steps):
        if not np.all(np.isfinite(z)) or abs(z[0])>50 or abs(z[1])>100 or abs(z[2])>200:
            status = 1
            break
        z, jac, rej, lim = advance(z, p, policy)
        tangent = jac@tangent
        norm = np.sqrt(np.sum(tangent*tangent))
        if not np.isfinite(norm) or norm < 1e-290:
            status = 2
            break
        tangent /= norm
        if i>=warm:
            completed += 1
            rates[(i-warm)//(steps//windows)] += np.log(norm)
            low = np.minimum(low,z)
            high = np.maximum(high,z)
            rejects += rej
            limits += lim
    rates /= steps//windows
    return rates, low, high, z, status, completed, rejects, limits


@njit(cache=True, parallel=True)
def batch(params, policies, initials, warm, steps, windows):
    n = len(params)
    rates = np.empty((n,windows))
    bounds = np.empty((n,2,3))
    final = np.empty((n,3))
    counters = np.empty((n,4),dtype=np.int64)
    for i in prange(n):
        r, lo, hi, z, status, count, reject, limit = run_one(params[i],int(policies[i]),initials[i],warm,steps,windows)
        rates[i] = r
        bounds[i,0] = lo
        bounds[i,1] = hi
        final[i] = z
        counters[i] = np.array([status,count,reject,limit])
    return rates, bounds, final, counters


def verify():
    rng = np.random.default_rng(2026110171)
    errors=[]
    for policy in (0,1):
        for _ in range(256):
            z = rng.uniform(-2,2,3)
            p = np.array([rng.uniform(.2,2),rng.uniform(.01,.8),rng.uniform(.5,.99)])
            _, jac, _, _ = advance(z,p,policy)
            finite = np.empty((3,3))
            for j in range(3):
                d = np.eye(3)[j]*1e-6
                finite[:,j] = (advance(z+d,p,policy)[0]-advance(z-d,p,policy)[0])/2e-6
            errors.append(np.max(np.abs(jac-finite)))
    maximum=float(max(errors))
    assert maximum < 1e-7,maximum
    # The documented PI coefficients k=(.8,-.31), embedded order p=1.
    p = np.array([.4,.155,.9])
    e = p[2]**(1/(p[0]-p[1]))
    fixed = np.log([2,e/2,e])
    value,jac,_,_=advance(fixed,p,0)
    assert np.max(np.abs(value-fixed))<1e-14
    assert abs(np.linalg.det(jac))<1e-14
    return {"finite_difference_checks":512,"max_absolute_jacobian_error":maximum,
            "default_fixed_eigenvalues":[[float(z.real),float(z.imag)] for z in np.linalg.eigvals(jac)]}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--power",type=int,default=12)
    ap.add_argument("--warm",type=int,default=10000)
    ap.add_argument("--steps",type=int,default=30000)
    ap.add_argument("--threads",type=int,default=12)
    args=ap.parse_args()
    assert args.steps%10==0
    args.out.mkdir(parents=True,exist_ok=False)
    set_num_threads(args.threads)
    check=verify()
    sobol=qmc.Sobol(3,scramble=True,seed=2026110170).random_base2(args.power)
    b=.001+1.499*sobol[:,1]
    params=np.column_stack([b+.025+1.975*sobol[:,0],b,.5+.49*sobol[:,2]])
    controls=np.array([[.4,.155,.9],[.5,0,.9],[.3175,.134,.9]])
    params=np.vstack([params,controls])
    policies=np.repeat(np.arange(2),len(params))
    params=np.tile(params,(2,1))
    initials=np.array([start(p[2]) for p in params])
    begin=time.time()
    rates,bounds,final,counters=batch(params,policies,initials,args.warm,args.steps,10)
    flags=(counters[:,0]==0)&(rates.min(axis=1)>.001)
    np.savez_compressed(args.out/"screen.npz",params=params,policies=policies,initials=initials,
                        rates=rates,bounds=bounds,final=final,counters=counters,flags=flags)
    report={"schema":1,"experiment":"Heun/Euler PI on scalar decay; raw and declared guarded maps, not an upstream solver run",
            "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "parameters":["a=current-error exponent","b=previous-error exponent","safety"],
            "state":["log(proposed h)","log(abs(y)/atol)","log(previous accepted normalized error)"],
            "policies":{"0":"raw PI, all steps accepted","1":"if error>1 recompute h=safety*sqrt(2/u), then PI factor clamped to [.1,5]"},
            "status_codes":{"0":"completed","1":"left declared log box or nonfinite","2":"tangent underflow/nonfinite"},
            "seed":2026110170,"sobol_power":args.power,"warm":args.warm,"steps":args.steps,"windows":10,
            "checks":check,"count":len(params),"flag_count":int(flags.sum()),
            "completed_by_policy":[int(np.sum((policies==j)&(counters[:,0]==0))) for j in range(2)],
            "flags_by_policy":[int(np.sum((policies==j)&flags)) for j in range(2)],
            "seconds":time.time()-begin}
    (args.out/"report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":
    main()
