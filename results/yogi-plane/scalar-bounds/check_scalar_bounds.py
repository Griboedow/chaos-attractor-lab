"""Exact algebra controls and finite stress checks for the scalar first-exit bound.

The written proof establishes infinite-time boundedness. These controls only
check algebraic hypotheses and implementations; their bounds are intentionally
loose and passing them is not a numerical proof of the theorem.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numba import njit, prange, set_num_threads


def log_bound(b, e, sigma, x, u):
    kappa = min(e, 1.) / 8
    length = max(1, math.ceil(math.log(kappa) / math.log(b))) if b else 1
    while b**length > kappa:
        length += 1
    logk = -math.log(sigma * kappa)
    radius = max(abs(x), abs(u))
    initial = math.log(2 * radius) + length * math.log1p(1/e) if radius else -math.inf
    return max(initial, math.log(8*length/e) + logk, math.log(4*(length+1)) + logk), length


@njit(parallel=True)
def run_stress(parameters, starts, steps):
    maxima = np.zeros((len(parameters), 3))
    tails = np.zeros((len(parameters), 512, 3))
    floor_errors = np.zeros(len(parameters))
    for i in prange(len(parameters)):
        b, c, e, kind, corrected = parameters[i]
        x, m, v = starts[i]
        maxima[i] = np.abs(starts[i])
        for t in range(1, steps+1):
            mn = b*m + (1-b)*x
            if kind == 0:  # Yogi, evaluated independently of production code.
                sign = 1 if v > x*x else (-1 if v < x*x else 0)
                vn = v - (1-c)*sign*x*x
                floor = min(c, 1-c)*x*x
                floor_errors[i] = max(floor_errors[i], (floor-vn)/(1+floor))
                den = math.sqrt(vn/(1-c**t) if corrected else vn) + e
            elif kind == 1:  # Adam.
                vn = c*v + (1-c)*x*x
                floor = (1-c)*x*x
                floor_errors[i] = max(floor_errors[i], (floor-vn)/(1+floor))
                den = math.sqrt(vn/(1-c**t) if corrected else vn) + e
            else:  # Adamax convention with epsilon inside the maximum.
                vn = max(c*v, abs(x)+e)
                floor_errors[i] = max(floor_errors[i], (abs(x)+e-vn)/(1+abs(x)+e))
                den = vn
            x = x - (mn/(1-b**t) if corrected else mn)/den
            m, v = mn, vn
            maxima[i,0] = max(maxima[i,0], abs(x))
            maxima[i,1] = max(maxima[i,1], abs(m))
            maxima[i,2] = max(maxima[i,2], abs(v))
            if t > steps-512:
                tails[i,t-(steps-512)-1] = (x,m,v)
    return maxima, tails, floor_errors


def exact_controls():
    rng = np.random.default_rng(2026110101)
    branch_counts = {-1:0, 0:0, 1:0}
    for i in range(1024):
        c = Q(int(rng.integers(1,1000)), 1000)
        b = Q(int(rng.integers(0,1000)), 1000)
        e = Q(int(rng.integers(1,1000)), int(rng.integers(1,1000)))
        x = Q(int(rng.integers(-1000,1001)), int(rng.integers(1,1000)))
        v = x*x if i % 3 == 0 else Q(int(rng.integers(0,1000)),100)
        sign = (v > x*x) - (v < x*x)
        vn = v - (1-c)*sign*x*x
        assert vn >= min(c,1-c)*x*x and vn >= c*v
        branch_counts[sign] += 1
        kappa = min(e, Q(1))/8
        assert kappa*(1+2/e) <= Q(3,8)
        assert Q(3,4)-Q(7,4)*kappa > 0
        t = int(rng.integers(1,32))
        a = b*(1-b**t)/(1-b**(t+1))
        assert 0 <= a <= b
        assert 1-a == (1-b)/(1-b**(t+1))
        # Nonzero initial moments are permitted by the corrected recurrence.
        m = Q(int(rng.integers(-1000,1001)), 73)
        assert (b*m+(1-b)*x)/(1-b**(t+1)) == a*m/(1-b**t)+(1-a)*x
    return dict(cases=1024, yogi_branches=branch_counts,
                arithmetic='fractions.Fraction; all asserted identities exact')


def main(args):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    set_num_threads(4)
    exact = exact_controls()
    rng = np.random.default_rng(2026110102)
    base = np.column_stack([1-10**rng.uniform(-3,0,64),
                            1-10**rng.uniform(-3,0,64),
                            10**rng.uniform(-8,2,64)])
    base[:4] = [[0,.5,.1],[.99,.01,1e-7],[.01,.99,10],[.999,.999,1e-8]]
    starts0 = np.column_stack([rng.choice([-1.,1.],64)*10**rng.uniform(-6,6,64),
                               rng.choice([-1.,1.],64)*10**rng.uniform(-6,6,64),
                               10**rng.uniform(-10,8,64)])
    parameters = np.array([list(p)+[kind,corrected] for kind in range(3)
                           for corrected in range(2) for p in base])
    starts = np.tile(starts0,(6,1))
    maxima,tails,errors = run_stress(parameters,starts,100000)
    assert np.isfinite(maxima).all() and np.isfinite(tails).all()
    assert errors.max() < 1e-12
    # The explicit B is for the autonomous start; corrected orbits start after
    # their first update, and are checked below with that corrected first moment.
    bounds=[]
    for p,s,maximum in zip(parameters,starts,maxima):
        b,c,e,kind,corrected=p; x,m,v=s
        sigma = math.sqrt(min(c,1-c)) if kind==0 else (math.sqrt(1-c) if kind==1 else 1.)
        if corrected:
            mn=b*m+(1-b)*x
            vn=v-(1-c)*np.sign(v-x*x)*x*x if kind==0 else (c*v+(1-c)*x*x if kind==1 else max(c*v,abs(x)+e))
            den=math.sqrt(vn/(1-c))+e if kind<2 else vn
            m=mn/(1-b); x=x-m/den
        logb,length=log_bound(b,e,sigma,x,m)
        # maxima includes the initial raw state; include it explicitly.
        raw_log=max(math.log(max(abs(s[0]),abs(s[1]),1e-300)),logb)
        assert math.log(max(maximum[0],maximum[1],1e-300)) <= raw_log+1e-12
        bounds.append([logb,length])
    np.savez_compressed(out/'histories.npz',parameters=parameters,starts=starts,
                        maxima=maxima,tails=tails,denominator_floor_errors=errors,
                        log_bounds_and_memory_lengths=bounds)
    result=dict(exact_controls=exact,finite_stress=dict(cases=len(parameters),steps=100000,
                kinds=['Yogi','Adam','Adamax'],corrected=[False,True],
                nonzero_initial_memories=True,max_floor_residual=float(errors.max())),
                scope='Implementation controls, not a substitute for the first-exit proof. Large explicit bounds are stored as natural logarithms to avoid overflow.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (out/'proof.md').write_bytes(Path(__file__).parents[1].joinpath('docs/scalar-optimizer-boundedness.md').read_bytes())
    (out/'hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest()
        for p in out.iterdir() if p.is_file()},indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
