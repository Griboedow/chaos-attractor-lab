"""Independent Cartesian equations, centered Jacobians and adaptive spectrum replay."""
import argparse
import json
import math
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from field import (make_case,make_ring,make_rotor,make_planes,velocity_tangent,ring_tangent,
                   rotor_tangent,plane_tangent,make_crossed,crossed_tangent)


@njit(cache=True)
def polar_ring(y,p):
    x, yy, z = y[:3]
    n = int(p[0])
    r = math.hypot(x, yy)
    theta = math.atan2(yy, x)
    if r == 0:
        return np.array([0., 0., -p[7]*z-p[8]])
    power = 2/(1+r*r)
    radial = r*(p[1]-r*r+p[2]*z) + p[6]*power**((n-2)/2)*r**(n-1)*math.cos(n*theta)
    angular = p[3]+p[4]*z+p[5]*r*r-p[6]*power**((n-2)/2)*r**(n-2)*math.sin(n*theta)
    axial = (-p[7]*z+p[8]*(r*r-1)/(1+r*r)
             + p[9]*power**(n/2)*r**n*math.cos(n*theta+p[10]))
    return np.array([radial*x/r-angular*yy, radial*yy/r+angular*x, axial])


@njit(cache=True)
def independent_rhs(y, p, family):
    x, yy, z = y[:3]
    if family == 3:
        rows = p[7:].reshape((int(p[0]),7))
        out = p[4:7]-p[1:4]*y[:3]
        for row in rows:
            signal = np.tanh(np.dot(row[:3],y[:3])+row[3])
            out += signal*row[4:]
        return out
    if family == 2:
        grad = p[6:9]+p[9:12]*y[:3]+p[12:15]*y[:3]**3+p[15]*np.array([yy*z,z*x,x*yy])
        return np.cross(y[:3],grad)/(1+np.dot(y[:3],y[:3])) + p[3:6] - p[:3]*y[:3]
    if family == 1:
        return polar_ring(y,p)
    if family == 4:
        rotation = p[22:31].reshape(3,3)
        shifted = rotation.T @ (y[:3]-p[31:34])
        return (1-p[34])*polar_ring(y,p[:11])+p[34]*(rotation @ polar_ring(shifted,p[11:22]))
    out = np.empty(3)
    for i in range(3):
        j, k = (i+1)%3, (i+2)%3
        d, a, b, c, phi, psi, chi, m, n = p[9*i:9*i+9]
        out[i] = (a*math.sin(y[j]+phi) + c*math.sin(m*y[j]+n*y[k]+chi)
                  + b*math.sin(y[k]+psi)) - d*y[i]
    return out


@njit(cache=True)
def full_rhs(t, state, p, family):
    position = state[:3]
    output = np.empty(13)
    output[:3] = independent_rhs(position, p, family)
    jac = np.empty((3, 3))
    for i in range(3):
        eps = 2e-6*(1+abs(position[i]))
        above = position.copy()
        below = position.copy()
        above[i] += eps
        below[i] -= eps
        jac[:,i] = (independent_rhs(above,p,family)-independent_rhs(below,p,family))/(2*eps)
    for i in range(3):
        for j in range(3):
            output[3+3*i+j] = jac[i,0]*state[3+j]+jac[i,1]*state[6+j]+jac[i,2]*state[9+j]
    output[12] = jac[0,0]+jac[1,1]+jac[2,2]
    return output


def replay(case, p, start, family, args):
    state = np.r_[start, np.eye(3).ravel(), 0.]
    direction = np.array([1., 2., -3.])
    analytic = (velocity_tangent,ring_tangent,rotor_tangent,plane_tangent,crossed_tangent)[family]
    probes = np.random.default_rng(189+case).uniform(-2,2,(12,3))
    checks = []
    for loc in probes:
        y = np.r_[loc, np.eye(3).ravel(),0.]
        reference = full_rhs(0.,y,p.ravel(),family)
        primary = analytic(np.r_[loc,direction],p)
        checks.append([max(abs(primary[:3]-reference[:3])),
                       max(abs(primary[3:]-reference[3:12].reshape(3,3)@direction))])
    rates = []
    divergence = []
    max_error = 0.
    nfev = 0
    for block in range(-1,args.blocks):
        elapsed = args.burn if block < 0 else args.window
        steps = round(elapsed / args.qr)
        logs = np.zeros(3)
        div = 0.
        for i in range(steps):
            result = solve_ivp(full_rhs,(0,args.qr),state,args=(p.ravel(),family),method='DOP853',
                               rtol=args.rtol,atol=args.rtol/100)
            if not result.success:
                raise RuntimeError(result.message)
            nfev += result.nfev
            state = result.y[:,-1].copy()
            q, r = np.linalg.qr(state[3:12].reshape(3,3))
            logs += np.log(np.abs(np.diag(r)))
            div += state[12]
            state[3:12] = q.ravel()
            state[12] = 0.
        if block >= 0:
            rates.append((logs/elapsed).tolist())
            divergence.append(div/elapsed)
            max_error = max(max_error, abs(logs.sum()-div)/elapsed)
    avg = np.mean(rates,axis=0)
    return dict(case=case, family=args.family, initial=start.tolist(),
                windows=rates, mean_spectrum=avg.tolist(), window_divergence=divergence,
                divergence_qr_maxerr=max_error, max_velocity_discrepancy=float(np.max(checks,axis=0)[0]),
                max_tangent_discrepancy=float(np.max(checks,axis=0)[1]), nfev=nfev)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--ids', required=True, nargs='+', type=int)
    ap.add_argument('--inputs',type=Path,help='Read already selected parameters and initial states')
    ap.add_argument('--seed', default=2026112123, type=int)
    ap.add_argument('--family', choices=['cartesian','ring','rotor','planes','crossed'], default='cartesian')
    ap.add_argument('--qr', type=float, default=2.)
    ap.add_argument('--window', type=float, default=160.)
    ap.add_argument('--burn', type=float, default=200.)
    ap.add_argument('--blocks', type=int, default=6)
    ap.add_argument('--rtol', type=float, default=1e-8)
    args = ap.parse_args()
    for val in [args.burn,args.window]:
        if abs(val/args.qr-round(val/args.qr))>1e-10:
            raise ValueError('Require complete QR intervals.')
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/'run.json').write_text(json.dumps({**vars(args),'out':str(args.out)},indent=2,default=str)+'\n')
    family = ('cartesian','ring','rotor','planes','crossed').index(args.family)
    generate = (make_case, make_ring, make_rotor, make_planes, make_crossed)[family]
    inputs = json.loads(args.inputs.read_text()) if args.inputs else {}
    for case in args.ids:
        parameters, start = generate(args.seed,case)
        if inputs:
            parameters=np.array(inputs[str(case)]['parameters']);start=np.array(inputs[str(case)]['initial'])
        result = replay(case,parameters,start,family,args)
        (args.out/f'case-{case}.json').write_text(json.dumps(result,indent=2)+'\n')
        print(case,result['mean_spectrum'],flush=True)


if __name__ == '__main__':
    main()
