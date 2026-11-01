"""Independent algebra, solver, initial-state and parameter controls."""
from pathlib import Path
import argparse
import json
import numpy as np
from numba import njit, set_num_threads
from scipy.integrate import solve_ivp
from gaussian_vortex import field_jacobian, advance, integrate

ROOT = Path(__file__).resolve().parents[1]


def complex_field(z, p):
    f = -p[0]*z+p[8:11]
    for n in range(3):
        r = z-p[7]*np.eye(3)[n]
        f += p[1+n]*np.exp(-np.dot(r, r)/(2*p[4+n]**2))*np.cross(np.eye(3)[(n+1) % 3], r)
    return f


@njit
def independent_field(z, p):
    axes = np.eye(3)
    f = -p[0]*z + p[8:11]
    jac = -p[0]*np.eye(3)
    for n in range(3):
        r = z-p[7]*axes[n]
        axis = axes[(n+1) % 3]
        c = np.cross(axis, r)
        g = p[n+1]*np.exp(-np.dot(r, r)/(2*p[n+4]**2))
        cross_matrix = np.column_stack((np.cross(axis, axes[0]),
                                       np.cross(axis, axes[1]),
                                       np.cross(axis, axes[2])))
        f += g*c
        jac += g*(cross_matrix-np.outer(c, r)/p[n+4]**2)
    return f, jac


def adaptive_spectrum(p, warm=1000., observed=6000., block=5.):
    def rhs(t, state):
        f, j = independent_field(state[:3], p)
        return np.r_[f, (j @ state[3:].reshape(3, 3)).ravel()]
    z = np.array([.1, .2, .3])
    q = np.eye(3)
    logs = []
    for i in range(round((warm+observed)/block)):
        sol = solve_ivp(rhs, (0., block), np.r_[z, q.ravel()],
                        method='DOP853', rtol=2e-11, atol=2e-13)
        assert sol.success
        z = sol.y[:3, -1]
        q, r = np.linalg.qr(sol.y[3:, -1].reshape(3, 3))
        if i >= round(warm/block):
            logs.append(np.log(np.abs(np.diag(r))))
    logs = np.array(logs)
    windows = logs.reshape(-1, round(500/block), 3).sum(1)/500
    return {'method': 'DOP853', 'rtol': 2e-11, 'atol': 2e-13,
            'qr_interval': block, 'warm_time': warm, 'observed_time': observed,
            'window_time': 500., 'windows': windows.tolist(),
            'spectrum': windows.mean(0).tolist(), 'end': z.tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=ROOT/'results/gaussian-vortex')
    args = parser.parse_args()
    set_num_threads(4)
    p = np.array(json.loads((args.out/'candidates.json').read_text())['candidates'][0]['parameters'])
    rng = np.random.default_rng(20261110)
    residuals = {'field': 0., 'jacobian_complex_step': 0., 'divergence': 0., 'cyclic_equivariance': 0.}
    for z in rng.uniform(-4, 4, (100, 3)):
        f, j = field_jacobian(z, p)
        f2, _ = independent_field(z, p)
        jc = np.column_stack([complex_field(z.astype(complex)+1e-25j*e, p).imag/1e-25 for e in np.eye(3)])
        residuals['field'] = max(residuals['field'], float(np.max(np.abs(f-f2))))
        residuals['jacobian_complex_step'] = max(residuals['jacobian_complex_step'], float(np.max(np.abs(j-jc))))
        residuals['divergence'] = max(residuals['divergence'], float(abs(np.trace(j)+3*p[0])))
        residuals['cyclic_equivariance'] = max(residuals['cyclic_equivariance'], float(np.max(np.abs(field_jacobian(np.roll(z, 1), p)[0]-np.roll(f, 1)))))
    assert max(residuals.values()) < 2e-13
    # The vortex-free field has a closed-form solution and three exponents -gamma.
    linear = p.copy(); linear[1:4] = 0.
    z0 = np.array([.1, .2, .3]); z = z0.copy(); q = np.eye(3)
    for i in range(1000):
        z, q = advance(z, q, linear, .01)
    exact = linear[8:11]/linear[0]+(z0-linear[8:11]/linear[0])*np.exp(-linear[0]*10)
    residuals['linear_state'] = float(np.max(np.abs(z-exact)))
    residuals['linear_tangent'] = float(np.max(np.abs(q-np.eye(3)*np.exp(-linear[0]*10))))
    assert residuals['linear_state'] < 1e-11 and residuals['linear_tangent'] < 1e-11
    # Compare the local RK4 state against a separate adaptive implementation.
    sol = solve_ivp(lambda t, z: independent_field(z, p)[0], (0, 2.), z0,
                    method='DOP853', rtol=1e-12, atol=1e-14)
    state_errors = []
    for h in [.025, .0125, .00625]:
        z = z0.copy(); q = np.eye(3)
        for i in range(round(2/h)):
            z, q = advance(z, q, p, h)
        state_errors.append(float(np.linalg.norm(z-sol.y[:, -1])))
    assert state_errors[2] < state_errors[1]/10 < state_errors[0]/100
    starts = rng.uniform(-3, 3, (16, 3))
    parameters = np.tile(p, (16, 1))
    w, tail, end, means, radius = integrate(parameters, starts, np.full(16, .0125), 1000., 6000., 500., 10000, 4)
    np.savez_compressed(args.out/'fresh-starts.npz', parameters=parameters, starts=starts,
                        stepsize=.0125, warm_time=1000., observed_time=6000., window_time=500.,
                        tail_stride=4, windows=w, tail=tail, end=end, means=means, max_radius=radius)
    print('fresh starts', w.mean(1)[:, 0].tolist(), flush=True)
    # Small general perturbations break the exact C3 symmetry; also test rounding.
    neighbors = np.tile(p, (13, 1))
    neighbors[:12] *= 1+rng.uniform(-.005, .005, (12, 11))
    neighbors[12] = np.round(p, 2)
    starts2 = np.tile(z0, (13, 1))
    w2, tail2, end2, means2, radius2 = integrate(neighbors, starts2, np.full(13, .0125), 1000., 6000., 500., 10000, 4)
    np.savez_compressed(args.out/'neighbors.npz', parameters=neighbors, starts=starts2,
                        stepsize=.0125, warm_time=1000., observed_time=6000., window_time=500.,
                        tail_stride=4, windows=w2, tail=tail2, end=end2, means=means2, max_radius=radius2)
    print('neighbors', w2.mean(1)[:, 0].tolist(), flush=True)
    adaptive = adaptive_spectrum(p)
    result = {'algebra_max_abs_residuals': residuals,
              'local_rk4_errors_h_025_0125_00625': state_errors,
              'independent_adaptive': adaptive,
              'fresh_start_seed': 20261110, 'fresh_start_box': [-3, 3],
              'neighbor_relative_box': [-.005, .005], 'rounded_neighbor_index': 12,
              'fresh_spectra': w.mean(1).tolist(), 'neighbor_spectra': w2.mean(1).tolist(),
              'bounded_nonlinear_speed': float(np.sum(np.abs(p[1:4])*p[4:7])*np.exp(-.5)),
              'absorbing_radius_threshold': float((np.linalg.norm(p[8:11])+np.sum(np.abs(p[1:4])*p[4:7])*np.exp(-.5))/p[0])}
    (args.out/'controls.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'algebra': residuals, 'adaptive_spectrum': adaptive['spectrum']}), flush=True)


if __name__ == '__main__':
    main()
