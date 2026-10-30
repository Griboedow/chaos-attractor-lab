"""Numerical periodic shooting in exact weighted motor-cluster restrictions.

Cluster weights are integer row counts divided by eight. The restriction is
checked against the full published N=8 field. No changed coupling is used.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from motor_rows import rhs as full_rhs, field_jac as full_field_jac

ROOT = Path(__file__).resolve().parents[1]


def unpack(state, weights):
    m = len(weights)
    c, s = state[:m], state[m:2*m]
    x = np.r_[state[2*m:], -np.dot(weights[:-1], state[2*m:])/weights[-1]]
    return c, s, x


def field_jac(state, parameters, weights):
    rho, r, b = parameters
    weights = np.asarray(weights)
    m = len(weights)
    c, s, x = unpack(state, weights)
    v = -rho*x+r*(s-np.dot(weights, s))
    field = np.r_[-c-v*s+1, -s+v*c+b, v[:-1]]
    jac = np.zeros((3*m-1, 3*m-1), dtype=np.result_type(field))
    eye = np.eye(m)
    projection = eye-np.tile(weights, (m, 1))
    embedding = np.vstack([np.eye(m-1), -weights[:-1]/weights[-1]])
    jac[:m, :m] = -eye
    jac[:m, m:2*m] = -np.diag(v)-r*s[:, None]*projection
    jac[:m, 2*m:] = rho*s[:, None]*embedding
    jac[m:2*m, :m] = np.diag(v)
    jac[m:2*m, m:2*m] = -eye+r*c[:, None]*projection
    jac[m:2*m, 2*m:] = -rho*c[:, None]*embedding
    jac[2*m:, m:2*m] = r*projection[:-1]
    jac[2*m:, 2*m:] = -rho*np.eye(m-1)
    return field, jac


def lift(state, groups):
    weights = np.array([len(q) for q in groups])/8
    c, s, x = unpack(state, weights)
    a = np.zeros((8, 3), dtype=np.result_type(state))
    for j, group in enumerate(groups):
        a[group] = c[j], s[j], x[j]
    return np.r_[a[:, 0], a[:, 1], a[:-1, 2]]


def restrict(state, groups):
    x = np.r_[state[16:], -np.sum(state[16:])]
    return np.r_[[np.mean(state[q]) for q in groups],
                 [np.mean(state[np.array(q)+8]) for q in groups],
                 [np.mean(x[q]) for q in groups[:-1]]]


def audit(groups, parameters):
    w = np.array([len(q) for q in groups])/8
    d = 3*len(w)-1
    e = np.column_stack([lift(q, groups) for q in np.eye(d)])
    rng = np.random.default_rng(2026103011)
    errors = dict(field=0., tangent=0., jacobian_complex=0.)
    for _ in range(12):
        x = rng.normal(size=d)
        f, j = field_jac(x, parameters, w)
        full, jf = full_field_jac(lift(x, groups), parameters, 8)
        errors['field'] = max(errors['field'], float(np.max(abs(e@f-full))))
        errors['tangent'] = max(errors['tangent'], float(np.max(abs(e@j-jf@e))))
        jc = np.column_stack([np.imag(field_jac(x.astype(complex)+1e-27j*q, parameters, w)[0])/1e-27 for q in np.eye(d)])
        errors['jacobian_complex'] = max(errors['jacobian_complex'], float(np.max(abs(j-jc))))
    if max(errors.values()) > 1e-10:
        raise ArithmeticError(errors)
    return errors


def integrate(state, period, p, weights, augmented=False, method='DOP853', rtol=2e-11, max_step=.025, dense=False):
    d = len(state)
    def fun(t, z):
        f, j = field_jac(z[:d], p, weights)
        return np.r_[f, (j@z[d:].reshape(d, d)).ravel()] if augmented else f
    z = np.r_[state, np.eye(d).ravel()] if augmented else state
    sol = solve_ivp(fun, (0, period), z, method=method, rtol=rtol, atol=rtol*.01,
                    max_step=max_step, dense_output=dense)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def shoot(state, period, p, weights, iterations=12):
    state = state.copy()
    d = len(state)
    phase = 2*len(weights)
    state[phase] = 0.
    free = np.delete(np.arange(d), phase)
    history = []
    for iteration in range(iterations):
        sol = integrate(state, period, p, weights, augmented=True)
        end = sol.y[:d, -1]
        f, j = field_jac(end, p, weights)
        m = sol.y[d:, -1].reshape(d, d)
        residual = end-state
        err = float(np.max(abs(residual)))
        history.append(dict(iteration=iteration, period=float(period), residual=err))
        if err < 2e-11:
            return state, float(period), sol, history
        a = np.column_stack([(m-np.eye(d))[:, free], f])
        correction = np.linalg.solve(a, -residual)
        if np.max(abs(correction[:-1])) > .5 or abs(correction[-1]) > period*.25:
            raise ArithmeticError('Newton correction exceeds bounded step: '+str(history))
        state[free] += correction[:-1]
        period += correction[-1]
    raise ArithmeticError('Shooting not converged: '+str(history))


def run(sample_path, groups, lag, output):
    sample_path, output = sample_path.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    source = np.load(sample_path)
    p = source['parameters']; weights = np.array([len(q) for q in groups])/8
    initial = restrict(source['section_states'][-1-lag], groups)
    period = source['section_time'][-1]-source['section_time'][-1-lag]
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), groups=groups,
                  weights=weights.tolist(), parameters=p.tolist(), section_lag=lag,
                  source=sample_path.relative_to(ROOT).as_posix(), source_sha256=hashlib.sha256(sample_path.read_bytes()).hexdigest(),
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  audit=audit(groups, p), status='running')
    began = time.perf_counter()
    try:
        state, period, sol, history = shoot(initial, period, p, weights)
        d = len(state)
        mon = sol.y[d:, -1].reshape(d, d)
        evals = np.linalg.eigvals(mon)
        # Independent implicit solver for cluster state and complete quotient monodromy.
        implicit = integrate(state, period, p, weights, augmented=True, method='Radau', rtol=2e-12, max_step=.0125)
        sample = integrate(state, period, p, weights, augmented=True, rtol=1e-12, max_step=.0125, dense=True)
        times = np.linspace(0, period, 1001); states = sample.sol(times)[:d].T
        lifted = np.array([lift(v, groups) for v in states])
        # Full independent 23D variational integration verifies restrictions and omitted directions.
        def fullaug(t, z):
            f, j = full_field_jac(z[:23], p, 8)
            return np.r_[f, (j@z[23:].reshape(23, 23)).ravel()]
        full = solve_ivp(fullaug, (0, period), np.r_[lift(state, groups), np.eye(23).ravel()],
                         method='DOP853', rtol=2e-12, atol=2e-14, max_step=.0125)
        if not full.success:
            raise RuntimeError(full.message)
        fm = full.y[23:, -1].reshape(23, 23)
        full_evals = np.linalg.eigvals(fm)
        order = np.argsort(-abs(full_evals)); full_evals = full_evals[order]
        orderc = np.argsort(-abs(evals)); evals = evals[orderc]
        record.update(status='complete', period=period, initial_state=state.tolist(), shooting_history=history,
            quotient_multipliers=[dict(real=float(v.real), imag=float(v.imag), modulus=float(abs(v))) for v in evals],
            full_multipliers=[dict(real=float(v.real), imag=float(v.imag), modulus=float(abs(v))) for v in full_evals],
            full_exponents=(np.log(abs(full_evals))/period).tolist(),
            full_unstable_count_at_1e_minus_6=int(np.sum(abs(full_evals)>1+1e-6)),
            full_neutral_count_at_1e_minus_6=int(np.sum(abs(abs(full_evals)-1)<1e-6)),
            closure_error=float(np.max(abs(sol.y[:d, -1]-state))),
            Radau_closure_error=float(np.max(abs(implicit.y[:d, -1]-state))),
            DOP853_Radau_variational_difference=float(np.max(abs(sol.y[d:, -1]-implicit.y[d:, -1]))),
            full_closure_error=float(np.max(abs(full.y[:23,-1]-lift(state, groups)))),
            full_flow_identity_error=float(np.max(abs(fm@full_rhs(0,lift(state, groups),p,8)-full_rhs(0,lift(state, groups),p,8)))),
            caveat='Floating-point periodic shooting and multiplier estimates; not an interval certificate.')
        np.savez_compressed(output/'cycle.npz', parameters=p, time=times, cluster_states=states,
                            full_states=lifted, quotient_monodromy=mon, full_monodromy=fm,
                            quotient_initial_state=state)
    except Exception as error:
        record.update(status='failed', error=repr(error))
        raise
    finally:
        record['wall_seconds']=time.perf_counter()-began
        (output/'analysis.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
        (output/'source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record, indent=2), flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', type=Path, required=True)
    parser.add_argument('--groups', required=True, help='JSON list of zero-based row-index groups')
    parser.add_argument('--lag', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    a=parser.parse_args();run(a.sample, json.loads(a.groups), a.lag, a.output)
