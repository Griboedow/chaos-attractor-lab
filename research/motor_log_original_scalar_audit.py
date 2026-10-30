"""Original-seed scalar audit with an explicit20tau integration horizon.

All algebra, Jacobian, gauge, and underflow checks are copied from the frozen
audit. The original200tau attempt failed and is preserved separately. The
imported scalar RHS is unchanged. This shorter guard does not validate long
pointwise trajectory agreement.
"""

import argparse
import cmath
import datetime
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp


from motor_log_scalar import rhs, push_field
import motor_log_scalar as evaluated_scalar


def audit(args):
    import motor_log_contrasts as log
    import motor_pair_coordinates as pair
    folder = args.output.resolve()
    if folder.exists():
        raise FileExistsError(folder)
    folder.mkdir(parents=True)
    data = np.load(args.source/'trajectory.npz')
    p = data['parameters']
    z = data['final_log_state'].copy() if 'final_log_state' in data else log.to_log(data['final_chart_state'])
    rng = np.random.default_rng(2026103053)
    errors = dict(field_scaled=0., jacobian_scaled=0., pushed_field_scaled=0.,
                  original_field_derivative_scaled=0., gauge_derivative_scaled=0.)
    for _ in range(32):
        chart = rng.normal(size=23)
        q = log.to_log(chart)
        gauge = np.exp(rng.uniform(-2, 2, 4))
        for offset in [11, 15, 19]:
            q[offset:offset+4] *= gauge
        q[23:27] -= np.log(gauge)
        f, baseline = rhs(0., q, p), log.rhs(0., q, p)
        errors['field_scaled'] = max(errors['field_scaled'], float(np.max(abs(f-baseline))/(1+np.max(abs(baseline)))))
        pushed = push_field(q, f)
        original = pair.EINV@pair.full_field_jac(pair.to_full(log.to_chart(q)), p, 8)[0]
        errors['pushed_field_scaled'] = max(errors['pushed_field_scaled'], float(np.max(abs(pushed-original))/(1+np.max(abs(original)))))
        h = 1e-28
        js = np.column_stack([np.imag(rhs(0., q.astype(complex)+1j*h*v, p))/h for v in np.eye(27)])
        jn = np.column_stack([np.imag(log.rhs(0., q.astype(complex)+1j*h*v, p))/h for v in np.eye(27)])
        errors['jacobian_scaled'] = max(errors['jacobian_scaled'], float(np.max(abs(js-jn))/(1+np.max(abs(jn)))))
        jp = np.column_stack([np.imag(push_field(q.astype(complex)+1j*h*v, rhs(0., q.astype(complex)+1j*h*v, p)))/h for v in np.eye(27)])
        reconstruction = np.zeros((23, 27))
        reconstruction[:11, :11] = np.eye(11)
        for offset in [11, 15, 19]:
            for j in range(4):
                reconstruction[offset+j, offset+j] = math.exp(q[23+j])
                reconstruction[offset+j, 23+j] = math.exp(q[23+j])*q[offset+j]
        jf = pair.full_field_jac(pair.to_full(log.to_chart(q)), p, 8)[1]
        expected_derivative = pair.EINV@jf@pair.E@reconstruction
        errors['original_field_derivative_scaled'] = max(errors['original_field_derivative_scaled'], float(np.max(abs(jp-expected_derivative))/(1+np.max(abs(expected_derivative)))))
        gauge_derivative = 2*(q[11:15]*f[11:15]+q[15:19]*f[15:19]+q[19:23]*f[19:23])
        errors['gauge_derivative_scaled'] = max(errors['gauge_derivative_scaled'], float(np.max(abs(gauge_derivative))/(1+np.max(log.gauge_norm_squared(q)))))
    if max(errors.values()) > 1e-12:
        raise ArithmeticError(errors)
    underflows = []
    for ell in [-370., -372.5, -372.6, -500., -1000.]:
        q = z.copy()
        q[23:27] = ell
        first, second = log.UnderflowTracker(), log.UnderflowTracker()
        a, b = rhs(0., q, p, first), log.rhs(0., q, p, second)
        assert first.records == second.records
        underflows.append(dict(ell=ell, maximum_field_difference=float(np.max(abs(a-b))),
                               exactly_equal_underflow_records=True, records=first.records))
    timing = []
    for name, field in [('numpy', log.rhs), ('scalar', rhs)]:
        tracker = log.UnderflowTracker()
        q = z.copy()
        q[24] = -500.
        started = time.perf_counter()
        for _ in range(20000):
            field(0., q, p, tracker)
        timing.append(dict(field=name, calls=20000, wall_seconds=time.perf_counter()-started,
                           tracking_active_underflow=True))
    start = float(data['time'][-1])
    checkpoints = np.linspace(start, start+20., 201)
    kwargs = dict(method='DOP853', max_step=.025, rtol=2e-12, atol=2e-14, t_eval=checkpoints)
    first_tracker, second_tracker = log.UnderflowTracker(), log.UnderflowTracker()
    began = time.perf_counter()
    first = solve_ivp(lambda t, q: log.rhs(t, q, p, first_tracker), (start, start+20.), z, **kwargs)
    first_wall = time.perf_counter()-began
    began = time.perf_counter()
    second = solve_ivp(lambda t, q: rhs(t, q, p, second_tracker), (start, start+20.), z, **kwargs)
    second_wall = time.perf_counter()-began
    assert first.success and second.success
    a, b = first.y.T, second.y.T
    control = dict(duration=20., maximum_mean_difference=float(np.max(abs(a[:, :11]-b[:, :11]))),
                   maximum_log10_norm_difference=float(np.max(abs(log.log10_norms(a)-log.log10_norms(b)))),
                   maximum_log_state_difference=float(np.max(abs(a-b))),
                   nfev=[first.nfev, second.nfev], wall_seconds=[first_wall, second_wall],
                   full_integrator_speed_ratio=first_wall/second_wall,
                   underflow_tracking_active=True)
    if max(control[k] for k in ['maximum_mean_difference', 'maximum_log10_norm_difference', 'maximum_log_state_difference']) > 1e-7:
        raise ArithmeticError(control)
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='passed',
                  source_data=str(args.source/'trajectory.npz'), source_data_sha256=pair.sha(args.source/'trajectory.npz'),
                  code_sha256=pair.sha(evaluated_scalar.__file__), audit_wrapper_sha256=pair.sha(__file__), log_code_sha256=pair.sha(log.__file__), pair_code_sha256=pair.sha(pair.__file__),
                  arbitrary_nonunit_gauges=32, errors=errors, underflow_checks=underflows,
                  timing=timing, rhs_speed_ratio=timing[0]['wall_seconds']/timing[1]['wall_seconds'],
                  short_numpy_scalar_control=control,
                  scope='Same real equations; independent scalar evaluation, analytic reconstruction derivative checked against original Cartesian Jacobian. Timing does not promise the same full-integrator speedup.')
    pair.write(folder/'audit.json', record)
    for name, module in [('audit-wrapper.py', __file__), ('source.py', evaluated_scalar.__file__), ('log-source.py', log.__file__), ('pair-source.py', pair.__file__)]:
        (folder/name).write_bytes(Path(module).read_bytes())
    print(__import__('json').dumps(record, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    audit(parser.parse_args())
