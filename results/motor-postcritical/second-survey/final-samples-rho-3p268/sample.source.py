"""Focused numerical follow-up of the sixfold transverse motor-cycle flip.

Exact published N=8 equations, original normalized time, fixed source gain.
Finite-time screens and periodic shooting are numerical evidence, not proofs.
Existing evidence is read-only. Each command refuses to overwrite its output.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from motor_rows import rhs, field_jac, source_start, unpack
from analyze_motor_cluster_cycle import Integrator, shoot, embedding

ROOT = Path(__file__).resolve().parents[1]
RHO_STAR = 3.3702183995790422
GAIN = 1 + .95/np.sqrt(.064002)
N = 8
D = 23


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def provenance():
    return {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), ROOT/'research/motor_rows.py',
                      ROOT/'research/analyze_motor_cluster_cycle.py')}


def write(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def prepare(folder):
    if folder.exists():
        raise FileExistsError('Preserve previous experiments; choose another output directory')
    folder.mkdir(parents=True)
    source = ROOT/'results/motor-rows/cluster-cycles/crossing-check.json'
    old = json.loads(source.read_text())
    clock = Integrator(600)
    offsets = [.02, .005, -.002, -.005, -.01, -.02, -.04, -.08, -.12, -.20, -.30, -.436127490488133]
    rng = np.random.default_rng(2026103010)
    contrasts = []
    for _ in range(3):
        q = rng.normal(size=(2, 4)); q -= q.mean(axis=1)[:, None]
        q /= np.linalg.norm(q, axis=1)[:, None]
        contrasts.append(q)
    contrasts.append(np.array([[-.5, -.5, .5, .5], [-.5, -.5, .5, .5]]))
    parameters, starts, records, cycles = [], [], [], []
    for index, offset in enumerate(offsets):
        rho = RHO_STAR+offset
        p = (rho, GAIN, 1.)
        guess, period, history = shoot(clock, old['initial_state'], old['period'], p)
        base = embedding(guess, p, 8)
        sol = clock.integrate(guess, period, p, augmented=True)
        eigvec = []
        eigenvalues = []
        for off in (30, 39):
            vals, vecs = np.linalg.eig(sol.y[off:off+9, -1].reshape(3, 3))
            j = np.argmin(vals.real)
            v = vecs[:, j].real
            v /= np.linalg.norm(v)
            if v[0] < 0:
                v *= -1
            eigvec.append(v); eigenvalues.append(float(vals[j].real))
        cycles.append(dict(rho=rho, offset=offset, period=period,
                           state=guess.tolist(), shooting=history,
                           transverse_flip_multipliers=eigenvalues))
        for seed_index, contrast in enumerate(contrasts):
            c, s, xx = unpack(base, 8)
            perturb = np.concatenate([contrast[k, :, None]*eigvec[k] for k in (0, 1)], axis=0)
            physical = np.c_[c, s, xx] + 1e-3*perturb
            start = np.r_[physical[:, 0], physical[:, 1], physical[:-1, 2]]
            assert np.max(physical[:, 0]**2+physical[:, 1]**2) < 2
            parameters.append(p); starts.append(start)
            records.append(dict(lane=len(starts)-1, rho_index=index,
                                initial_kind='generic_cycle_contrast' if seed_index < 3 else 'two_plus_two_in_each_cluster',
                                seed_index=seed_index, perturbation_amplitude=1e-3,
                                contrast=contrast.tolist(), unstable_row_vectors=[v.tolist() for v in eigvec]))
        start, positions = source_start(8, 2026103010)
        parameters.append(p); starts.append(start)
        records.append(dict(lane=len(starts)-1, rho_index=index,
                            initial_kind='source_style', seed=2026103010,
                            absolute_positions_nm=positions.tolist()))
        print('Prepared', index, rho, period, eigenvalues, flush=True)
    np.savez_compressed(folder/'initial.npz', parameters=np.array(parameters), initial_states=np.array(starts))
    write(folder/'design.json', dict(created_utc=now(), rho_star_reference=RHO_STAR,
          rho_star_status='Earlier approximate value; interval proof is documented separately',
          gain=GAIN, gain_exact_definition='1+(19/20)/sqrt(32001/500000)', b=1, N=8,
          parameter_variation='rho only; original normalized tau', cycles=cycles, lanes=records,
          source_file=source.relative_to(ROOT).as_posix(), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
          code_sha256=provenance(), solver_calls=clock.calls))
    (folder/'prepare.source.py').write_bytes(Path(__file__).read_bytes())


def survey(folder, warm=2000., observe=2000., wall=600., name='survey', previous=None):
    from motor_rows_gpu import MotorRowsGPU, save
    out = folder/(name+'.npz')
    if out.exists() or out.with_suffix('.json').exists():
        raise FileExistsError(out)
    source = folder/(previous+'.npz') if previous else folder/'initial.npz'
    old = np.load(source)
    starts = old['final_states'] if previous else old['initial_states']
    tangents = old['final_tangents'] if previous else None
    if previous and np.any(old['status'] != 0):
        raise ValueError('Continuation requires all previous lanes complete')
    result = MotorRowsGPU(8).scan(old['parameters'], starts, tangents=tangents,
             transient=warm, duration=observe, rtol=1e-10, step_factor=.5,
             max_steps=30000000, max_wall=wall, time_mode='original')
    result.update(design_file='design.json', source_file=source.relative_to(ROOT).as_posix(),
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  created_utc=now(), code_sha256=provenance(),
                  continuation_note='Every lane is carried without filtering' if previous else 'All declared starts included')
    save(out, result)
    (folder/(name+'.source.py')).write_bytes(Path(__file__).read_bytes())


def sample(folder, source_name='survey', duration=300., lanes=None, max_step=.04, name='samples'):
    out = folder/name
    if out.exists():
        raise FileExistsError(out)
    out.mkdir()
    source = folder/(source_name+'.npz')
    old = np.load(source)
    selected = list(range(len(old['parameters']))) if lanes is None else lanes
    summary = dict(created_utc=now(), source=source.relative_to(ROOT).as_posix(),
                   source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   code_sha256=provenance(), duration=duration, max_step=max_step,
                   tolerance=dict(rtol=2e-11, atol=2e-13), lanes=[])
    for lane in selected:
        began = time.perf_counter()
        if old['status'][lane] != 0:
            summary['lanes'].append(dict(lane=lane, status='skipped incomplete source'))
            continue
        p, initial = old['parameters'][lane], old['final_states'][lane]
        def section(t, y):
            return y[16]
        section.direction = 1
        times = np.linspace(0, duration, int(round(duration/.02))+1)
        sol = solve_ivp(lambda t, y: rhs(t, y, p, 8), (0, duration), initial,
                        method='DOP853', rtol=2e-11, atol=2e-13, max_step=max_step,
                        t_eval=times, events=section)
        states = sol.y.T
        section_states = sol.y_events[0]
        delta_t = np.diff(sol.t_events[0])
        # Return lag distances distinguish period 1,2,... on the X1=0 upward section.
        recurrence = []
        for lag in range(1, 17):
            if len(section_states) <= lag+5:
                continue
            distances = np.max(abs(section_states[lag:]-section_states[:-lag]), axis=1)
            recurrence.append(dict(lag=lag, max_last20=float(np.max(distances[-20:])),
                                   median_last20=float(np.median(distances[-20:]))))
        full_x = np.c_[states[:, 16:], -states[:, 16:].sum(axis=1)]
        pair_metrics = []
        for a in range(8):
            for b in range(a+1, 8):
                diff = (states[:, a]-states[:, b])**2+(states[:, 8+a]-states[:, 8+b])**2+(full_x[:, a]-full_x[:, b])**2
                pair_metrics.append(dict(rows=[a,b], rms=float(np.sqrt(np.mean(diff))), maximum=float(np.sqrt(np.max(diff)))))
        record = dict(lane=lane, parameters=p.tolist(), status='complete' if sol.success else 'failed',
                      message=sol.message, nfev=sol.nfev, section_count=len(section_states),
                      last_20_return_intervals=delta_t[-20:].tolist(), recurrence=recurrence,
                      row_pair_distances=pair_metrics,
                      max_disk_radius_squared=float(np.max(states[:, :8]**2+states[:, 8:16]**2)),
                      wall_seconds=time.perf_counter()-began)
        np.savez_compressed(out/f'lane-{lane:03d}.npz', parameters=p, time=sol.t, states=states,
                            section_time=sol.t_events[0], section_states=section_states,
                            initial_state=initial, final_state=states[-1])
        write(out/f'lane-{lane:03d}.json', record)
        summary['lanes'].append(record)
        write(out/'summary.json', summary)
        print(lane, record['status'], record['section_count'], 'lag minima',
              sorted([(x['max_last20'], x['lag']) for x in recurrence])[:3], flush=True)
    (out/'sample.source.py').write_bytes(Path(__file__).read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare', 'survey', 'sample'])
    parser.add_argument('--output', type=Path, default=ROOT/'results/motor-postcritical/initial-survey')
    parser.add_argument('--warm', type=float, default=2000.)
    parser.add_argument('--observe', type=float, default=2000.)
    parser.add_argument('--wall', type=float, default=600.)
    parser.add_argument('--name')
    parser.add_argument('--previous')
    parser.add_argument('--source', default='survey')
    parser.add_argument('--lanes', type=int, nargs='+')
    parser.add_argument('--duration', type=float, default=300.)
    parser.add_argument('--max-step', type=float, default=.04)
    args = parser.parse_args()
    args.output = args.output.resolve()
    if args.command == 'prepare':
        prepare(args.output)
    elif args.command == 'survey':
        survey(args.output, args.warm, args.observe, args.wall, args.name or 'survey', args.previous)
    else:
        sample(args.output, args.source, args.duration, args.lanes, args.max_step, args.name or 'samples')


if __name__ == '__main__':
    main()
