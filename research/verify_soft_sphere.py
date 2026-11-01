"""Replay the declared sphere search and account for QR and rendering data."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from search_soft_sphere import screen, follow, step, checks, DEPENDENCIES
from follow_spherical_constellations import qr_step
from render_soft_sphere import orbit


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    set_num_threads(4); begin = time.time(); rows = []
    for name in DEPENDENCIES + ['render_soft_sphere.py', Path(__file__).name]:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    sources = {}
    for folder in [args.scout, args.follow, args.long, args.dense]:
        sources[str(folder)] = {p.name: sha(p) for p in sorted(folder.iterdir()) if p.is_file()}
    with np.load(args.scout/'screen.npz') as f:
        inputs = {k: f[k].copy() for k in ['params', 'points', 'counts', 'flags', 'base_ids']}
        differential = checks(inputs['params'], inputs['points'], inputs['counts'])
        values = screen(inputs['params'], inputs['points'], inputs['counts'], 5000, 15000)
        for key, value in zip(['rates', 'tails', 'status', 'minima', 'soft_stats', 'boundratios'], values):
            assert np.array_equal(value, f[key], equal_nan=True), ('screen', key)
        flags = (f['status'][:, 0] == 0) & (f['rates'].min(1) > .005)
        assert np.array_equal(flags, f['flags'])
    short = json.loads((args.follow/'report.json').read_text(encoding='utf-8'))
    assert short['ids'] == np.flatnonzero(inputs['flags']).tolist()
    short_summary = []
    for row in short['cases']:
        with np.load(args.follow/f"case-{row['id']}.npz") as f:
            assert f['status'].tolist() == row['status']
            assert np.array_equal(f['rates'].mean(0), row['spectrum_qr_order'])
            assert float(f['rates'][:, 0].min()) == row['minimum_leading_window']
            z = f['states']; scale = max(float(np.linalg.norm(np.ptp(z, axis=0))), 1.)
            residuals = [float(np.max(np.linalg.norm(z[-2048+lag:]-z[-2048:-lag], axis=1))) for lag in range(1, 513)]
            assert residuals == row['return_residuals']
            periods = [i+1 for i, r in enumerate(residuals) if r < 1e-8*scale]
            assert (min(periods) if periods else None) == row['detected_period']
            short_summary.append(dict(id=row['id'], positive_every_window=row['minimum_leading_window']>.005, detected_period=row['detected_period']))
    print(json.dumps(dict(stage='scout and complete follow-up accounting', seconds=time.time()-begin)), flush=True)
    report = json.loads((args.long/'report.json').read_text(encoding='utf-8'))
    ids = report['ids']; p = inputs['params'][ids]; sites = inputs['points'][ids]; counts = inputs['counts'][ids]
    values = follow(p, sites, counts, report['warm'], report['observed'])
    dense_report = json.loads((args.dense/'report.json').read_text(encoding='utf-8'))
    for j, ident in enumerate(ids):
        with np.load(args.long/f'case-{ident}.npz') as f:
            for key, value in zip(['rates', 'states', 'status', 'minima', 'soft_stats'], values):
                assert np.array_equal(value[j], f[key]), (ident, key)
        assert values[2][j].tolist() == [0, report['observed']]
        beta, tau = p[j]; points = sites[j, :counts[j]]; x = values[1][j, 0].copy()
        Q = np.eye(3); L = np.eye(3); s = np.zeros(3); t = np.zeros(3); volume = 0.
        for it in range(2048):
            x, J, r, peak, entropy, flag = step(x, beta, points, tau); assert not flag
            Q, a, flag = qr_step(J, Q); assert not flag
            L, R = np.linalg.qr(J@L); d = np.diag(R); L *= np.where(d<0, -1., 1.)
            s += a; t += np.log(abs(d)); sign, logdet = np.linalg.slogdet(J); assert sign != 0
            volume += logdet
        row = dict(id=ident, long_replay='bitwise equal', lapack_mean_rate_error=float(abs(s-t).max()/2048), mean_logdet_error=float(abs(s.sum()-volume)/2048), parameter_hex=[float(v).hex() for v in p[j]], point_hex=[[float(v).hex() for v in b] for b in points])
        assert row['lapack_mean_rate_error'] < 1e-8 and row['mean_logdet_error'] < 1e-8
        if ident in dense_report['ids']:
            meta = json.loads((args.dense/f'case-{ident}.json').read_text(encoding='utf-8'))
            z, minimum = orbit(beta, tau, points, dense_report['warm'], dense_report['count'])
            with np.load(args.dense/f'case-{ident}-states.npz') as f: assert np.array_equal(z, f['states'])
            xy = (z-np.array(meta['center']))@np.array(meta['basis'])@np.array(meta['projection']).T
            h = np.histogram2d(xy[:, 0], xy[:, 1], bins=1600, range=[[-meta['halfwidth'], meta['halfwidth']]]*2)[0].T[::-1].astype(np.uint32)
            with np.load(args.dense/f'case-{ident}-histogram.npz') as f: assert np.array_equal(h, f['histogram'])
            assert h.sum() == len(z) and minimum == meta['minimum_radius']
            row['dense_replay'] = 'all consecutive states and histogram bitwise equal'
        rows.append(row)
        print(json.dumps(dict(id=ident, seconds=time.time()-begin)), flush=True)
    result = dict(status='passed', source_files=sources, differential_checks=differential, scout_replay='429 cases, all output arrays bitwise equal', follow_accounting=short_summary, cases=rows, seconds=time.time()-begin, scope='Finite numerical reproduction and independent local arithmetic checks, not an interval certificate or historical priority claim.')
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['scout', 'follow', 'long', 'dense', 'out']: p.add_argument('--'+name, type=Path, required=True)
    main(p.parse_args())
