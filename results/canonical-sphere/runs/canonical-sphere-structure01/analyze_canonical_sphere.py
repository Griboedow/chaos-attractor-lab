"""Finite basin, full-spectrum, and posterior-only controls for canonical sites."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads
from search_soft_sphere import step, DEPENDENCIES
from follow_spherical_constellations import qr_step


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p, v): p.write_text(json.dumps(v, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def period(z):
    scale = max(float(np.linalg.norm(np.ptp(z, axis=0))), 1.)
    residuals = [float(np.max(np.linalg.norm(z[-2048+k:]-z[-2048:-k], axis=1))) for k in range(1, 513)]
    found = [k+1 for k, v in enumerate(residuals) if v < 1e-8*scale]
    return min(found) if found else None


@njit(cache=True, parallel=True)
def full_cohort(beta, tau, points, starts):
    rates = np.full((len(starts), 20, 3), np.nan)
    tails = np.full((len(starts), 4096, 3), np.nan)
    status = np.zeros((len(starts), 2), np.int64)
    for k in prange(len(starts)):
        x = starts[k].copy(); Q = np.eye(3); sums = np.zeros((20, 3))
        for it in range(600000):
            x, J, r, peak, entropy, flag = step(x, beta, points, tau)
            if flag or not np.isfinite(x).all(): status[k, 0] = flag if flag else 3; break
            Q, logs, flag = qr_step(J, Q)
            if flag: status[k, 0] = 4; break
            if it >= 100000: sums[(it-100000)//25000] += logs; status[k, 1] += 1
            if it >= 600000-4096: tails[k, it-(600000-4096)] = x
        if status[k, 0] == 0: rates[k] = sums/25000
    return rates, tails, status


def main(args):
    args.out.mkdir(parents=True, exist_ok=False); set_num_threads(4)
    for name in [Path(__file__).name]+DEPENDENCIES:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.scout/'screen.npz') as f:
        params = f['params'].copy(); points = f['points'].copy(); counts = f['counts'].copy()
    protocol = json.loads((args.cohort/'protocol.json').read_text(encoding='utf-8'))
    starts = np.array(protocol['starts']); rows = []
    for ident in protocol['ids']:
        with np.load(args.cohort/f'case-{ident}.npz') as f:
            periods = [period(z) for z in f['tails']]
            minima = f['rates'].min(1).tolist()
        beta, tau = params[ident]; sites = points[ident, :counts[ident]]; mean_controls = []
        for start in starts:
            x = start.copy()
            for k in range(10000):
                logits = -np.sum((sites-x)**2, axis=1)/tau
                weights = np.exp(logits-logits.max()); weights /= weights.sum()
                y = weights@sites; residual = float(np.linalg.norm(y-x)); x = y
                if residual < 1e-12: break
            mean_controls.append(dict(steps=k+1, residual=residual, stopped=residual < 1e-12))
        rows.append(dict(id=ident, periods=periods, minimum_windows=minima, mean_only=mean_controls))
    full = []; begin = time.time()
    for ident in [148, 653]:
        beta, tau = params[ident]; sites = points[ident, :counts[ident]]
        rates, tails, status = full_cohort(beta, tau, sites, starts)
        assert np.all(status[:, 0] == 0) and np.all(status[:, 1] == 500000)
        np.savez_compressed(args.out/f'full-{ident}.npz', starts=starts, params=params[ident], points=sites, rates=rates, tails=tails, status=status)
        row = dict(id=ident, mean_spectra=rates.mean(1).tolist(), minimum_windows=rates.min(1).tolist(), periods=[period(z) for z in tails], all_leading_windows_positive=int((rates[:, :, 0].min(1)>0).sum()), all_second_windows_positive=int((rates[:, :, 1].min(1)>0).sum()), minimum_second_window=float(rates[:, :, 1].min()))
        full.append(row); save(args.out/f'full-{ident}.json', row)
        print(json.dumps({k:v for k,v in row.items() if k not in ['mean_spectra','minimum_windows','periods']}), flush=True)
    save(args.out/'report.json', dict(source_sha256=sha(args.scout/'screen.npz'), cohort_source_sha256=sha(args.cohort/'report.json'), cohort_controls=rows, full_spectrum=dict(ids=[148,653], starts=starts.tolist(), warm=100000, observed=500000, windows=20, cases=full), seconds=time.time()-begin, scope='Finite sample controls, not a global basin, invariant measure, or asymptotic hyperchaos certificate.'))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['scout', 'cohort', 'out']: p.add_argument('--'+name, type=Path, required=True)
    main(p.parse_args())
