"""Actual discrete approaches to the two nonuniform posterior shell candidates."""
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable
from scipy.spatial import cKDTree
from numba import njit
from search_product_alphabet import step


@njit(cache=True)
def advance(starts, beta, tau, levels, order, count):
    states = np.empty((len(starts), count + 1, 3))
    for k in range(len(starts)):
        x = starts[k].copy()
        states[k, 0] = x
        for n in range(count):
            x, _, _, _, _, flag = step(x, beta, tau, levels, order)
            if flag:
                raise ValueError('Invalid iterate')
            states[k, n + 1] = x
    return states


def main(args):
    root = Path(__file__).resolve().parents[1]
    args.out.mkdir(parents=True, exist_ok=False)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    protocol = json.loads(args.cohort.read_text(encoding='utf-8'))
    starts = np.array(protocol['starts'])
    shown = np.unique(np.rint(np.geomspace(1, 2001, 151)).astype(int) - 1)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), gridspec_kw={'width_ratios': [1, 1, 1.15]})
    rows = []
    norm = LogNorm(1, 2001)
    for row, ident in enumerate([6, 7]):
        source = args.source / f'case-{ident}-states.npz'
        metadata = args.source / f'case-{ident}.json'
        with np.load(source) as f:
            beta, tau = f['params']; levels = f['levels']; order = int(f['order'])
            reference = f['states'][::20].copy()
        meta = json.loads(metadata.read_text(encoding='utf-8'))
        C = np.array(meta['basis']) @ np.array(meta['projection']).T
        center = np.array(meta['center'])
        states = advance(starts, beta, tau, levels, order, 2000)
        assert np.isfinite(states).all()
        # Independent direct weights and nested reflections at every start and final preimage.
        errors = []
        for k in range(len(starts)):
            for n in [0, 1999]:
                x = states[k, n]; m = np.empty(3)
                for j in range(3):
                    ell = levels[j, :order]
                    logits = -(x[j] - ell) ** 2 / tau
                    weights = np.exp(logits - logits.max()); weights /= weights.sum()
                    m[j] = weights @ ell
                reflected = 2 * m - x
                sphere_reflected = 2 * reflected / np.linalg.norm(reflected) - reflected
                y = beta / 2 * (sphere_reflected + x) + (1 - beta) * m
                errors.append(float(np.max(np.abs(y - states[k, n + 1]))))
        assert max(errors) < 1e-11
        scale = np.linalg.norm(np.ptp(reference, axis=0))
        distances = cKDTree(reference).query(states.reshape(-1, 3), workers=1)[0].reshape(33, 2001) / scale
        quantiles = np.quantile(distances, [.1, .5, .9], axis=0)
        xy = (reference - center) @ C
        p = (states[:, shown] - center) @ C
        for col in [0, 1]:
            ax = axes[row, col]
            ax.scatter(xy[:, 0], xy[:, 1], s=.12, c='#164e63', alpha=.3, rasterized=True)
            ax.scatter(p[:, :, 0].ravel(), p[:, :, 1].ravel(), c=np.tile(shown + 1, 33),
                       s=3.5, cmap='plasma_r', norm=norm, alpha=.7, edgecolors='none')
            if col == 1:
                half = float(np.abs(xy).max() * 1.25)
                ax.set(xlim=(-half, half), ylim=(-half, half))
            ax.set_aspect('equal')
            ax.set(title=f'Case {ident}: ' + ('33 initial states' if col == 0 else 'central zoom'),
                   xlabel='View coordinate 1', ylabel='View coordinate 2')
        ax = axes[row, 2]
        n = np.arange(2001)
        ax.fill_between(n, np.maximum(quantiles[0], 1e-12), quantiles[2], color='#70a9b7', alpha=.35, label='10th–90th percentile')
        ax.plot(n, quantiles[1], color='#9d4939', lw=1.2, label='Median')
        ax.set(xscale='symlog', yscale='log', xlabel='Iteration', ylabel='Distance / reference-cloud diagonal',
               title=f'Case {ident}: distance to reference points')
        ax.grid(alpha=.15); ax.legend(frameon=False, fontsize=8)
        compact = args.out / f'case-{ident}-approach.npz'
        np.savez_compressed(compact, starts=starts, iterations=shown, states=states[:, shown], distance_quantiles=quantiles)
        assert compact.stat().st_size < 256 * 1024
        rows.append(dict(id=ident, source_sha256=sha(source), metadata_sha256=sha(metadata),
                         compact_sha256=sha(compact), reference_states=len(reference), reference_stride=20,
                         independent_step_checks=len(errors), maximum_step_error=max(errors),
                         distance_quantiles_at_2000=quantiles[:, -1].tolist()))
    fig.subplots_adjust(left=.06, right=.98, top=.94, bottom=.16, wspace=.38, hspace=.42)
    cax = fig.add_axes([.08, .075, .50, .018])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap='plasma_r'), cax=cax, orientation='horizontal')
    cb.set_ticks([1, 2, 11, 101, 1001, 2001], labels=['0', '1', '10', '100', '1000', '2000'])
    cb.set_label('Iteration of the plotted discrete state (logarithmic sampling)')
    fig.text(.5, .018, 'Actual iterates, not continuous paths. Equal spatial units in each view; zooms clip remote points. Finite-cloud distance is not a proof of attraction.', ha='center', fontsize=8.5)
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, dpi=160, facecolor='white'); plt.close(fig)
    subprocess.run([sys.executable, str(root / 'tools/render_dates.py'), 'record',
                    args.figure.resolve().relative_to(root).as_posix(), '--role', 'diagnostic'], check=True)
    report = dict(date='2026-11-01', cases=rows, starts=33, updates=2000, displayed_iterations=shown.tolist(),
                  cohort_protocol_sha256=sha(args.cohort), figure_sha256=sha(args.figure),
                  source_script_sha256=sha(Path(__file__)),
                  scope='Finite approaches from the existing 33-start cohort to a stride-20, 100000-point reference. No basin, invariant-set, or asymptotic-attraction certificate.')
    (args.out / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(report))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['source', 'cohort', 'out', 'figure']:
        p.add_argument('--' + name, type=Path, required=True)
    main(p.parse_args())
