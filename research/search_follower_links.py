"""Predeclared three-link follower-force pilot; retain every outcome."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from follower_link import pilot, matrices, field, augmented3


def write(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n', encoding='utf8')


def algebra_controls():
    rng = np.random.default_rng(17346)
    worst = 0.
    smallest = 1e99
    for _ in range(20):
        x = rng.normal(size=3)
        v = rng.normal(size=3)
        k, drag, length = np.exp(rng.uniform(-1, 1, (3, 3)))
        force = 7.
        m, _ = matrices(x, force, k, drag, length)
        difference = np.eye(3)-np.eye(3, k=-1)
        elastic = difference.T@np.diag(k)@difference
        tail = np.cumsum(drag[::-1])[::-1]
        base = tail[np.maximum.outer(np.arange(3), np.arange(3))]*np.outer(length, length)
        got = augmented3(np.r_[x, v], force, elastic, base, length)
        f = field(x, force, k, drag, length)
        jv = field(x+1e-25j*v, force, k, drag, length).imag/1e-25
        worst = max(worst, float(np.max(abs(got-np.r_[f, jv]))))
        smallest = min(smallest, float(np.linalg.eigvalsh(m).min()))
    # Source N=2 characteristic polynomial: lambda^2+2(3-Sigma)lambda+1.
    polyerror = 0.
    for force in [1., 2.5, 3., 3.5, 5.]:
        j = np.column_stack([field(1e-25j*np.eye(2)[i], force, np.ones(2), np.ones(2), np.ones(2)).imag/1e-25 for i in range(2)])
        polyerror = max(polyerror, float(np.max(abs(np.poly(j)-[1., 2*(3-force), 1.]))))
    assert worst < 1e-10 and smallest > 0 and polyerror < 1e-12
    return dict(directional_derivative_max_error=worst, smallest_sampled_resistance_eigenvalue=smallest,
                source_two_link_characteristic_polynomial_error=polyerror)


def run(job):
    started = time.perf_counter()
    trace, lle, windows, final, minimum, maximum = pilot(job['force'], np.array(job['stiffness']),
        np.ones(3), np.ones(3), np.array(job['initial']), job['step'], job['warm'], job['duration'], job['stride'])
    path = Path(job['output'])/f"case-{job['id']:04}.npz"
    if not np.isfinite(lle):
        result = dict(id=job['id'], status='nonfinite; rejected', job=job)
    else:
        np.savez_compressed(path, states=trace, final=final, force=job['force'], stiffness=job['stiffness'],
                            step=job['step'], sample_step=job['step']*job['stride'], warm=job['warm'], duration=job['duration'])
        result = dict(id=job['id'], force=job['force'], stiffness=job['stiffness'], leading_rate=float(lle),
            windows=windows.tolist(), min_angle_norm=float(minimum), max_absolute_angle=float(maximum),
            file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), wall_seconds=time.perf_counter()-started,
            status='fixed-step finite-time pilot; no attractor or novelty claim')
    write(path.with_suffix('.json'), result)
    return result


def main(args):
    out=args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    controls=algebra_controls()
    jobs=[]
    forces=np.geomspace(.6, 60., 48)
    springs=[[1.,1.,1.]] if args.equal else [[1., a, b] for a in [.3, 1., 3.] for b in [.3, 1., 3.]]
    for k in springs:
        for force in forces:
            jobs.append(dict(id=len(jobs),force=float(force),stiffness=k,initial=[.1,.02,-.03],
                             step=args.step,warm=200.,duration=800.,stride=int(round(.05/args.step)),output=str(out)))
    write(out/'protocol.json',dict(jobs=jobs, controls=controls,
        source='De Canio, Lauga, Goldstein (2017), arXiv:1707.02346, Sec. IV',
        extension='Three links instead of two; unit lengths and bead drags; elastic joint ratios varied as declared.',
        prior_warning='Chaos in continuous follower-force filaments is already known (Ling2018, Clarke2024).',
        method='RK4 state plus one analytic variational direction; .1 tau renormalization; four observed windows.',
        limitations='No hydrodynamic coupling, self-avoidance, contact or real-plant claim. Screen needs independent integration.'))
    for filename in ['follower_link.py','search_follower_links.py']:
        (out/filename).write_bytes((Path(__file__).parent/filename).read_bytes())
    records=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            record=future.result();records.append(record)
            write(out/'analysis.json',dict(completed=len(records),declared=len(jobs),runs=sorted(records,key=lambda r:r['id'])))
            print(record['id'], record.get('leading_rate'), record.get('windows'), flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=3);p.add_argument('--step',type=float,default=.002)
    p.add_argument('--equal',action='store_true');main(p.parse_args())
