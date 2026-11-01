"""Independent arbitrary-precision, finite-separation Adamax checks."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np


def run_one(job):
    ident, raw_parameters, digits, output = job
    mp.mp.dps = digits
    b, c, e = [mp.mpf(float(value)) for value in raw_parameters]
    delta = mp.mpf('1e-30' if digits == 70 else '1e-40')

    def advance(z, correction=1):
        x, m, u = z
        next_m = b*m + (1-b)*x
        next_u = max(c*u, abs(x)+e)
        return (x-next_m/(correction*next_u), next_m, next_u)

    z = (mp.mpf(float(.73)), mp.mpf(0), mp.mpf(0))
    power = mp.mpf(1)
    before = time.perf_counter()
    for _ in range(100000):
        power *= b
        z = advance(z, 1-power)
    nearby = (z[0]+delta, z[1], z[2])
    totals = []
    total = mp.mpf(0)
    tail = []
    resets = 0
    min_margin = mp.inf
    for t in range(105000):
        if t >= 5000:
            margin = abs(z[0])+e-c*z[2]
            resets += int(margin >= 0)
            min_margin = min(min_margin, abs(margin)/(1+z[2]))
        z = advance(z)
        nearby = advance(nearby)
        difference = [v-w for v,w in zip(nearby,z)]
        distance = mp.sqrt(sum(v*v for v in difference))
        if not distance:
            raise RuntimeError('Finite separation collapsed at arbitrary precision')
        nearby = tuple(w+delta*v/distance for w,v in zip(z,difference))
        if t >= 5000:
            total += mp.log(distance/delta)
            n = t-5000+1
            if n % 10000 == 0:
                totals.append(float(total/10000))
                total = mp.mpf(0)
            if n > 99000:
                tail.append([float(v) for v in z])
    row = dict(id=int(ident), digits=digits,
        parameters_binary64=[float(v) for v in raw_parameters],
        parameters_exact_promoted=[mp.nstr(v,digits) for v in (b,c,e)],
        separation=mp.nstr(delta,digits), leading_windows=totals,
        leading_mean=float(np.mean(totals)), minimum_window=min(totals),
        reset_count=resets, minimum_scaled_switch_margin=float(min_margin),
        final_state=[mp.nstr(v,digits) for v in z], seconds=time.perf_counter()-before)
    stem = Path(output)/f'case-{ident}-{digits}digits'
    np.savez_compressed(stem.with_suffix('.npz'), tail=np.array(tail), leading_windows=np.array(totals))
    row['tail_sha256'] = hashlib.sha256(stem.with_suffix('.npz').read_bytes()).hexdigest()
    stem.with_suffix('.json').write_text(json.dumps(row,indent=2)+'\n')
    return row


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source)
    jobs=[(ident,source['parameters'][ident].tolist(),digits,str(out)) for ident in args.ids for digits in (70,90)]
    protocol=dict(source=str(args.source.resolve()),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        ids=args.ids,digits=[70,90],cold_start=[.73,0,0],bias_corrected_warmup=100000,
        autonomous_alignment=5000,autonomous_observed=100000,window=10000,
        method='Independent scalar max recurrence, Euclidean finite separation renormalized every step; no analytic Jacobian or production map import.',
        parameter_semantics='Exact promotion of the selected binary64 coefficients; same real coefficients at both precisions.',
        selection='Visually distinct 847 and 4099, nearby crown 239, and standard default pair 4096; all passed the preceding million-step check.',
        scope='Finite-time numerical check, not proof of asymptotic chaos or novelty.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(run_one,job) for job in jobs]
        for future in as_completed(futures):
            row=future.result();rows.append(row)
            print(row['id'],row['digits'],row['leading_mean'],row['minimum_window'],flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=sorted(rows,key=lambda r:(r['id'],r['digits'])),completed=len(rows),declared=len(jobs)),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--ids',type=int,nargs='+',required=True)
    main(parser.parse_args())
