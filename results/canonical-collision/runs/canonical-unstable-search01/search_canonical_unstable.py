"""Search for a short reflection-plane crossing from a saddle's unstable direction.

The starting line is not asserted to be the exact unstable manifold. A later
interval graph-transform enclosure is required before using a crossing as proof.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from canonical_cube import parameters, model


def main(args):
    sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps = 250
    args.out.mkdir(parents=True, exist_ok=False)
    for name in [Path(__file__).name, 'canonical_cube.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    par = parameters(args.source, 653, mp.mp)
    def equation(s, t):
        x = [s, s, t]; y = model(x, par, mp.mp, False)
        return y[0]-s, y[2]-t
    s, t = mp.findroot(equation, ('-.5744', '-.1238'), tol=mp.mpf('1e-235'))
    P = [s, s, t]; _, J, radius = model(P, par, mp.mp)
    lam = J[0][0]-J[0][1]
    dec = lambda v: mp.nstr(v, 260)
    protocol = dict(source=str(args.source), source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(), id=653, precision=250, fixed_point=[dec(v) for v in P], multiplier=dec(lam), initial_direction=[1,-1,0], epsilon='1e-80', phases=args.phases, maximum_steps=1200, selection='First x-y sign reversal, ranked by number of steps after distance from the saddle exceeds 0.01.', scope='Numerical guide only: the initial line differs from the true unstable manifold by an unvalidated higher-order term.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    rows = []; begin = time.time()
    for k in range(args.phases):
        eps = mp.mpf('1e-80')*(1+(lam-1)*k/(args.phases-1))
        x = [s+eps,s-eps,t]; escape = None; row = dict(phase=k, epsilon=dec(eps))
        for n in range(1, 1201):
            y = model(x, par, mp.mp, False)
            if escape is None and mp.sqrt(sum((a-b)**2 for a,b in zip(y,P))) > mp.mpf('.01'): escape=n
            if y[0]-y[1] < 0:
                row.update(crossing_iteration=n, escape_iteration=escape, steps_after_escape=n-escape, before=[dec(v) for v in x], after=[dec(v) for v in y], final_difference=dec(y[0]-y[1]))
                break
            x = y
        else: row.update(crossing_iteration=None,escape_iteration=escape)
        rows.append(row)
        (args.out/'partial.json').write_text(json.dumps(dict(protocol,cases=rows),indent=2)+'\n',encoding='utf-8')
        print(json.dumps({key:value for key,value in row.items() if key not in ['epsilon','before','after','final_difference']}),flush=True)
    valid = [r for r in rows if r['crossing_iteration'] is not None]
    best = min(valid,key=lambda r:r['steps_after_escape']) if valid else None
    (args.out/'report.json').write_text(json.dumps(dict(protocol,cases=rows,selected=best,seconds=time.time()-begin),indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mpmath-path',type=Path,required=True);p.add_argument('--phases',type=int,default=32);main(p.parse_args())
