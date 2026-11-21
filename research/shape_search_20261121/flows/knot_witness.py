"""Build or replay exact integer-polygon knot witnesses."""
import argparse
import json
from pathlib import Path

import numpy as np
from cycle_scan import resample
from exact_knot import exact_diagram,prepare


def check(cert):
    for view in cert['views']:
        got=exact_diagram(cert['integer_vertices'], view['integer_matrix'])
        if got['coefficients']!=view['coefficients'] or got['crossings']!=view['crossings']:
            raise AssertionError('Stored polygon witness failed')
    return True


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path)
    ap.add_argument('--out',type=Path)
    ap.add_argument('--verify',type=Path)
    ap.add_argument('--ids',nargs='+',type=int)
    args=ap.parse_args()
    if args.verify:
        for file in sorted(args.verify.glob('case-*.json')):
            check(json.loads(file.read_text()))
        print('All exact polygon witnesses verified.')
        return
    args.out.mkdir(parents=True,exist_ok=False)
    inputs=args.ids or sorted(int(p.stem.split('-')[1]) for p in args.source.glob('case-*.npz'))
    for case in inputs:
        info=json.loads((args.source/f'case-{case}.json').read_text())
        points=np.load(args.source/f'case-{case}.npz')['points']
        arc=resample(points,640)
        cameras=sorted(set(v['camera_seed'] for v in info['views']))
        views=[]
        for camera in cameras:
            vertices,mat,shift,scale=prepare(arc,camera)
            try:result=exact_diagram(vertices,mat)
            except ValueError:continue
            if result['coefficients']!=info['alexander_coefficients']:
                raise ValueError(f'Polygon discrepancy in case {case}')
            result.pop('minimum_depth',None)
            views.append(dict(integer_matrix=mat.tolist(),**result))
        if len(views)<2:
            continue
        cert=dict(case=case,quantization_per_rms=10**9,center=shift.tolist(),rms=float(scale),
                  integer_vertices=vertices.tolist(),views=views,polynomial=info['alexander_coefficients'])
        path=args.out/f'case-{case}.json'
        path.write_text(json.dumps(cert,separators=(',',':'))+'\n')
        check(json.loads(path.read_text()))
        print(case,cert['polynomial'],flush=True)


if __name__=='__main__':
    main()
