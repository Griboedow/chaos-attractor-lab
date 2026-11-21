"""Verify linking directly on the stored float64 pair polylines as exact dyadic rationals."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from exact_knot import prepare
from symmetry_links import link_integer


def coordinates(arrays):
    parts=[[[float(x).as_integer_ratio() for x in row] for row in z] for z in arrays]
    denominator=max(d for rows in parts for row in rows for _,d in row)
    integers=[[tuple(n*(denominator//d) for n,d in row) for row in rows] for rows in parts]
    return integers,denominator


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pairs',type=Path,required=True)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();rows=json.loads(a.pairs.read_text());out=[]
    for row in rows:
        path=a.cache/row['cloud_relative'];blob=hashlib.sha256(path.read_bytes()).hexdigest()
        if blob!=row['sha256']:raise ValueError('Cloud hash mismatch')
        with np.load(path) as z: first,second=np.array_split(z['points'],2)
        ints,den=coordinates([first,second]);views=[]
        for camera in (17,35,81):
            _,matrix,_,_=prepare(first[::64],camera)
            value,cross=link_integer(*ints,matrix)
            views.append(dict(camera=matrix.tolist(),linking=value,crossings=cross))
        if len({v['linking'] for v in views})!=1:raise ValueError('Exact projections disagree')
        out.append(dict(order=row['child_period_order'],source_sha256=blob,cloud_relative=row['cloud_relative'],
                   vertices_each=len(first),exact_coordinate_denominator=str(den),views=views))
        a.out.write_text(json.dumps(out,indent=2)+'\n')
        print(row['child_period_order'],views[0]['linking'],flush=True)


if __name__=='__main__':
    main()
