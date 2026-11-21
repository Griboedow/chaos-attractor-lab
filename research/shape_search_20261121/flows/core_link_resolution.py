"""Resolve mesh-sensitive linking between the final order-two and order-four cores."""
import argparse
import json
from pathlib import Path
import numpy as np
from dyadic_links import coordinates
from exact_knot import prepare
from relative_cores import sector_cloud
from replay_case import orbit
from symmetry_links import link_integer


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,required=True);ap.add_argument('--sizes',type=int,nargs='+',default=[65536,131072])
    a=ap.parse_args();r=a.packet;rows={x['period_order']:x for x in json.loads((r/'nested-skeletons.json').read_text())['cores']}
    p=np.array(rows[2]['parameters']);out=[]
    for n in a.sizes:
        first=orbit(rows[2]['initial'],rows[2]['period'],p,1,n,3e-12)
        second,_=sector_cloud(np.array(rows[4]['initial']),rows[4]['sector_period'],p,rows[4]['relative_rotation_index'],n)
        polys,den=coordinates([first,second]);views=[]
        for seed in (17,35):
            _,m,_,_=prepare(first[::64],seed);lk,cross=link_integer(*polys,m)
            views.append(dict(matrix=m.tolist(),linking=lk,crossings=cross))
        out.append(dict(count=n,exact_dyadic_denominator=str(den),views=views))
        (r/'nested-link-refinement.json').write_text(json.dumps(out,indent=2)+'\n')
        print(n,views[0]['linking'],flush=True)


if __name__=='__main__':
    main()
