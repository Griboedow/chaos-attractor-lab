"""Compute exact pairwise linking of successive unstable cores at the final chaotic parameters."""
import argparse
import json
from pathlib import Path
import numpy as np
from dyadic_links import coordinates
from exact_knot import prepare
from symmetry_links import link_integer


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True)
    a=ap.parse_args();r=a.packet
    base=next(row for row in json.loads((r/'manifest.json').read_text())['knotted_periodic_skeletons'] if row['case']==180844)
    data=json.loads((r/'nested-skeletons.json').read_text())['cores'];previous=base;out=[]
    for row in data:
        cloud=[np.load(a.cache/x['cloud_relative'])['points'] for x in (previous,row)]
        poly,den=coordinates(cloud);views=[]
        for seed in (17,35,81):
            _,mat,_,_=prepare(cloud[0][::64],seed)
            num,cross=link_integer(poly[0],poly[1],mat)
            views.append(dict(matrix=mat.tolist(),linking=num,crossings=cross))
        assert len({v['linking'] for v in views})==1
        out.append(dict(parent=previous.get('period_order',1),child=row['period_order'],pair_shas=[previous['sha256'],row['sha256']],
                     dyadic_denominator=str(den),vertex_counts=[len(z) for z in cloud],views=views))
        (r/'nested-core-links.json').write_text(json.dumps(out,indent=2)+'\n')
        print(row['period_order'],views[0]['linking'],flush=True);previous=row


if __name__=='__main__':
    main()
