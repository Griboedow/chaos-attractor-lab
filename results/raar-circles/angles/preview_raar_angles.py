"""Compare rigid cameras on retained RAAR gallery clouds without resimulation."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render_twist import camera,density,colorize


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    records=json.loads((args.source/'report.json').read_text())['render'];outputs=[]
    angles=[(.6,.4),(1.4,.4),(2.2,.4),(3.0,.4),(.6,1.),(1.4,1.),(2.2,1.),(3.,1.)]
    for ident in args.ids:
        path=args.source/f'case-{ident}.npz'
        with np.load(path) as d:z=d['states']
        row=next(v for v in records if v['id']==ident);q=(z-np.array(row['center']))@np.array(row['basis'])
        half=np.linalg.norm(q,axis=1).max()*1.04
        fig,axes=plt.subplots(2,4,figsize=(16,9),facecolor='#071320');counts=[]
        for ax,(az,el) in zip(axes.ravel(),angles):
            hist=density(q,camera(az,el),650,half);assert int(hist.sum())==len(z);counts.append(int(hist.sum()))
            im,_=colorize(hist);ax.imshow(im,origin='upper');ax.set_axis_off();ax.set_title(f'azimuth={az:.1f}, elevation={el:.1f}',color='#dfeaf2',fontsize=10)
        fig.text(.025,.955,f'RAAR / R{ident} / RIGID CAMERA STUDY',color='#e7eef6',fontsize=20,weight='bold')
        fig.text(.025,.025,'Same 200,000 actual states in every panel. Equal spatial units; no clipping, deformation or connecting lines.',color='#a8bdce',fontsize=11)
        fig.subplots_adjust(left=.01,right=.99,top=.91,bottom=.06,hspace=.12,wspace=.02)
        fig.savefig(args.out/f'case-{ident}.png',dpi=130,facecolor='#071320');plt.close(fig)
        outputs.append(dict(id=ident,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),angles=angles,point_counts=counts))
    (args.out/'report.json').write_text(json.dumps(dict(source=str(args.source),runs=outputs),indent=2))
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())
