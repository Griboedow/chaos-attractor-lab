"""Lay out the two expanding prior cases and the two selected twist geometries."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render_twist import colorize


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    meta=json.loads((args.source/'analysis.json').read_text());selected=[1,2,3,4]
    bg='#071320';fg='#e8f0f5';muted='#9ab0c1';plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg})
    fig,axes=plt.subplots(2,2,figsize=(12,10.7),facecolor=bg)
    for ident,ax in zip(selected,axes.ravel()):
        row=meta['runs'][ident];h=np.load(args.source/f'panel-{ident}.npz')['histogram'];im,_=colorize(h,row['log_density_vmax'])
        ax.imshow(im);ax.set_axis_off();ax.set_title(row['label'],color=fg,fontsize=14,pad=9)
        ax.text(.025,.035,'λ: '+', '.join(f'{v:.3f}' for v in row['spectrum']),transform=ax.transAxes,color=muted,fontsize=10)
    fig.text(.04,.96,'ONE FAMILY, DISTINCT GEOMETRIES',fontsize=22,weight='bold')
    fig.text(.04,.922,'Top: known pumped-top parameter cases. Bottom: selected three-axis twist cases.',fontsize=12,color=muted)
    fig.text(.04,.075,'One million actual states per panel. Rigid orthographic views; each panel fits its own extent and density scale.',fontsize=10,color=muted)
    fig.text(.04,.047,'The 1988 map is an exact affine slice of this family. Different pictures do not establish a new dynamical class.',fontsize=10,color=muted)
    fig.text(.04,.02,'A third OCR-caption parameter case settled to fixed points; all outcomes and that reproduction gap are retained.',fontsize=9,color=muted)
    fig.subplots_adjust(left=.025,right=.975,top=.87,bottom=.11,hspace=.18,wspace=.025)
    fig.savefig(out/'twist-prior-comparison.png',dpi=170,facecolor=bg);plt.close(fig)
    (out/'render.json').write_text(json.dumps(dict(source=str(args.source),analysis_sha256=hashlib.sha256((args.source/'analysis.json').read_bytes()).hexdigest(),
        selected_panels=selected,selection='Both prior parameter cases with positive observed expansion and the two main visual candidates.',
        omitted_panel0='The OCR-caption A=-5,b=.5 case reached a numerical fixed point; documented separately with all33fresh-start outcomes.',
        omitted_panel5='Thin auxiliary candidate307, still included in the complete six-case comparison packet.',
        histogram_sha256={str(i):hashlib.sha256((args.source/f'panel-{i}.npz').read_bytes()).hexdigest() for i in selected}),indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
