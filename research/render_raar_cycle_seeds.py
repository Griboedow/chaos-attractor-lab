"""Show actual tangent-seed evolution and all declared radius/resolution controls."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont
from PIL import Image,ImageDraw,ImageFont
from render_twist import colorize


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();report=json.loads((args.source/'report.json').read_text());histmeta=report['histogram'];sources=[]
    histogram_path=args.histograms or args.source/'evolution-histograms.npz'
    with np.load(histogram_path) as d:hist=d['histograms'].copy();returns=d['returns'].copy();clipped=d['clipped'].copy()
    assert not clipped.any();density_area_factor=hist.shape[1]*hist.shape[2]/(700*900);vmax=float(np.log1p(80));font=ImageFont.truetype(findfont('DejaVu Sans'),22);small=ImageFont.truetype(findfont('DejaVu Sans'),17);title=ImageFont.truetype(findfont('DejaVu Sans'),30)
    with np.load(args.cycles) as d:cycle=d['base-unstable-period-7'];tree=cKDTree(cycle)
    yy,xx=np.nonzero(hist.any(0));crop=[max(0,int(xx.min())-15),max(0,int(yy.min())-15),min(hist.shape[2],int(xx.max())+16),min(hist.shape[1],int(yy.max())+16)];left,top,right,bottom=crop
    projection=np.array(histmeta['projection']);center=np.array(histmeta['camera_center']);uv=(cycle-center)@projection;xb=histmeta['x_bounds'];yb=histmeta['y_bounds'];cycle_pixels=np.c_[(uv[:,0]-xb[0])/(xb[1]-xb[0])*hist.shape[2]-left,(yb[1]-uv[:,1])/(yb[1]-yb[0])*hist.shape[1]-top]
    scale=min(1.,1000/(right-left),740/(bottom-top));size=(round((right-left)*scale),round((bottom-top)*scale));origin=((1100-size[0])//2,105+(700-size[1])//2);frames=[]
    for h,k in zip(hist,returns):
        canvas=Image.new('RGB',(1100,900),'#071320');draw=ImageDraw.Draw(canvas);draw.text((45,20),'A TINY SEGMENT UNFOLDS INTO THE ARCH',font=title,fill='#e6edf5');draw.text((45,66),f'Actual iterations {7*k:03d}–{7*k+6:03d} / {k} seven-step returns',font=font,fill='#f5c181')
        im,_=colorize(h[top:bottom,left:right]*density_area_factor,vmax);canvas.paste(im.resize(size,Image.Resampling.NEAREST),origin)
        if k<=12:
            for u,v in cycle_pixels:
                a=origin[0]+u*scale;b=origin[1]+v*scale;draw.ellipse((a-5,b-5,a+5,b+5),outline='#73b9ce',width=1)
        draw.text((45,822),'16,384 tangent seeds / radius 10⁻⁷ / fixed rigid camera / equal spatial units',font=small,fill='#a7bbcc');draw.text((45,846),'Early blue rings locate the tiny initial clusters; they are annotations.',font=small,fill='#a7bbcc');draw.text((45,870),'Finite evolution; no certified invariant manifold or historical-priority claim.',font=small,fill='#a7bbcc');frames.append(canvas)
        if k in [0,12,16,20,32,64]:canvas.save(args.out/f'evolution-frame-{k:02}.png')
    duration=[130]*65;duration[0]=1000;duration[-1]=2000;frames[0].save(args.out/'raar-cycle-unfolding.gif',save_all=True,append_images=frames[1:],duration=duration,loop=0,optimize=False)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(12,7.2),facecolor='#071320')
    for ax,k in zip(axes.flat,[0,12,16,20,32,64]):
        im,_=colorize(hist[k,top:bottom,left:right]*density_area_factor,vmax);ax.imshow(im,interpolation='nearest');ax.set_title(f'Iterations {7*k}–{7*k+6}',color='#e6edf5',fontsize=11);ax.set_axis_off()
        if k<=12:ax.scatter(cycle_pixels[:,0],cycle_pixels[:,1],s=35,facecolors='none',edgecolors='#73b9ce',linewidths=.6)
    fig.suptitle('Finite evolution of a tiny segment near the unstable seven-cycle',color='#e6edf5',fontsize=16);fig.subplots_adjust(left=.015,right=.985,bottom=.02,top=.90,wspace=.04,hspace=.12);fig.savefig(args.out/'raar-cycle-unfolding.png',dpi=220,facecolor=fig.get_facecolor());fig.savefig(args.out/'raar-cycle-unfolding.pdf',dpi=300,facecolor=fig.get_facecolor());plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,3.6));records=[];colors=['#ab5b0a','#087f9c','#7956aa']
    for row in report['runs']:
        path=args.source/(row['key']+'.npz');sources.append(dict(path=str(path),sha256=sha(path)))
        with np.load(path) as d:
            k=d['returns'];z=d['snapshots'];distance=d['reference_distances'];spread=np.array([np.median(tree.query(v.reshape(-1,3))[0]) for v in z]);med=np.median(distance,axis=1)
        color=colors[report['radii'].index(row['radius'])];style='-' if row['seed_count']==16384 else '--';label=f'ε={row["radius"]:.0e}, n={row["seed_count"]:,}'
        axes[0].semilogy(7*k,spread,style,c=color,lw=1.2,label=label);axes[1].semilogy(7*k,med,style,c=color,lw=1.2)
        records.append(dict(key=row['key'],iterations=(7*k).tolist(),median_distance_to_cycle=spread.tolist(),median_distance_to_reference=med.tolist()))
    axes[0].set(xlabel='Iterations (first phase of seven-state window)',ylabel='Median distance to the seven-cycle',title='The segment leaves the unstable cycle');axes[0].legend(fontsize=7,frameon=False,ncol=2)
    floor=report['runs'][0]['reference_split_quantiles'][0];axes[1].axhline(floor,c='#555555',lw=.8,ls=':',label='Reference split sampling median');axes[1].set(xlabel='Iterations (first phase of seven-state window)',ylabel='Median distance to finite reference cloud',title='Later states approach its sampling scale');axes[1].legend(fontsize=7,frameon=False)
    fig.subplots_adjust(left=.08,right=.98,bottom=.18,top=.87,wspace=.28);fig.savefig(args.out/'raar-cycle-seed-controls.png',dpi=250);fig.savefig(args.out/'raar-cycle-seed-controls.pdf');plt.close(fig)
    for name in ['render_raar_cycle_seeds.py','render_twist.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    assert np.array_equal(hist.sum((1,2)),hist[:,top:bottom,left:right].sum((1,2)))
    analysis=dict(source=str(args.source),source_report_sha256=sha(args.source/'report.json'),histogram_path=str(histogram_path),histogram_sha256=sha(histogram_path),cycle_source_sha256=sha(args.cycles),sources=sources,controls=records,frames=65,frame_point_count=7*16384,log_density_vmax=vmax,fixed_pixel_crop=crop,uniform_image_scale=scale,early_locator_rings='Annotations at certified cycle positions through return 12; not trajectory samples.',clipped=int(clipped.sum()),scope='Images contain every stored phase-window point. Camera and density scale are fixed. A common empty-margin crop and uniform nearest-neighbor image scaling preserve aspect and retain every point. Discrete points are never connected. Tangent seeds are not certified unstable-manifold coordinates. Reference distances compare finite samples of unequal densities and do not prove asymptotic attraction or set equality.')
    analysis['density_area_factor']=density_area_factor
    analysis['scope']='Every phase-window point contributes to the histogram. No histogram counts are clipped by the common crop. Camera, density per unit area and spatial scale are fixed. Discrete points are not connected. Tangent seeds are not certified unstable-manifold coordinates. Reference distances compare finite samples of unequal densities and do not prove asymptotic attraction or set equality.'
    (args.out/'analysis.json').write_text(json.dumps(analysis,indent=2),encoding='utf-8');print(json.dumps({k:analysis[k] for k in ['frames','frame_point_count','clipped']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--cycles',type=Path,default=Path('results/raar-cycle-structure/figures/cycle-coordinates.npz'));p.add_argument('--histograms',type=Path);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
