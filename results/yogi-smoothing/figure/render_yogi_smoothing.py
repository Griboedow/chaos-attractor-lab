"""Show every declared relative-smoothing control, including regular outcomes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from check_yogi_smoothing import smooth
from yogi_map import step
from render_twist import camera,density,colorize


@njit
def cloud(p,start,h):
    z=start.copy();out=np.empty((200000,3))
    for i in range(200000):
        if h==0:z=step(z,p)
        else:z,_,_=smooth(z,p,h)
        out[i]=z
    return out


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    analysis=json.loads((args.source/'analysis.json').read_text());rows=analysis['runs']
    with np.load(args.source/'histories.npz') as a:d={k:a[k] for k in a.files}
    fig=plt.figure(figsize=(19,8.5),facecolor='#071320');grid=fig.add_gridspec(2,6,width_ratios=[1.2,1.2,1,1,1,1])
    fg='#e8f0f5';muted='#9ab0c1';colors=['#ffbd69','#65d4f5','#a58de0'];records=[];histories=[];indices=[]
    for row,ident in enumerate([670,125]):
        ax=fig.add_subplot(grid[row,:2]);ax.set_facecolor('#071320')
        for kind,style,alpha in [('cold','-',1.),('original_late_state','--',.65)]:
            chosen=[i for i,r in enumerate(rows) if r['grid_id']==ident and r['start_kind']==kind]
            h=d['widths'][chosen];w=d['windows'][chosen]
            for j,color in enumerate(colors):
                ax.plot(h,w[:,:,j].mean(1),style,color=color,alpha=alpha,lw=1.4,label=f'rate {j+1}' if kind=='cold' else None)
                if kind=='cold':ax.fill_between(h,w[:,:,j].min(1),w[:,:,j].max(1),color=color,alpha=.15)
        ax.set_xscale('symlog',linthresh=1e-6,linscale=.3);ax.axhline(0,color=muted,lw=.8);ax.set_ylim(-.24,.27)
        ax.set_xlabel('Relative switch width h (0 = original sign)',color=muted,fontsize=10)
        ax.set_title(f'G{ident:04d} / finite growth rates',color=fg,fontsize=13);ax.tick_params(colors=muted,labelsize=8)
        for spine in ax.spines.values():spine.set_color('#2b4358')
        if row==0:ax.legend(frameon=False,labelcolor=fg,fontsize=9,loc='upper left')
        ax.text(.01,.025,'Solid: cold start; dashed: original late state.\nShading: 20 window minima/maxima. Third rate can leave view.',transform=ax.transAxes,color=muted,fontsize=8)
        for column,target in enumerate([0.,.01,.03162277660168379,1.]):
            i=next(i for i,r in enumerate(rows) if r['grid_id']==ident and r['start_kind']=='cold' and abs(r['width']-target)<1e-12)
            z=cloud(d['parameters'][i],d['ends'][i],d['widths'][i]);histories.append(z);indices.append(i)
            center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(z-center)@basis;shift=(q.min(0)+q.max(0))/2;q-=shift;center+=basis@shift
            half=max(np.linalg.norm(q,axis=1).max()*1.04,1e-12);view=camera(.4);hist=density(q,view,550,half)
            assert hist.sum()==len(z);picture,vmax=colorize(hist)
            panel=fig.add_subplot(grid[row,column+2]);panel.imshow(picture,origin='lower');panel.set_axis_off()
            period=rows[i]['period_up_to128'];label=f'period {period} detected' if period else 'no return up to 128'
            panel.set_title(f'h={target:.3g}\n{label}',color=fg,fontsize=10)
            records.append(dict(row=i,grid_id=ident,width=target,center=center.tolist(),basis=basis.tolist(),projection=view.tolist(),halfwidth=float(half),points=len(z),clipped=0))
    fig.text(.035,.96,'YOGI / A DECLARED SMOOTHING CONTROL',color=fg,fontsize=22,weight='bold')
    fig.text(.035,.91,'The relative tanh switch changes the equation. This control is distinct from the pre-existing fixed-scale tanh in TensorFlow Addons.',color=muted,fontsize=10)
    fig.text(.035,.035,'56 fresh runs: 14 widths x 2 starts x 2 parameter sets. Each has 100,000 warmup + 1,000,000 observed steps. Portraits: 200,000 further states.',color=muted,fontsize=10)
    fig.subplots_adjust(left=.045,right=.98,top=.83,bottom=.12,wspace=.28,hspace=.45)
    fig.savefig(out/'yogi-smoothing-control.png',dpi=180,facecolor=fig.get_facecolor());plt.close(fig)
    np.savez_compressed(out/'clouds.npz',source_rows=indices,states=histories)
    (out/'render.json').write_text(json.dumps(dict(source=str(args.source),source_sha256=hashlib.sha256((args.source/'histories.npz').read_bytes()).hexdigest(),portraits=records),indent=2)+'\n')
    for n in ['render_yogi_smoothing.py','check_yogi_smoothing.py','yogi_map.py','render_twist.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    for ident in [670,125]:
        print('ID',ident)
        for r in rows:
            if r['grid_id']==ident and r['start_kind']=='cold':print(r['width'],r['spectrum'],r['window_minima'],r['period_up_to128'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
