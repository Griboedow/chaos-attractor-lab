"""Parameter slices and independent-library cloud comparison, without deformations."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from render_twist import camera,colorize
from render_yogi_arc import histogram


def read(path):
    with np.load(path) as f:return {k:f[k] for k in f.files}


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);base=args.external
    d=read(base/'raar-circles-relaxation01/scan.npz');report=json.loads((base/'raar-circles-relaxation01/report.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(3,1,figsize=(12,9),sharey=True)
    for row,(ident,ax) in enumerate(zip([1051,1483,820],axes)):
        sl=slice(row*257,(row+1)*257);p=d['params'][sl];r=d['rates'][sl];x=p[:,0]
        for k,color in [(0,'#007d92'),(1,'#a64f00')]:
            ax.plot(x,r[:,:,k].mean(1),color=color,lw=.9,label=f'Rate {k+1}: mean and window range')
            ax.fill_between(x,r[:,:,k].min(1),r[:,:,k].max(1),color=color,alpha=.16)
        regular=d['period'][sl]>0
        ax.scatter(x[regular],np.full(regular.sum(),-.20),marker='|',s=50,c='#7d568f',label='Detected period 1–512')
        ax.axhline(0,color='#888888',lw=.7);ax.axvline(x[128],color='#777777',ls=':',lw=.8)
        g=report['groups'][row];ax.set(title=f'R{ident} / {g["positive"]}/257 leading-positive, {g["two_positive"]} two-positive, {g["periodic"]} periodic',xlabel='Relaxation beta',ylabel='Finite growth / step',ylim=(-.22,.4),xlim=(x[0],x[-1]))
        if row==0:ax.legend(loc='upper left',ncol=3,fontsize=8,frameon=False)
    fig.suptitle('Changing one parameter exposes alternating irregular and periodic regimes',fontsize=16)
    fig.text(.06,.015,'Fresh start at every setting; 100,000 warm-up + 400,000 measured steps. Bands are window ranges, not confidence intervals.',fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.96));fig.savefig(args.out/'raar-relaxation.png',dpi=200);fig.savefig(args.out/'raar-relaxation.pdf');plt.close(fig)
    renders=json.loads((base/'raar-circles-render01/analysis.json').read_text())['runs'];comparisons=[]
    fig,axes=plt.subplots(2,2,figsize=(13,9),facecolor='#071320')
    for row,ident in enumerate([1051,1483]):
        render=next(r for r in renders if r['id']==ident);basis=np.array(render['basis']);center=np.array(render['center']);proj=camera(render['azimuth'],render['elevation'])
        ref=read(base/f'raar-circles-gallery01/case-{ident}.npz')['states'];lib=read(base/f'raar-circles-upstream01/case-{ident}.npz')['states'][-200000:]
        allxy=np.concatenate([(z-center)@basis@proj.T for z in [ref,lib]])
        unit=max(2*np.abs(allxy[:,0]).max()/800,2*np.abs(allxy[:,1]).max()/520)*1.03
        h=[histogram((z-center)@basis,proj,800,520,unit) for z in [ref,lib]];vmax=max(float(np.log1p(q).max()) for q in h)
        for col,(q,label) in enumerate(zip(h,['Transcription','Official ProxToolbox evaluator'])):
            assert q.sum()==200000;im,_=colorize(q,vmax=vmax);ax=axes[row,col];ax.imshow(im);ax.set_axis_off();ax.set_title(f'R{ident} / {label}',color='#e6edf5',fontsize=12)
        distance=cKDTree(ref).query(lib,workers=4)[0];baseline=cKDTree(ref[:100000]).query(ref[100000:],workers=4)[0]
        comparisons.append(dict(id=ident,points=200000,units_per_pixel=unit,shared_density_scale=vmax,nearest_reference_median=float(np.median(distance)),nearest_reference_max=float(distance.max()),reference_split_median=float(np.median(baseline))))
        np.savez_compressed(args.out/f'case-{ident}-comparison.npz',histograms=h,distance=distance,baseline=baseline)
    fig.suptitle('The same forms survive original-library execution',color='#e6edf5',fontsize=18)
    fig.text(.035,.018,'Actual states; identical rigid camera, equal units and density scale within each row. Cloud agreement is finite, not exact orbit matching.',color='#b0c4d6',fontsize=9)
    fig.subplots_adjust(top=.91,bottom=.06,left=.02,right=.98,hspace=.13,wspace=.02)
    fig.savefig(args.out/'raar-library-clouds.png',dpi=180,facecolor=fig.get_facecolor());fig.savefig(args.out/'raar-library-clouds.pdf',facecolor=fig.get_facecolor());plt.close(fig)
    (args.out/'report.json').write_text(json.dumps(dict(comparisons=comparisons,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),slice_source_sha256=hashlib.sha256((base/'raar-circles-relaxation01/scan.npz').read_bytes()).hexdigest()),indent=2),encoding='utf-8')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--external',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
