"""Cycle-structure figures: actual coordinates, frozen outcomes, fixed rigid camera."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from continue_raar_cycles import iterate
from render_twist import camera
from verify_rk23_controller import sha,read


def figures(root,base,out):
    out.mkdir(parents=True,exist_ok=True);sources=[]
    def load(path):
        sources.append(dict(path=str(path),sha256=sha(path)))
        return json.loads(path.read_text()) if path.suffix=='.json' else read(path)
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=300,facecolor='white')
        fig.savefig(out/(name+'.pdf'),dpi=300,facecolor='white');plt.close(fig)
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False})
    scan=load(root/'fine-slice/scan.npz');events=load(root/'multiplier-proof/report.json')['cases'];cycles=load(root/'cycle-proof/report.json')['cases']
    beta=scan['params'][:,0];tails=scan['tails'];palette=['#087f9c','#b64d71','#bd7a14','#6b50a0']
    fig,axes=plt.subplots(3,1,figsize=(7,6.1),gridspec_kw={'height_ratios':[1.6,1.6,1]})
    for ax,limits,title in [(axes[0],(.9608,.9664),'A  Fresh-start outcomes: seven bands split'),(axes[1],(.9646,.9664),'B  Lower band, enlarged')]:
        mask=(beta>=limits[0])&(beta<=limits[1]);x=np.repeat(beta[mask],tails.shape[1]);y=tails[mask,:,0].ravel()
        ax.scatter(x,y,s=.25,color='#075d78',alpha=.6,rasterized=True,linewidths=0)
        ax.set(xlim=limits,ylabel='Actual state coordinate x',title=title)
        if ax==axes[1]:ax.set_ylim(-.43,-.375)
        for j,e in enumerate(events):
            b=float(e['beta_display'])
            if limits[0]<b<limits[1]:ax.axvline(b,c=palette[j],lw=.7,ls='--',label=f'p = {e["period"]}')
        ax.legend(loc='center' if ax==axes[0] else 'lower left',frameon=False,ncol=4,fontsize=6)
    rates=scan['rates'][:,:,0];axes[2].fill_between(beta,rates.min(1),rates.max(1),color='#d6e9ed',lw=0,label='20-window range');axes[2].plot(beta,rates.mean(1),c='#087f9c',lw=.7,label='Mean')
    axes[2].axhline(0,c='#444444',lw=.6);axes[2].set(xlim=(.9608,.9664),xlabel='Relaxation β',ylabel='Finite growth / step',title='C  Full-state leading growth estimate');axes[2].legend(frameon=False,fontsize=6)
    fig.subplots_adjust(left=.11,right=.985,bottom=.08,top=.96,hspace=.43);save(fig,'figure-1-bifurcation')
    fig,axes=plt.subplots(1,4,figsize=(7,2.2))
    for ax,e,col in zip(axes,events,palette):
        period=e['period'];d=load(root/f'continuation/period-{period}-branch.npz');b0=float(e['beta_display']);width={7:.0025,14:.00075,28:.00018,56:.00004}[period]
        mask=np.abs(d['beta']-b0)<width
        for j in range(3):ax.plot((d['beta'][mask]-b0)*1e6,d['multipliers'][mask,j].real,lw=.8,c=col)
        ax.axhline(-1,c='#333333',ls='--',lw=.6);ax.axvline(0,c='#999999',lw=.6)
        ax.scatter([0],[-1],s=20,c=col,zorder=3);ax.set(title=f'Period {period}',xlabel='β offset × 10⁶',ylim=(-1.5,1.15))
        if ax==axes[0]:ax.set_ylabel('Numerical Floquet multipliers')
    fig.subplots_adjust(left=.085,right=.99,bottom=.24,top=.83,wspace=.36);save(fig,'figure-2-multipliers')
    render=load(base/'render/analysis.json')['runs'][0];cloud=load(base/'render/case-1051-animation.npz')['states'];basis=np.array(render['basis']);center=np.array(render['center']);cam=camera(render['azimuth'],render['elevation']);xy=(cloud-center)@basis@cam.T
    cycle_data=[]
    for row in cycles:
        p=np.array(row['parameters']);x=np.array([float(v) for v in row['center_decimal']]);end,M,axis,orbit,status=iterate(x,p,row['period']);assert status==0
        cycle_data.append(dict(label=row['label'],period=row['period'],beta=p[0],orbit=orbit,projection=(orbit-center)@basis@cam.T,multipliers=np.linalg.eigvals(M)))
    unstable=next(r for r in cycle_data if r['label']=='base-unstable-period-7');q=unstable['projection'];dist=cKDTree(cloud).query(unstable['orbit'])[0]
    allxy=np.concatenate([xy]+[r['projection'] for r in cycle_data]);lo=allxy.min(0);hi=allxy.max(0);mid=(lo+hi)/2;half=max(hi-lo)*.55
    fig=plt.figure(figsize=(7,5.5));gs=fig.add_gridspec(2,4,height_ratios=[2.3,1],hspace=.12,wspace=.02);ax=fig.add_subplot(gs[0,:])
    ax.scatter(xy[:,0],xy[:,1],s=.04,c='#087c96',alpha=.65,linewidths=0,rasterized=True)
    ax.scatter(q[:,0],q[:,1],s=26,c='#df8628',edgecolors='white',linewidths=.5,zorder=4)
    for k,point in enumerate(q):ax.annotate(str(k+1),point,xytext=(5,5),textcoords='offset points',fontsize=8,color='#8b4c0c')
    ax.set_title('The base arch and its certified unstable 7-cycle',pad=7)
    ax.set(xlim=(xy[:,0].min()-.1,xy[:,0].max()+.1),ylim=(xy[:,1].min()-.12,xy[:,1].max()+.12),aspect='equal');ax.set_axis_off()
    for j,(period,col) in enumerate(zip([7,14,28,56],palette)):
        r=next(r for r in cycle_data if r['label']==f'stable-period-{period}');ax=fig.add_subplot(gs[1,j]);z=r['projection'];ax.scatter(z[:,0],z[:,1],s=8,c=col,linewidths=0)
        ax.set(xlim=(mid[0]-half,mid[0]+half),ylim=(mid[1]-half,mid[1]+half),aspect='equal',title=f'{period} states\nβ = {r["beta"]:.6f}');ax.set_axis_off()
    fig.subplots_adjust(left=.01,right=.99,bottom=.025,top=.94);save(fig,'figure-3-cycle-geometry')
    np.savez_compressed(out/'cycle-coordinates.npz',**{r['label']:r['orbit'] for r in cycle_data},base_cycle_reference_distances=dist)
    report=dict(sources=sources,camera=dict(center=center.tolist(),basis=basis.tolist(),orthographic=cam.tolist()),unstable_cycle_multipliers=[[float(v.real),float(v.imag)] for v in unstable['multipliers']],unstable_cycle_reference_distances=dist.tolist(),scope='Cycle dots are actual discrete states, not a continuous closed trajectory. The million-point reference is finite; nonzero cycle distances do not prove attractor membership. Insets share camera and spatial scale; the main panel uses the same camera with a closer extent. Lines in multiplier plots join numerical continuation samples; vertical event marks identify certified tiny parameter brackets, not a proved nondegenerate flip.')
    (out/'sources.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({k:report[k] for k in ['unstable_cycle_multipliers','unstable_cycle_reference_distances']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--base',type=Path,default=Path('results/raar-circles'));p.add_argument('--out',type=Path,required=True);a=p.parse_args();figures(a.root,a.base,a.out)
