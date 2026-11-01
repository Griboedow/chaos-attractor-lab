"""Show the sampled core and every recorded threshold exceedance at equal scale."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def main(a):
    a.out.mkdir(parents=True,exist_ok=False)
    r=json.loads((a.source/'report.json').read_text(encoding='utf-8'))
    camera=json.loads(a.camera.read_text(encoding='utf-8'))
    with np.load(a.source/'uniform-sample.npz') as f:z=f['states']
    with np.load(a.source/'all-threshold-exceedances.npz') as f:events=f['states'];times=f['indices']
    basis=np.array(camera['basis']);projection=np.array(camera['projection']);center=np.array(camera['center'])
    assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12
    assert np.max(abs(projection@projection.T-np.eye(2)))<1e-12
    xy=(z-center)@basis@projection.T;ev=(events-center)@basis@projection.T
    half=float(max(abs(xy).max(),abs(ev).max())*1.08)
    hist=np.histogram2d(xy[:,0],xy[:,1],bins=1300,range=[[-half,half]]*2)[0].T
    value=np.log1p(hist);cap=float(np.percentile(value[value>0],99.5))
    ink,muted,blue,orange='#233f4b','#71848b','#2b7c9a','#c66c31'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':ink,'axes.labelcolor':ink,'axes.edgecolor':'#cbd4d7','xtick.color':muted,'ytick.color':muted})
    fig=plt.figure(figsize=(14,7.8),facecolor='white')
    grid=fig.add_gridspec(2,2,left=.045,right=.965,bottom=.19,top=.79,width_ratios=[1.08,1],hspace=.57,wspace=.26)
    ax=fig.add_subplot(grid[:,0]);bx=fig.add_subplot(grid[0,1]);cx=fig.add_subplot(grid[1,1])
    cmap=LinearSegmentedColormap.from_list('posterior',['#ffffff','#c7dfe5','#559ab1','#244d63'])
    ax.imshow(value,origin='lower',extent=[-half,half,-half,half],vmin=0,vmax=cap,cmap=cmap,interpolation='nearest')
    ax.scatter(ev[:,0],ev[:,1],s=22,facecolors='none',edgecolors=orange,linewidths=1.1,label='All 56 threshold exceedances')
    ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal',xlabel='Orthogonal camera coordinate 1',ylabel='Orthogonal camera coordinate 2')
    ax.set_title('A   The dense wings and rare exits',loc='left',fontsize=13,pad=18)
    ax.legend(loc='lower left',frameon=False,fontsize=9)
    lower=np.array(r['threshold_lower']);upper=np.array(r['threshold_upper'])
    excess=np.maximum(np.max(events-upper,axis=1),np.max(lower-events,axis=1))
    bx.vlines(times/1e6,0,excess,color=orange,linewidth=1.2)
    bx.scatter(times/1e6,excess,s=10,color=orange,zorder=3)
    bx.set(xlim=(0,100),ylim=(0,.51),xlabel='Observed iterations (millions)',ylabel='Excess beyond envelope')
    bx.set_title('B   56 exits in 100 million updates',loc='left',fontsize=13,pad=18)
    minima=np.array([q['minimum_radius'] for q in r['rows']])
    cx.plot((np.arange(50)+.5)*2,minima,'o',color=blue,markersize=3.7)
    cx.step(np.arange(50)*2,np.minimum.accumulate(minima),where='post',color=orange,linewidth=1.5)
    cx.set(xlim=(0,100),ylim=(0,1.65),xlabel='Observed iterations (millions)',ylabel=r'Minimum $\|2m(x)-x\|$')
    cx.set_title('C   Smaller normalization denominators',loc='left',fontsize=13,pad=18)
    cx.text(.98,.08,'Dots: two-million-step blocks\nLine: cumulative minimum',transform=cx.transAxes,ha='right',fontsize=9,color=muted)
    for axis in [ax,bx,cx]:
        axis.spines[['top','right']].set_visible(False)
    for axis in [bx,cx]:axis.grid(axis='y',color='#e8edef',linewidth=.7);axis.set_axisbelow(True)
    fig.text(.045,.942,'Beyond the four-wing portrait',fontsize=23)
    fig.text(.045,.884,'Product alphabet · Case 988 · the visible core is not a forward-invariant envelope',fontsize=12,color=muted)
    fig.text(.045,.105,'Blue density: every 50th state over the full run. Orange rings: every exit, shown separately without joining lines.',fontsize=10,color=muted)
    fig.text(.045,.060,'The threshold is the original two-million-state coordinate box enlarged by 0.02. Three short escapes are interval-certified.',fontsize=10,color=muted)
    fig.savefig(a.image,dpi=180,facecolor='white');plt.close(fig)
    root=Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',a.image.resolve().relative_to(root).as_posix(),'--role','diagnostic'],check=True)
    np.savez_compressed(a.out/'figure-data.npz',histogram=hist,events=ev,event_times=times,event_excess=excess,block_minima=minima,halfwidth=half,log_saturation=cap,center=center,basis=basis,projection=projection)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (a.out/'figure.png').write_bytes(a.image.read_bytes());(a.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (a.out/'camera.json').write_bytes(a.camera.read_bytes())
    result=dict(status='rendered',inputs={str(p):sha(p) for p in [a.source/'report.json',a.source/'uniform-sample.npz',a.source/'all-threshold-exceedances.npz',a.camera]},image_sha256=sha(a.image),histogram_count=int(hist.sum()),event_count=len(events),halfwidth=half,scope='Equal-unit orthogonal projection with no clipping of the plotted sample or events; log-density saturation declared. Uniformly thinned sample is explicitly identified; all threshold exceedances are overlaid separately. No continuous trajectories, invariant support, basin probability or asymptotic classification is inferred.')
    (a.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','camera','out','image']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
