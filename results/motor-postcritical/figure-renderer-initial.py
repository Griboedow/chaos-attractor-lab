"""Render physical cluster traces and finite-time saddle-connection evidence."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from motor_postcritical_escape import physical

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/motor-postcritical'
COLORS=['#0072B2','#D55E00','#009E73','#CC79A7']


def render(output,manifest):
    if manifest.exists():raise FileExistsError(manifest)
    output.mkdir(parents=True,exist_ok=True)
    inputs=[BASE/'paired-cycle-rho-3p3652184-retry/cycle.npz',BASE/'paired-branch/analysis.json',
            BASE/'paired-second-crossing/analysis.json',ROOT/'results/motor-rows/cluster-cycles/analysis.json',
            BASE/'saddle-escape/trajectory.npz',BASE/'saddle-escape/diagnostics.json',BASE/'saddle-escape/summary.json']
    cycle=np.load(inputs[0]);branch=json.loads(inputs[1].read_text());second=json.loads(inputs[2].read_text())
    parent=json.loads(inputs[3].read_text());escape=np.load(inputs[4]);diagnostic=json.loads(inputs[5].read_text());meta=json.loads(inputs[6].read_text())
    plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.labelsize':10,'legend.fontsize':8,
                         'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
    fig,axes=plt.subplots(1,3,figsize=(15.5,4.8),layout='constrained')
    q=physical(cycle['full_states']);period=cycle['time'][-1]
    for k in range(4):
        xx=q[:,2*k,2]
        axes[0].plot(np.r_[cycle['time'],cycle['time'][1:]+period]/period,np.r_[xx,xx[1:]],
                     color=COLORS[k],lw=1.7,label=f'Rows {2*k+1}-{2*k+2}')
    axes[0].set(xlabel='Time / doubled-cycle period',ylabel='Physical row displacement X',
                title='Four synchronized pairs: a stable periodic state')
    axes[0].legend(ncol=2,loc='upper right');axes[0].grid(alpha=.15)
    axes[0].text(.02,.02,'rho = 3.3652184; T = 6.6904415\nEach pair follows a quarter-period shifted trajectory',
                 transform=axes[0].transAxes,fontsize=8,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
    selected=[c for c in branch['cases'] if c['label']=='toward_chaos']
    delta=np.array([c['delta'] for c in selected])
    axes[1].plot(delta,[c['largest_quotient_nonflow_modulus'] for c in selected],color=COLORS[0],label='Paired quotient, excluding flow')
    axes[1].plot(delta,[c['largest_transverse_modulus'] for c in selected],color=COLORS[1],label='Within-pair contrasts')
    pd=np.array([(c['parameters'][1]-1)**2/c['parameters'][0]**2-1 for c in parent['cases']])
    axes[1].plot(pd,[max(v['modulus'] for v in c['transverse_A']['multipliers']) for c in parent['cases']],
                 ':',color='#555555',label='Parent 4+4 contrasts')
    axes[1].axhline(1,color='black',lw=.8)
    axes[1].axvline(.2414738454054,color='#7B4B94',ls='--',lw=1.3,label='First crossing: interval certificate')
    axes[1].axvline(second['delta'],color='#B45F06',ls='-.',lw=1.3,label='Second crossing: numerical')
    axes[1].set(xlim=(.19,.36),ylim=(0,2),xlabel='Source activity excess delta',ylabel='Largest nonflow multiplier modulus',
                title='The doubled state later loses transverse stability')
    axes[1].legend(loc='upper left',fontsize=7.5);axes[1].grid(alpha=.15)
    near=[c for c in branch['cases'] if c['label']=='near_onset']
    eps=np.array([c['offset_below_first_crossing'] for c in near]);amp=np.array([c['split_pair_rms_distance'] for c in near])
    axes[2].loglog(eps,amp,'o',color=COLORS[0],ms=6,label='Continued periodic states')
    x=np.geomspace(eps.min(),eps.max(),100)
    axes[2].loglog(x,np.sqrt(3.6093*x),'--',color=COLORS[1],label='Reference slope 1/2')
    axes[2].set(xlabel='rho at first crossing minus rho',ylabel='Pair-splitting RMS amplitude',
                title='Near-onset splitting follows square-root scaling')
    axes[2].legend(loc='upper left');axes[2].grid(which='both',alpha=.15)
    fig.suptitle('Published eight-row motor model: local periodic branches and cluster breakup',fontsize=15)
    image1=output/'motor-postcritical-branches.png';fig.savefig(image1,dpi=180);plt.close(fig)

    fig=plt.figure(figsize=(14,9),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,1.05])
    ax0=fig.add_subplot(gs[0,0]);ax1=fig.add_subplot(gs[0,1]);ax2=fig.add_subplot(gs[1,:]);ax3=fig.add_subplot(gs[2,:])
    q=physical(escape['states']);t=escape['time'];strobe=np.arange(0,len(t),40)
    for case,color in [(0,COLORS[1]),(2,COLORS[0])]:
        distance=np.linalg.norm(q[:,case,6]-q[:,case,7],axis=1)
        ax0.semilogy(t[strobe],distance[strobe],color=color,label=f"epsilon = {meta['cases'][case]['epsilon']:.0e}")
        first=distance[0];fit=first*np.exp(meta['unstable_exponent']*t)
        mask=(fit<.015)&(t<2200)
        ax0.semilogy(t[mask],fit[mask],'--',color=color,alpha=.65,lw=1)
    ax0.set(xlim=(0,2300),ylim=(1e-10,.1),xlabel='Original normalized time tau',ylabel='Distance between rows 7 and 8',
            title='Leaving the initial periodic saddle')
    ax0.legend();ax0.grid(alpha=.15)
    ax0.text(.03,.96,'Dashed: Floquet growth 0.00869627 per tau',va='top',transform=ax0.transAxes,fontsize=8)
    for case,color in [(0,COLORS[1]),(2,COLORS[0])]:
        m=diagnostic['cases'][case]['orbit_matches']
        ax1.semilogy([v['tau'] for v in m],[max(1e-12,v['distance_to_permuted_saddle']) for v in m],
                     color=color,label=f"epsilon = {meta['cases'][case]['epsilon']:.0e}")
    ax1.set(xlabel='Original normalized time tau',ylabel='Distance after phase and row matching',
            title='Approaching a symmetry-related periodic saddle',ylim=(1e-11,.1))
    ax1.legend();ax1.grid(alpha=.15)
    ax1.text(.03,.03,'Endpoint floor is limited by interpolating the reference orbit',transform=ax1.transAxes,fontsize=8)
    basepairs=[(0,1),(2,3),(4,5),(6,7)]
    distance=np.array([np.linalg.norm(q[:,0,i]-q[:,0,j],axis=1) for i,j in basepairs])
    colors=ax2.imshow(np.log10(np.maximum(distance[:,::4],1e-14)),origin='lower',aspect='auto',
                      extent=(t[0],t[-1],-.5,3.5),vmin=-14,vmax=0,cmap='viridis')
    ax2.set(yticks=np.arange(4),yticklabels=['Rows 1-2','Rows 3-4 (preserved)','Rows 5-6','Rows 7-8'],
            xlabel='Original normalized time tau',title='One connection: the old expanding pair splits; a new pair synchronizes')
    fig.colorbar(colors,ax=ax2,label='log10 row-state distance (zero clipped to -14)',shrink=.8)
    # Endpoint traces display which rows coincide, without a nonlinear projection.
    end=t>t[-1]-3*meta['period']
    for row in range(8):
        ax3.plot(t[end]-t[end][0],q[end,0,row,2],lw=1.4 if row%2==0 else .9,
                 color=COLORS[row//2],ls='-' if row%2==0 else '--',label=f'Row {row+1}')
    ax3.set(xlabel='Time from start of final three-period window',ylabel='Physical displacement X',
            title='Destination: rows 3-4 are the expanding pair, held equal by the invariant restriction')
    ax3.legend(ncol=4,loc='upper right',fontsize=7);ax3.grid(alpha=.15)
    fig.suptitle('Finite-time connection evidence at rho = 3.268; no infinite-time connection theorem',fontsize=15)
    image2=output/'motor-postcritical-saddle-connection.png';fig.savefig(image2,dpi=180);plt.close(fig)
    rec=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_files={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        images={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [image1,image2]},
        matplotlib=matplotlib.__version__,numpy=np.__version__,
        interpretation='Physical traces and numerical branch/connection evidence. Certified first crossing is distinguished from numerical subsequent results.')
    manifest.parent.mkdir(parents=True,exist_ok=True);manifest.write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8')
    print(image1);print(image2)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/assets')
    parser.add_argument('--manifest',type=Path,default=BASE/'figure-manifest.json')
    a=parser.parse_args();render(a.output.resolve(),a.manifest.resolve())
