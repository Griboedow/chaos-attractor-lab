"""Show six full-space trajectories approaching the certified periodic orbit."""
import argparse,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation,PillowWriter
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
COLORS=['#d45a36','#2c83a0','#8b63af','#bc8a1f','#328b70','#d05b8e']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'results/motor-cycle-attraction-visual/local-0p03')
    parser.add_argument('--output',type=Path,default=ROOT/'docs/assets/motor-cycle-attraction')
    parser.add_argument('--animate',action='store_true')
    args=parser.parse_args();folder=args.source.resolve();stem=args.output.resolve()
    summary=json.loads((folder/'analysis.json').read_text())
    assert summary['all_six_approach_same_reference_orbit'] and summary['all_integration_controls_passed']
    protocol=json.loads((folder/'protocol.json').read_text());period=protocol['period']
    paths=[folder/f'trajectory-{i:02}.npz' for i in range(6)]
    for p,r in zip(paths,summary['runs']):assert sha(p)==r['trajectory_sha256']
    histories=[np.load(p) for p in paths]
    reference=np.load(folder/'reference.npz')
    coordinates=[16,17,20] # X0, X1, X4; includes a singleton as well as paired clusters.
    orbit=reference['states'][:,coordinates]
    trajectories=[h['strict_states'][:,coordinates] for h in histories]
    joined=np.vstack([orbit]+trajectories)
    low,high=joined.min(axis=0),joined.max(axis=0);margin=(high-low)*.12
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#23354a',
                         'axes.labelcolor':'#344a60','xtick.color':'#526579','ytick.color':'#526579'})
    fig=plt.figure(figsize=(12,6.8),facecolor='#f8fafc')
    fig.text(.5,.952,'Six trajectories, one attracting cycle',ha='center',fontsize=22,fontweight='bold')
    fig.text(.5,.907,'Nearby initial states evolve freely in the full 23-dimensional system',ha='center',fontsize=11,color='#596d80')
    ax=fig.add_axes([.025,.17,.56,.63],projection='3d',facecolor='#f8fafc')
    ax.set(xlabel=r'$X_0$',ylabel=r'$X_1$',zlabel=r'$X_4$',
           xlim=(low[0]-margin[0],high[0]+margin[0]),ylim=(low[1]-margin[1],high[1]+margin[1]),
           zlim=(low[2]-margin[2],high[2]+margin[2]))
    ax.set_box_aspect(high-low);ax.set_proj_type('ortho');ax.view_init(elev=24,azim=-48)
    ax.tick_params(labelsize=8,pad=0)
    for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
        axis.set_pane_color((.93,.95,.97,1))
        axis._axinfo['grid'].update(color=(.70,.76,.82,.4),linewidth=.6)
    colored=[];heads=[]
    for i,(z,color) in enumerate(zip(trajectories,COLORS)):
        colored.append(ax.plot(*z.T,color=color,lw=1.25,alpha=.78)[0])
        ax.scatter(*z[0],s=45,facecolors='white',edgecolors=color,linewidths=1.8,depthshade=False)
        heads.append(ax.plot(*z[-1:,].T,marker='o',ms=4,color=color,linestyle='')[0])
    closed=np.vstack([orbit,orbit[:1]])
    ax.plot(*closed.T,color='#182432',lw=1.5,alpha=.78,ls=(0,(4,2)))
    ax.set_title('Physical-coordinate projection',fontsize=12,pad=7)
    bx=fig.add_axes([.665,.245,.30,.49],facecolor='white')
    distance_lines=[];distance_heads=[]
    for i,(h,color) in enumerate(zip(histories,COLORS)):
        t=h['distance_time']/period;d=h['full_state_orbit_distance']
        distance_lines.append(bx.plot(t,d,color=color,lw=1.5,label=str(i+1))[0])
        distance_heads.append(bx.plot(t[-1:],d[-1:],marker='o',ms=3,color=color,linestyle='')[0])
    bx.set(yscale='log',xlim=(0,6),ylim=(1e-8,.5),xlabel='Elapsed time / cycle period',
           ylabel='Distance to the orbit · full state')
    bx.set_title('Convergence beyond the 3D projection',fontsize=11,pad=12)
    bx.grid(which='major',alpha=.18);bx.spines[['top','right']].set_visible(False)
    bx.legend(title='Start',ncol=6,loc='upper center',bbox_to_anchor=(.5,-.19),frameon=False,
              columnspacing=.9,handlelength=1.1,fontsize=8,title_fontsize=8)
    handles=[Line2D([],[],color='#182432',lw=1.5,ls=(0,(4,2)),label='Attracting cycle'),
             Line2D([],[],marker='o',mfc='white',mec='#526579',ls='',label='Initial states'),
             Line2D([],[],marker='o',color='#526579',ls='',ms=4,label='States after six periods')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.30,.145),ncol=3,frameon=False,fontsize=9,
               columnspacing=1.2,handlelength=1.4)
    fig.text(.5,.050,'They approach the same orbit, not necessarily the same point at the same time.',
             ha='center',fontsize=11,color='#344a60')
    fig.text(.5,.019,'Six specified nearby starts · numerical illustration of local attraction · distances minimize phase on the full-state orbit',
             ha='center',fontsize=8.5,color='#6c7d8f')
    stem.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(stem.with_suffix('.png'),dpi=180,facecolor=fig.get_facecolor())
    outputs=[stem.with_suffix('.png')]
    if args.animate:
        # Real evolution is shown at uniform accelerated time; the camera stays fixed.
        frames=121
        time_text=fig.text(.31,.855,'',ha='center',fontsize=11,color='#344a60')
        handles[-1].set_label('Current states')
        fig.legends[0].remove()
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.30,.145),ncol=3,frameon=False,fontsize=9,
                   columnspacing=1.2,handlelength=1.4)
        def update(frame):
            fraction=frame/(frames-1)
            for h,z,line,head,dl,dh in zip(histories,trajectories,colored,heads,distance_lines,distance_heads):
                k=int(round(fraction*(len(z)-1)));j=int(round(fraction*(len(h['distance_time'])-1)))
                line.set_data_3d(*z[:k+1].T);head.set_data_3d(*z[k:k+1].T)
                dl.set_data(h['distance_time'][:j+1]/period,h['full_state_orbit_distance'][:j+1])
                dh.set_data(h['distance_time'][j:j+1]/period,h['full_state_orbit_distance'][j:j+1])
            time_text.set_text(f'Elapsed time: {fraction*6:.2f} periods')
            return []
        FuncAnimation(fig,update,frames=frames,interval=90).save(stem.with_suffix('.gif'),writer=PillowWriter(fps=12),dpi=100)
        outputs.append(stem.with_suffix('.gif'))
        gif=Image.open(outputs[-1]);assert gif.n_frames==frames and gif.info.get('loop')==0
    plt.close(fig)
    record=dict(source_sha256=sha(__file__),inputs={p.relative_to(ROOT).as_posix():sha(p) for p in paths+[folder/'reference.npz',folder/'analysis.json',folder/'protocol.json']},
                projection=['X0','X1','X4'],full_state_distance='Euclidean norm of all24(c,s,X) physical components with orbit phase minimized.',
                no_rescaling_of_trajectories=True,all_six_starts_displayed=True,
                animation='Fixed camera; actual saved evolution from zero to six periods, uniformly accelerated.' if args.animate else None,
                outputs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in outputs})
    stem.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record['outputs']))

if __name__=='__main__':main()
