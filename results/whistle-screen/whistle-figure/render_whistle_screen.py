"""Scientific summary of the negative exploratory whistle screen."""
import argparse,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    base=args.external;cycle=np.load(base/'whistle-source-cycle/cycle.npz');tail=np.load(base/'whistle-returns/case-0426.npz')
    scout=np.load(base/'whistle-pilot/complete-or-budgeted.npz');good=scout['status']==0
    bg='#071320';fg='#e5edf5';muted='#a4b4c5';teal='#43d7ce';amber='#ffba72'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,
        'ytick.color':muted,'axes.edgecolor':'#456078','axes.facecolor':bg})
    fig=plt.figure(figsize=(14,9),facecolor=bg)
    histories=[cycle['states'],tail['states'][tail['time']>tail['time'][-1]-6*5.66069]]
    for column,(x,title,color) in enumerate(zip(histories,['Source B: one attracting cycle','Case426: budget endpoint, six final turns'],[teal,amber])):
        points=np.column_stack([x[:,4],x[:,2],x[:,0]-x[:,1]])
        ax=fig.add_subplot(2,2,column+1,projection='3d',facecolor=bg)
        ax.plot(*points.T,color=color,lw=1.0,alpha=.9)
        ax.set_xlabel('vortex displacement q',labelpad=2);ax.set_ylabel('pressure / supply',labelpad=2);ax.set_zlabel('(Q1 - Q2) / jet flow',labelpad=2)
        ax.set_title(title,fontsize=12,pad=12);ax.view_init(22,-50);ax.grid(False)
        span=np.ptp(points,axis=0);ax.set_box_aspect(np.maximum(span,.08))
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.pane.set_edgecolor('#24374a');axis.line.set_color('#456078')
        ax.tick_params(labelsize=8,pad=0)
    ax=fig.add_subplot(223)
    period=float(cycle['period']);t=cycle['time']/period
    ax.plot(t,cycle['states'][:,2]*cycle['scale'][2],color=teal,lw=2,label='source cycle')
    ax.axhline(93/2,color=muted,lw=1,ls='--',label='half supply pressure')
    ax.set_xlabel('time / period');ax.set_ylabel('cavity pressure [Pa]');ax.set_title('Literal equations: mean44.951Pa; period398.848Hz',fontsize=11,pad=12)
    ax.set_xlim(0,1);ax.legend(frameon=False,labelcolor=fg,fontsize=9);ax.grid(alpha=.12)
    ax=fig.add_subplot(224)
    ids=np.flatnonzero(good);score=scout['lle_windows'][good].min(1);order=np.argsort(score)
    ax.scatter(np.arange(len(ids)),score[order],s=7,color=teal,alpha=.7)
    ax.axhline(.0005,color=amber,ls='--',lw=1.4,label='weak retention threshold')
    ax.axhline(0,color=muted,lw=.8);ax.set_xlabel('completed lanes, sorted by minimum window rate');ax.set_ylabel('minimum of four leading-rate windows')
    ax.set_title('510 completed lanes: none passes weak retention',fontsize=11,pad=12);ax.grid(alpha=.12);ax.legend(frameon=False,labelcolor=fg,fontsize=9)
    fig.text(.05,.963,'LONG-CAVITY WHISTLE | NO RETAINED CHAOS CANDIDATE',fontsize=19)
    fig.text(.05,.924,'Literal six-state empirical model. Independent implicit integration closes all four budget endpoints\nand the two strongest weak completed cases to numerical period-one returns.',fontsize=11,color=muted)
    fig.text(.05,.025,'514 exploratory lanes; four original GPU intervals remain incomplete. Physical coefficient sensitivity, not experimental calibration.',fontsize=10,color=muted)
    fig.subplots_adjust(top=.84,bottom=.105,left=.075,right=.96,wspace=.28,hspace=.42)
    fig.savefig(out/'whistle-screen.png',dpi=170,facecolor=bg);plt.close(fig)
    files=['whistle-source-cycle/cycle.npz','whistle-returns/case-0426.npz','whistle-pilot/complete-or-budgeted.npz']
    (out/'render.json').write_text(json.dumps(dict(inputs=[dict(file=f,sha256=hashlib.sha256((base/f).read_bytes()).hexdigest()) for f in files],
        projection='q,pc/Ps,(Q1-Q2)/(A Uj),literal linear coordinates; aspect proportional to coordinate spans.',
        selection='Source one-period shooting control;case426 is one of four unresolved scout intervals,continued independently;all510complete scout scores.'),indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
