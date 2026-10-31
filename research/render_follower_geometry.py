"""Plot saved physical angle projections without embedding or coordinate warps."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from check_follower_chain import coefficients, contact_audit


def main(args):
    args.output.mkdir(parents=True,exist_ok=False)
    data = np.load(args.source)
    x,p = data['states'],data['parameters']
    n = x.shape[1]
    length = coefficients(p)[3]
    pairs = list(itertools.combinations(range(n),3))[:4]
    fig = plt.figure(figsize=(16,9),facecolor='#081321')
    plt.rcParams.update({'text.color':'#dce7f4','axes.labelcolor':'#aebfd0',
                         'xtick.color':'#8195aa','ytick.color':'#8195aa','font.size':10})
    for panel,ids in enumerate(pairs):
        ax = fig.add_subplot(2,3,panel+1,projection='3d',facecolor='#081321')
        y = x[:,ids]
        ax.plot(*y.T,color='#69dfc6',lw=.23,alpha=.47,rasterized=True)
        for k in range(3):
            getattr(ax,'set_'+'xyz'[k]+'label')(f'angle {ids[k]+1} (rad)',labelpad=3)
            getattr(ax,'set_'+'xyz'[k]+'lim')(y[:,k].min(),y[:,k].max())
        ax.set_box_aspect(np.maximum(np.ptp(y,axis=0),1e-6))
        ax.view_init(22,-48)
        ax.grid(False)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.pane.fill=False
            axis.pane.set_edgecolor('#23344a')
            axis.set_tick_params(labelsize=7,pad=1)
        ax.set_title('Physical angle projection',fontsize=11)
    ax = fig.add_subplot(2,3,5,facecolor='#081321')
    section = data['section_states']
    if len(section):
        ax.scatter(section[:,1],section[:,2],s=5,c='#f4ae68',alpha=.7)
    ax.set_xlabel('angle 2 (rad)');ax.set_ylabel('angle 3 (rad)')
    ax.set_title('Upward angle-1 = 0 crossings')
    ax.spines[['top','right']].set_visible(False)
    ax = fig.add_subplot(2,3,6,facecolor='#081321')
    directions = np.stack((np.cos(x),np.sin(x)),axis=2)*length[None,:,None]
    position = np.concatenate((np.zeros((len(x),1,2)),np.cumsum(directions,axis=1)),axis=1)
    ax.plot(position[:,-1,0],position[:,-1,1],lw=.2,c='#7b8fea',alpha=.4)
    for i in np.linspace(0,min(len(x)-1,5000),12,dtype=int):
        ax.plot(position[i,:,0],position[i,:,1],lw=.7,alpha=.45,c='#f3bf76')
    ax.scatter([0],[0],c='white',s=20)
    ax.set_aspect('equal');ax.set_xlabel('horizontal position');ax.set_ylabel('vertical position')
    ax.set_title('Physical chain and tip path')
    ax.spines[['top','right']].set_visible(False)
    fig.suptitle(f'{n}-link follower-force chain | geometry audit',fontsize=22,x=.045,ha='left',y=.97)
    fig.text(.045,.025,'Saved trajectory only. Projections are not a 3D embedding. Model ignores contact. Novelty unestablished.',fontsize=10,color='#aabed2')
    fig.subplots_adjust(left=.035,right=.97,bottom=.1,top=.89,wspace=.22,hspace=.23)
    fig.savefig(args.output/'geometry.png',dpi=160,facecolor=fig.get_facecolor())
    plt.close(fig)
    record = dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                  samples=len(x),projections=[list(i) for i in pairs],equal_units_per_3d_axis=True,
                  contact=contact_audit(x,length),
                  warning='Visual screening of saved coordinates; no asymptotic or historical novelty claim.')
    (args.output/'figure.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(json.dumps(record['contact']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
