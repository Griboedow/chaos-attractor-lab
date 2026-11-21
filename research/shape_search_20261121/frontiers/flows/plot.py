"""Return-section geometry, resolved map, and attempts to reach the stable cycles."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from native_section import RAW,RESULT

def main():
    fig,axes=plt.subplots(1,3,figsize=(12,3.65),constrained_layout=True)
    e=np.load(RAW/'s000-h0.02-t90000.npz')['events']
    s=np.load(RAW/'stable-orbit-section.npy')
    ax=axes[0];ax.scatter(e[:,1],e[:,2],s=.8,color='#067489',label='Observed chaotic bands')
    ax.scatter(s[:,1],s[:,2],s=12,facecolors='none',edgecolors='#ad4627',lw=.8,label='Seven cycles')
    ax.set(xlabel='Native x',ylabel='Native z',title='Positive physical section')
    ax.legend(frameon=False,fontsize=7)
    a=np.load(RAW/'reduced-small.npz');z=a['z'][:,0];y=a['returns'][:,:,2]
    ax=axes[1];ax.fill_between(z,y.min(1),y.max(1),color='#067489',alpha=.45);ax.plot(z,y[:,8],color='#067489',lw=1.4)
    ax.plot(z,z,ls=':',color='.45');ax.set(xlabel='Initial z',ylabel='Returned z',title='One symmetry-reduced return')
    data=json.loads((RESULT/'fates.json').read_text(encoding='utf-8'))
    ax=axes[2]
    for h,c,m in [(.02,'#ad4627','o'),(.01,'#067489','x')]:
        g=[q for q in data['transverse_stress'] if q['h']==h and q['radius']>0]
        ax.semilogx([q['radius'] for q in g],[q['cycle']/q['total'] for q in g],marker=m,color=c,lw=1.2,label=f'RK4 step {h}')
    ax.set(xlabel='Section offset in x',ylabel='Fraction reaching a stable cycle',ylim=(-.04,1.04),title='Transverse escape attempts')
    ax.legend(frameon=False,fontsize=7);ax.grid(alpha=.18)
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8)
    fig.savefig(RESULT/'return-diagnostic.png',dpi=150)
    plt.close(fig)

if __name__=='__main__':main()
