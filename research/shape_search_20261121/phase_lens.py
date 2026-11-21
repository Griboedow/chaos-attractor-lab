"""Draw measured nested return phases with a fixed equal-unit camera within each case."""
import argparse
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from inspect_views import VIEWS
from render_portraits import camera

COLORS=LinearSegmentedColormap.from_list("phases",["#007888","#478f51","#d59a35","#c04b4f","#6a56a2","#126e93"])


def score(z,stride,period,c):
    q=z[::stride]@c.T
    labels=np.arange(len(q))%period
    center=np.array([q[labels==k].mean(0) for k in range(period)])
    radius=np.array([np.quantile(np.linalg.norm(q[labels==k]-center[k],axis=1),.9) for k in range(period)])
    ds=np.linalg.norm(center[:,None,:]-center[None,:,:],axis=2)
    np.fill_diagonal(ds,np.inf)
    return np.percentile(ds.min(1)/(radius+1e-12),10)


def draw(cache,out):
    fig,axes=plt.subplots(2,3,figsize=(12.5,7.2))
    fig.set_facecolor("white")
    fig.subplots_adjust(top=.835,bottom=.1,left=.035,right=.985,hspace=.43,wspace=.14)
    fig.text(.05,.94,"Nested phases in two twist maps",size=23,weight="bold",color="#182f39")
    fig.text(.052,.889,"Each zoom follows one return class. Bar lengths are in the original map coordinates.",
             size=11,color="#46616b")
    designs=[(20002293,[(1,34),(34,5),(170,1)]),(20800299,[(1,13),(13,6),(78,7)])]
    for row,(ident,levels) in enumerate(designs):
        with np.load(cache/f"maps/clouds/{ident}.npz") as data:
            z=data["points"].astype(float)
        cams=[camera(e,a) for e in (0,25,50,70) for a in range(0,180,20)]
        c=max(cams,key=lambda q:min(score(z,stride,p,q) for stride,p in levels if p>1))
        for ax,(stride,n) in zip(axes[row],levels):
            q=z[::stride]@c.T
            mid=(q.min(0)+q.max(0))/2
            q-=mid
            half=np.abs(q).max()*1.11
            labels=np.arange(len(q))%n
            ix=np.arange(len(q))
            np.random.default_rng(171).shuffle(ix)
            ax.scatter(q[ix,0],q[ix,1],c=labels[ix]/max(1,n-1),cmap=COLORS,
                       s=2.8 if n>1 else 1.8,alpha=.85,edgecolors="none",vmin=0,vmax=1)
            ax.set(xlim=(-half,half),ylim=(-half,half),aspect="equal")
            desired=.62*half
            base=10.**np.floor(np.log10(desired))
            bar=max(b for b in (.1*base,.2*base,.5*base,base,2*base,5*base) if b<=desired)
            x,y=-.90*half,-.90*half
            ax.plot([x,x+bar],[y,y],color="#293e49",lw=2.3)
            ax.text(x,y+.045*half,f"{bar:g}",fontsize=8,color="#384a54")
            ax.set_axis_off()
            title=f"Case {ident}" if stride==1 else f"Every {stride}th iterate"
            ax.set_title(title,fontsize=12,loc="left",color="#152f3b",weight="bold")
            desc=(f"{n} return classes" if n>1 else "One of the 170 classes")
            ax.text(0,-.07,desc,fontsize=10,ha="left",transform=ax.transAxes,color="#4a616b")
    fig.savefig(out,dpi=150,facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    draw(args.cache,args.out)
