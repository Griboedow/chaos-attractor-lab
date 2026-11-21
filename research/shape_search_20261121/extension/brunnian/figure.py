"""Draw coexisting solenoids in two Brunnian configurations from their native XYZ states."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

from sample import centers

PALETTE=["#007D80","#BF542E","#55448C"]


def draw(cache,designs,out):
    selected=[designs[0],designs[6]]
    fig=plt.figure(figsize=(12.6,7.0),facecolor="#fbfaf7")
    layout=fig.add_gridspec(2,2,left=.055,right=.97,top=.85,bottom=.10,hspace=.29,height_ratios=[2.3,1])
    fig.text(.055,.95,"Three chaotic attractors linked only as a triple",size=23,weight="bold",va="center",color="#212c34")
    fig.text(.056,.895,"Each color follows one trajectory in a different trapping tube. Both views use the same native spatial scale.",size=11,color="#465661")
    view=np.array([[1.,0.,0.],[0.,.72,.693974],[0.,-.693974,.72]])
    for j,row in enumerate(selected):
        ax=fig.add_subplot(layout[0,j])
        with np.load(cache/f"case-{row['case']}.npz") as f:
            points=f["points"]@view.T
            labels=f["component"]
        order=np.argsort(points[:,2])
        from matplotlib.colors import to_rgb
        colors=np.asarray([to_rgb(c) for c in PALETTE])
        ax.scatter(points[order,0],points[order,1],s=.38,c=colors[labels[order]],linewidths=0,alpha=.85,rasterized=True)
        ax.set(xlim=(-4.42,4.42),ylim=(-3.3,3.3),aspect="equal")
        ax.axis("off")
        ax.set_title(f"Case {row['case']}    |μ₁₂₃| = {abs(row['mu_123'])}",loc="left",fontsize=13,fontweight="bold",color="#24343c")
        ax=fig.add_subplot(layout[1,j],facecolor="#fbfaf7")
        theta=np.linspace(0,2*np.pi,3200,endpoint=True)
        segs=[];weights=[];inks=[]
        for strand in range(3):
            u,v=centers(theta,row["word"],strand)
            xy=np.column_stack((theta/(2*np.pi),u))
            segs.extend(np.stack((xy[:-1],xy[1:]),axis=1))
            weights.extend((v[1:]+v[:-1])/2)
            inks.extend([PALETTE[strand]]*(len(theta)-1))
        order=np.argsort(weights)
        ax.add_collection(LineCollection(np.asarray(segs)[order],colors=np.asarray(inks)[order],linewidths=2.1))
        ax.set(xlim=(0,1),ylim=(-1.2,1.2),xticks=[0,.25,.5,.75,1],yticks=[])
        ax.set_title("Trapping-tube cores, cut along a meridian",loc="left",fontsize=10,color="#465661")
        for spine in ("top","left","right"): ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#b7bfc3")
        ax.tick_params(labelsize=8,colors="#55646d")
        ax.set_xlabel("One turn around the braid axis",fontsize=9,color="#465661")
    fig.savefig(out,dpi=155,facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("designs",type=Path);ap.add_argument("--cache",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    draw(args.cache,json.loads(args.designs.read_text())["cases"],args.out)
