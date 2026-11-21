"""Render four measured responses of one cubic oscillator at fixed forcing."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

COLORS=np.array(["#087d80","#bb4f2a","#6952a0"])


def draw(file,report,out):
    assert hashlib.sha256(file.read_bytes()).hexdigest()==report["sha256"]
    with np.load(file) as a:
        p=a["points"];off=a["offsets"];gains=a["gain"];orbits=a["orbit"]
    fig=plt.figure(figsize=(15,5.8),facecolor="#fbfaf7")
    gs=fig.add_gridspec(2,4,top=.72,bottom=.11,left=.035,right=.985,wspace=.13,hspace=.31,height_ratios=[2.3,1])
    fig.text(.038,.937,"A cubic oscillator changes the way its responses are linked",size=23,weight="bold",color="#183342")
    fig.text(.040,.863,"Case 48000001  ·  Same periodic forcing throughout. Each color is a distinct orbit reached in the 15-start sample.",size=11,color="#455962")
    titles=["Three separate cycles","One observed cycle","One three-turn cycle","Three linked cycles"]
    for col,(g,label) in enumerate(zip([.03,.15,.23,.34],titles)):
        ax=fig.add_subplot(gs[0,col]);cut=fig.add_subplot(gs[1,col])
        idx=np.flatnonzero(gains==g)
        unique=sorted(set(orbits[idx]))
        segments=[];depth=[];colors=[]
        for i in idx:
            q=p[off[i]:off[i+1]]
            x,y,z=q.T
            xy=np.column_stack((x,.7*y+.71414284*z))
            segments.extend(np.stack((xy[:-1],xy[1:]),axis=1))
            depth.extend(-.71414284*y[:-1]+.7*z[:-1])
            color=COLORS[unique.index(orbits[i])];colors.extend([color]*(len(x)-1))
            cut.plot(np.linspace(0,1,len(x)),np.hypot(x,y)-3.25,c=color,lw=1.15)
        order=np.argsort(depth)
        ax.add_collection(LineCollection(np.asarray(segments)[order],colors=np.asarray(colors)[order],linewidths=1.45))
        ax.set(xlim=(-6,6),ylim=(-4.5,4.5),aspect="equal");ax.axis("off")
        ax.set_title(f"Gain {g:g}\n{label}",fontsize=12,color="#183342",weight="bold")
        cut.set(xlim=(0,1),ylim=(-2.6,2.6),xticks=[0,.5,1],yticks=[-2,0,2])
        cut.set_facecolor("#fbfaf7");cut.spines[["top","right"]].set_visible(False)
        cut.spines[["bottom","left"]].set_color("#c6c9c6")
        cut.tick_params(labelsize=8,colors="#52636a")
        cut.set_xlabel("One forcing turn",size=9,color="#52636a")
        if col==0:cut.set_ylabel("Re(z)",size=10)
    fig.savefig(out,dpi=145,facecolor=fig.get_facecolor());plt.close(fig)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("trace",type=Path);ap.add_argument("packet",type=Path)
    ap.add_argument("--out",type=Path,required=True);a=ap.parse_args()
    draw(a.trace,json.loads(a.packet.read_text()),a.out)
