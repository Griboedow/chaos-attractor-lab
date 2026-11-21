"""Plot aligned excursions, physical geometry and the discrete halfspace observable."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def main(a):
    read=lambda path: json.loads(path.read_text())
    rows=[read(p) for p in sorted(a.output.glob("cohort-*.json"))]
    meta=read(a.output/"episodes.json")
    ep=meta["episodes"]; rows0=[v for v in ep if v["cohort"]==0]
    data=read(a.output/"cohort-00.json"); trace=np.load(a.cache/data["trace"],mmap_mode="r")
    with np.load(a.cloud) as f: ref=f["points"].astype(float)
    n=np.asarray(data["normal"]);c=data["cut"];rms=data["rms"]
    navy="#123451"; grey="#91a3ad"; red="#b33e2d"; teal="#007d77"; pale="#ced9dd"
    plt.rcParams.update({"font.size":9,"axes.labelsize":9,"axes.titlesize":11,"font.family":"DejaVu Sans",
                         "axes.spines.top":False,"axes.spines.right":False,"axes.edgecolor":"#84949c",
                         "xtick.color":navy,"ytick.color":navy,"text.color":navy,"axes.labelcolor":navy})
    fig,axs=plt.subplots(2,2,figsize=(12,8.4),layout="constrained")
    ax=axs[0,0]
    for row in rows:
        y=row["cohort"]
        ax.hlines(y,0,row["observed"]/1e6,color=pale,lw=2)
    for row in ep:
        ax.scatter(row["first"]/1e6,row["cohort"],s=8+8*np.log10(row["total"]+1),color=(red if row["phase_flipped"] else teal),zorder=3,lw=0)
    ax.set(xlabel="Iterations after burn-in (millions)",ylabel="Initial condition",yticks=[0,4,8,12,16,20,24])
    ax.invert_yaxis()
    ax.legend(handles=[Line2D([0],[0],ls="",marker="o",color=x,label=y,ms=4)for x,y in [(red,"Phase reverses"),(teal,"Phase resumes")]],loc="lower right",frameon=False)
    ax.set_title("A   Episodes across independent starts",loc="left",fontweight="bold")
    ax=axs[0,1]
    with np.load(a.cache/data["windows"]) as f: windows=f["xyz"]
    center=ref.mean(0); rel=ref-center;rem=rel-(rel@n)[:,None]*n
    tang=np.linalg.eigh(rem.T@rem)[1][:,-1]
    def plane(z): return (z-center)@tang/rms,(z@n-c)/rms
    ax.scatter(*plane(ref[::11]),s=1.2,c=navy,alpha=.16,rasterized=True,lw=0,label="Reference cloud")
    ax.scatter(*plane(windows[:,16]),s=7,c=red,alpha=.6,rasterized=True,lw=0,label="Same-side pair starts")
    ax.axhline(0,color=grey,ls="--",lw=.8)
    ax.set(xlabel="Leading in-plane coordinate (RMS)",ylabel="Signed plane distance (RMS)")
    ax.legend(loc="lower left",markerscale=1.5,frameon=False)
    ax.set_title("B   Pair starts leave the dense lobes",loc="left",fontweight="bold")
    ax=axs[1,0]; exemplar=min(rows0,key=lambda v:v["total"])
    left=exemplar["first"];right=exemplar["last"];z=trace[left-20:right+22];x=np.arange(len(z))-20
    proj=(z@n-c)/rms;same=(proj[:-1]>0)==(proj[1:]>0)
    ax.plot(x,proj,color=navy,lw=.9,marker="o",ms=2)
    for j in np.flatnonzero(same):ax.plot(x[j:j+2],proj[j:j+2],color=red,lw=2.7,marker="o",ms=3)
    ax.axhline(0,color=grey,ls="--",lw=.8)
    ax.set(xlabel=f"Iteration relative to first pair at {left:,}",ylabel="Signed plane distance (RMS)")
    ax.set_title("C   One observed phase reversal",loc="left",fontweight="bold")
    ax=axs[1,1];lags=np.asarray(meta["time_offsets"])
    for key,color,label in [("primary",red,"Original start"),("independent",teal,"Other starts")]:
        a0=np.asarray(meta[key]["median_nearest_rms_on"]);a1=np.asarray(meta[key]["median_nearest_rms_off"])
        mask=lags<0;ax.plot(lags[mask],a0[mask],".-",c=color,ms=4,lw=1.3,label=label)
        mask=lags>0;ax.plot(lags[mask],a1[mask],".-",c=color,ms=4,lw=1.3)
    ax.axhline(meta["cloud_distance_threshold_rms"],color=grey,ls=":",lw=1.2,label="Reference 99.9% cutoff")
    ax.set_xscale("symlog",linthresh=4);ax.set_yscale("log");ax.legend(frameon=False,loc="upper right",fontsize=8)
    ax.set(xticks=[-1024,-64,-4,4,64,1024],xticklabels=["−1024","−64","−4","+4","+64","+1024"],xlabel="Steps before first pair / after last pair",ylabel="Median distance to reference (RMS)")
    ax.set_title("D   The trajectory returns to dense lobes",loc="left",fontweight="bold")
    fig.suptitle("Case 2298665   •   Intermittent excursions from two alternating lobes",fontsize=16,fontweight="bold",x=.035,ha="left")
    a.figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.figure,dpi=155,facecolor="white");plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("cloud","cache","output","figure"):
        p.add_argument("--"+k,required=True,type=Path)
    main(p.parse_args())
