"""Map the dynamics sampled on a diagonal of the certified entropy region."""
import argparse
import csv
from pathlib import Path

import numpy as np

from confirm import native_spectrum


def parameters(t):
    a,b=.034+.001*t,.024+.001*t
    return np.array([4.48+.04*t,a,b,1-a-b,.00045+.0001*t])


def period(z):
    for p in range(1,257):
        if np.max(np.abs(z[-512:]-z[-512-p:-p])) < 1e-8:
            return p
    return 0


def run(out,points,window):
    out.mkdir(parents=True,exist_ok=True)
    kept=[]
    with (out/"line.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f)
        writer.writerow(["t","r","a","b","c","s","mean_l1","mean_l2","mean_l3",
                         "min_l1","min_l2","max_l1","period","determinant_error"])
        for t in np.linspace(0,1,points):
            par=parameters(t)
            rates,logs,z=native_spectrum(par,np.array([.59,1.12,.47]),40000,6,window,1024)
            error=np.max(np.abs(rates.sum(axis=1)-logs))
            assert error<1e-6
            per=period(z) if rates[:,0].max()<0 else 0
            writer.writerow([t,*par,*rates.mean(0),rates[:,0].min(),rates[:,1].min(),
                             rates[:,0].max(),per,error])
            kept.append(z[-180:,0])
    np.savez_compressed(out/"tails.npz",abundance=np.asarray(kept))
    print(points,"grid points,",6*window,"measured steps per point")


def plot(line,out):
    import matplotlib.pyplot as plt
    x=np.genfromtxt(line/"line.csv",names=True,delimiter=",")
    with np.load(line/"tails.npz") as f:
        tail=f["abundance"]
    fig,(a,b)=plt.subplots(2,1,figsize=(9,6.4),sharex=True,
                           height_ratios=[1.5,1],layout="constrained")
    for ax in (a,b):
        ax.spines[["top","right"]].set_visible(False)
    a.scatter(np.broadcast_to(x["t"][:,None],tail.shape).ravel(),tail.ravel(),s=.23,
              color="#114c5c",alpha=.35,edgecolors="none",rasterized=True)
    a.set_ylabel("Abundance")
    b.axhline(0,color="#253d48",lw=.7)
    b.plot(x["t"],x["mean_l1"],lw=1,color="#005d73",label="First exponent")
    b.plot(x["t"],x["mean_l2"],lw=1,color="#b25629",label="Second exponent")
    stable=x["period"]>0
    b.scatter(x["t"][stable],x["mean_l1"][stable],s=6,facecolors="none",edgecolors="#243138",
              linewidths=.5,label="Numerical periodic return")
    b.legend(frameon=False,ncols=3,fontsize=9,loc="lower left")
    b.set(xlabel="Position along the certified parameter segment",ylabel="Rate per update")
    a.set_xlim(0,1)
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=170,facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--count",type=int,default=801)
    ap.add_argument("--window",type=int,default=8000)
    ap.add_argument("--plot",type=Path)
    args=ap.parse_args()
    if args.count:
        run(args.out,args.count,args.window)
    if args.plot:
        plot(args.out,args.plot)
