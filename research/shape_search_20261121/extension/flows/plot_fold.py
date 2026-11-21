"""Plot the observed fold and multiplier refinement as a small diagnostic figure."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap=argparse.ArgumentParser();ap.add_argument('packet',type=Path);a=ap.parse_args()
    rows=[r for r in csv.DictReader((a.packet/'continuation.csv').open()) if r['direction']=='asym-down']
    d=json.loads((a.packet/'bifurcations.json').read_text());ev=list(d.values())[-1]
    fold=ev['fold']['e'];width=(ev['flip']['e']-fold)*1e6
    e=(np.array([float(r['e']) for r in rows])-fold)*1e6;y=np.array([float(r['period']) for r in rows]);mu=np.array([float(r['multiplier']) for r in rows])
    m=e<36
    fig,ax=plt.subplots(1,2,figsize=(8.4,3.05),layout='constrained')
    for face,v,label in zip(ax,[(y-ev['fold']['period'])*1e6,mu],['Period relative to the fold (×10⁻⁶)','Transverse multiplier']):
        face.plot(e[m],v[m],'.-',color='#42445c',lw=1.2,ms=5)
        face.axvspan(0,width,color='#62b99d',alpha=.2)
        face.set_xlim(-1,35);face.set_xlabel('Axial coefficient above fold (×10⁻⁶)');face.set_ylabel(label)
        face.grid(alpha=.22)
    ax[1].axhline(1,color='#64748b',ls='--',lw=.7);ax[1].axhline(-1,color='#64748b',ls='--',lw=.7)
    ax[0].set_title('Two asymmetric branches',loc='left',fontsize=12);ax[1].set_title('Stable window in green',loc='left',fontsize=12)
    fig.savefig(a.packet/'saddle-node-diagnostic.png',dpi=150,metadata={'Software':'Matplotlib'})
    plt.close(fig)


if __name__=='__main__':
    main()
