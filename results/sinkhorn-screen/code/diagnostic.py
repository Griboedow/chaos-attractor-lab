"""Draw scalar growth diagnostics, not an attractor portrait."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(12,4.1),layout='constrained')
special={(3,43):('#70a9b7','3 × 3, Case 43'),(4,55):('#54606b','4 × 4, Case 55'),(4,64):('#164e63','4 × 4, Case 64')}
for n in (3,4):
    first=json.loads(Path(f'followup01/n{n}.json').read_text())['runs'];last=json.loads(Path(f'extension01/n{n}.json').read_text())['runs']
    for a,b in zip(first,last):
        blocks=a['blocks']+b['blocks'];x=[r['stop_step']/1e6 for r in blocks];y=[r['rate'] for r in blocks]
        style=special.get((n,a['case']))
        if style:axs[0].plot(x,y,color=style[0],label=style[1],linewidth=1.6,zorder=5)
        else:axs[0].plot(x,y,color='#ced9dd',linewidth=.6,zorder=1)
axs[0].axhline(0,color='#122f42',linewidth=.6);axs[0].set(xlabel='Iteration (millions)',ylabel='Finite-window growth per iteration',title='All 40 flagged trajectories')
axs[0].legend(frameon=False,fontsize=9)
rows=json.loads(Path('spectra02/summary.json').read_text())['runs']
for r in rows:
    style=special.get((r['shape'][0],r['case']))
    if style:axs[1].plot(np.arange(1,len(r['mean_rates'])+1),r['mean_rates'],marker='o',color=style[0],label=style[1])
axs[1].axhline(0,color='#122f42',linewidth=.6);axs[1].set(xlabel='QR column',ylabel='Finite-window growth per iteration',title='200,000-step spectra after step 5,540,000',xticks=range(1,7))
axs[1].legend(frameon=False,fontsize=9)
args.output.parent.mkdir(exist_ok=True,parents=True);fig.savefig(args.output,dpi=180)
