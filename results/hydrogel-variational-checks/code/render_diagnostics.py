"""Finite growth and slow drift: a diagnostic, not an attractor portrait."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args()
def read(name):return json.loads((root/name).read_text(encoding='utf-8'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(11,4.2))
colors={117:'#88acc5',124:'#145a86',127:'#377d9b',130:'#182f40'}
for name,label,color in [('spectrum117-n64','Field 117',colors[117]),
                         ('spectrum124-n64','Field 124',colors[124]),
                         ('spectrum-rotation','Rotating control','#75838c')]:
    data=read(name+'/report.json')['observations']
    if name=='spectrum124-n64':data+=read('spectrum124-extension/report.json')['observations']
    data=[r for r in data if r['exponents'] is not None and r['observed_duration']>=24]
    axes[0].plot([r['observed_duration'] for r in data],[r['exponents'][0] for r in data],label=label,color=color,lw=1.4)
axes[0].axhline(0,color='.65',lw=.7)
axes[0].set(xlabel='Observed duration after burn',ylabel='First QR growth rate',title='A. Finite growth (64 segments)')
axes[0].legend(frameon=False,fontsize=9)
rows=read('drift.json')['runs']
for field,names in [(117,['spectrum117-n64']),(124,['spectrum124-n64','spectrum124-extension']),
                    (127,['spectrum127-n64']),(130,['field130-continuation'])]:
    windows=[w for row in rows if row['run'] in names for w in row['windows'] if w['end']-w['start']>=23.99]
    axes[1].plot([(w['start']+w['end'])/2 for w in windows],[w['mean'][2] for w in windows],
        color=colors[field],label=f'Field {field}',marker='o',ms=3,lw=1.1)
axes[1].set(xlabel='Dimensionless time',ylabel='Mean tip height / filament length',title='B. Consecutive 24-unit windows')
axes[1].legend(frameon=False,fontsize=9)
fig.tight_layout(pad=1.5)
fig.savefig(args.output,dpi=170,facecolor='white');plt.close(fig)
