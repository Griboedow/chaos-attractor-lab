"""Make a labeled figure from the frozen finite return sample."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
if args.output.exists():
    raise FileExistsError('Choose a new figure folder')
args.output.mkdir(parents=True)
data=np.load(args.input/'sample.npz');times=data['section_time'];intervals=np.diff(times)
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
axes[0].plot(times[1:],intervals,'.',ms=4)
axes[0].set(xlabel='Time from the archived spectrum endpoint',ylabel='Time between section crossings',title='All 124 positive displacement crossings')
axes[1].scatter(intervals[:-1],intervals[1:],s=12,alpha=.7)
axes[1].set(xlabel='Return interval n',ylabel='Return interval n+1',title='Consecutive return intervals')
for ax in axes:
    ax.grid(alpha=.15)
fig.suptitle('Five-cluster candidate: a separate 1,000-tau return check\nFinite evidence; no proof of aperiodicity or attraction',fontsize=12)
fig.savefig(args.output/'returns.png',dpi=160);plt.close(fig)
(args.output/'source.py').write_bytes(Path(__file__).read_bytes())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(args.output/'analysis.json').write_text(json.dumps(dict(status='complete',data_unchanged=True,
    sample_sha256=sha(args.input/'sample.npz'),source_sha256=sha(Path(__file__)),image_sha256=sha(args.output/'returns.png'),
    change='Label spacing and readability only; every recorded section crossing retained.'),indent=2)+'\n',encoding='utf-8')
