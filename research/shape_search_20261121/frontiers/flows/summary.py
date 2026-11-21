"""Summarize verified positive-orientation crossings, captures and sampled separation."""
import argparse
import hashlib
import json
import re
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from native_section import RAW,RESULT,orbit_section
from ensemble import starts

def metrics(events,ref):
    c=[cKDTree(ref[ref[:,0]==k,1:]).query(events[:,1:])[0] for k in range(7)]
    d=np.array(c).T
    eps=2e-5
    tail_window=128
    tail=np.max(d[-tail_window:],axis=0)
    copy=int(np.argmin(tail))
    capture=None
    if len(d)>=tail_window and tail[copy]<eps:
        failures=np.flatnonzero(d[:,copy]>=eps)
        first=0 if not len(failures) else int(failures[-1])+1
        capture=float(events[first,0])
    accepted=(events[:,0]>300)
    a=np.min(d[accepted],axis=1)
    tailpts=events[-256:,1:]
    repeats=[float(np.max(np.linalg.norm(tailpts[i:]-tailpts[:-i],axis=1))) for i in range(1,33)]
    return dict(events=len(events),capture_time=capture,captured_copy=copy if capture is not None else None,
        captured_tail_worst=float(tail[copy]),tail_closest_period=int(np.argmin(repeats))+1,
        tail_closest_period_worst=min(repeats),elapsed=float(events[-1,0]),
        cycle_distances_quantiles=np.quantile(a,[0,.001,.01,.5,.99,1]).tolist())

def summarize():
    rs=starts()
    ref=orbit_section()
    np.save(RAW/'stable-orbit-section.npy',ref)
    rows=[]
    for f in sorted(RAW.glob('s*-h*-t*.npz')):
        m=re.fullmatch(r's(\d+)-h([.\deE+-]+)-t([.\deE+-]+)\.npz',f.name)
        if not m:continue
        i=int(m[1])
        try:
            with np.load(f) as z:e=z['events']
        except (ValueError,EOFError,OSError):continue
        r=dict(index=i,kind=rs[i]['kind'],source_index=rs[i]['source_index'],step=float(m[2]),duration=float(m[3]),**metrics(e,ref))
        rows.append(r)
    print('orbits section counts',np.unique(ref[:,0],return_counts=True)[1])
    for key in sorted(set((r['kind'],r['step'],r['duration']) for r in rows)):
        sel=[r for r in rows if (r['kind'],r['step'],r['duration'])==key]
        cap=[r for r in sel if r['capture_time'] is not None]
        survivors=[r for r in sel if r['capture_time'] is None]
        print(key,'n',len(sel),'capture',len(cap),'times',[round(r['capture_time']) for r in cap][:18],
          'uncaptured least distance',min([r['cycle_distances_quantiles'][0] for r in survivors],default=None),flush=True)
    (RAW/'ensemble-summary.json').write_text(json.dumps(dict(section_points=ref.tolist(),runs=rows),indent=2)+'\n')
    return rows

if __name__=='__main__':summarize()
