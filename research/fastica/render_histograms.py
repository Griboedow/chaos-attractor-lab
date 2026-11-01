"""Count every consecutive state in a declared native-coordinate projection."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from paths import output_root,ROOT
sys.path.insert(0,str(ROOT/'tools'))
from render_portraits import camera
def make(source,output):
    output.mkdir(parents=True,exist_ok=True);rows=[]
    for case in [204,351]:
        path=source/f'case-{case}.npz'
        with np.load(path) as data:x=data['cloud']
        assert x.shape==(1000000,4) and np.isfinite(x).all()
        z=x[:,[0,1,3]];center=(z.min(0)+z.max(0))/2;C=camera();xy=(z-center)@C.T;half=float(abs(xy).max()*1.06)
        h=np.histogram2d(xy[:,0],xy[:,1],bins=760,range=[[-half,half]]*2)[0].T[::-1].astype(np.uint32)
        assert int(h.sum())==len(x)
        target=output/f'case-{case}-histogram.npz';np.savez_compressed(target,histogram=h)
        rows.append({'case':case,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'state_count':len(x),'columns_zero_based':[0,1,3],'center':center.tolist(),'camera':C.tolist(),'half_width':half,'clipped':0,'histogram_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'sampling':'Every consecutive state after 10000 discarded steps; no sign folding, symmetry replication, interpolation or decimation.'})
    (output/'render-inputs.json').write_text(json.dumps({'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows},indent=2)+'\n',encoding='utf-8',newline='\n')
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=output_root()/'histograms');a=p.parse_args();make(output_root()/'portraits-01',a.output)
