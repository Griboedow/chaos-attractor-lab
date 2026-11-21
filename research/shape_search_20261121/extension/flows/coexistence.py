"""Long lifetime, independent spectrum and replay of a second state at the stable-loop coefficient."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from branches import polar_ring
from adaptive import replay
from field import ring_tangent,integrate


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('data','basins','out'):ap.add_argument('--'+k,type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    row=json.loads(a.data.read_text())['representatives']['stable'];p=np.array(row['parameters'])
    starts=json.loads(a.basins.read_text())['broad'];report=dict(e=row['e'],parameters=p.tolist(),long_runs=[])
    write=lambda:(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    for idx,step in [(0,.0125),(1,.01),(4,.008)]:
        x=np.array(starts[idx]['initial'])
        rates,z=integrate(p,x,step,round(1500/step),round(4000/step),15,round(.2/step),ring_tangent)
        report['long_runs'].append(dict(start_index=idx,initial=x.tolist(),step=step,burn=1500.,window=4000.,rates=rates.tolist(),
              rms=float(np.linalg.norm(np.std(z,axis=0))),extrema=np.array([z.min(0),z.max(0)]).tolist()))
        if idx==0:
            path=a.out/'coexisting-chaotic-cloud.npz';np.savez_compressed(path,points=z)
            report['cloud']=dict(cloud_relative=path.resolve().relative_to(a.out.resolve().parents[2]).as_posix(),path=str(path.resolve()),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),shape=list(z.shape),point_key='points',point_dtype='float64',
                parameters=p.tolist(),initial=x.tolist(),step=step,burn=1500.,window=4000.,blocks=15,stride=round(.2/step))
        write();print(idx,min(rates),max(rates),flush=True)
    cfg=SimpleNamespace(seed=2026112123,family='ring',qr=2.,window=240.,burn=800.,blocks=8,rtol=1e-9)
    report['independent_adaptive_spectrum']=replay(180844,p,np.array(starts[0]['initial']),1,cfg)
    write();print(report['independent_adaptive_spectrum']['mean_spectrum'],flush=True)


if __name__=='__main__':
    main()
