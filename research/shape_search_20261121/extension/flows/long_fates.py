"""Distinguish timestep-dependent basin selection from a long chaotic transient."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from branches import rotation
from field import integrate,ring_tangent


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ['data','basins','chaos','out']:ap.add_argument('--'+key,type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    stable=json.loads(a.data.read_text())['representatives']['stable'];starts=json.loads(a.basins.read_text())['broad']
    p=np.array(stable['parameters']);z=np.load(stable['path'])['points'];rms=np.linalg.norm(np.std(z,axis=0));tree=cKDTree(z)
    ref=np.load(a.chaos/'coexisting-chaotic-cloud.npz')['points'];cloud=cKDTree(ref);out=[]
    for idx,h in [(1,.0125),(4,.0125),(0,.00625),(1,.01)]:
        x=np.array(starts[idx]['initial']);rates,pts=integrate(p,x,h,round(1200/h),round(4000/h),15,round(.2/h),ring_tangent)
        d=[];c=[]
        for seg in np.array_split(pts,15):
            tail=seg[::16]
            d.append(float(min(np.quantile(tree.query(tail@rotation(2*np.pi*k/7))[0],.95) for k in range(7))/rms))
            c.append(float(np.quantile(cloud.query(tail)[0],.95)/rms))
        out.append(dict(start_index=idx,step=h,initial=x.tolist(),burn=1200.,window=4000.,rates=rates.tolist(),
            p95_nearest_seven_copy_family_each_window=d,p95_to_chaotic_reference_each_window=c))
        (a.out/'report.json').write_text(json.dumps(out,indent=2)+'\n')
        print(idx,h,rates[0],rates[-1],d[-1],c[-1],flush=True)


if __name__=='__main__':
    main()
