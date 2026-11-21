"""Test an observed alternation of two separated point sets over long native orbits."""
import argparse,json
from pathlib import Path
import numpy as np
from numba import njit
from affine_model import step
from confirm import from_record


@njit(cache=True)
def probe(raw,x,burn,length,normal,threshold):
    state=x.copy();lo=np.full(2,np.inf);hi=np.full(2,-np.inf);same=0;prev=False
    for i in range(burn+length):
        state,_=step(state,raw)
        if i>=burn:
            v=np.dot(state,normal);idx=(i-burn)%2
            lo[idx]=min(lo[idx],v);hi[idx]=max(hi[idx],v)
            side=v>threshold
            if i>burn and side==prev:same+=1
            prev=side
    return lo,hi,same


def main(args):
    rec=json.loads(args.record.read_text());raw=from_record(rec)
    with np.load(args.cache_root/rec['cloud_relative']) as f:z=f['points'].astype(float)
    center=z.mean(0);rad=np.sqrt(np.var(z,axis=0).sum())
    normal=z[1::2].mean(0)-z[::2].mean(0);normal/=np.linalg.norm(normal)
    primary=z@normal;cut=float((primary[1::2].min()+primary[::2].max())/2)
    tests=[]
    for j,row in enumerate(rec['cohort']):
        length=args.length if j==0 else args.other
        low,high,same=probe(raw,np.array(row['initial']),100000,length,normal,cut)
        gap=float(max(low[1]-high[0],low[0]-high[1]))
        tests.append(dict(cohort_index=j,observed=length,same_halfspace_successive_pairs=int(same),low=low.tolist(),high=high.tolist(),parity_gap_in_reference_rms=gap/rad))
    doc=dict(case=rec['case'],normal=normal.tolist(),cut=cut,source_mean=center.tolist(),source_rms=rad,source_parity_gap_in_rms=float((primary[1::2].min()-primary[::2].max())/rad),discard=100000,tests=tests)
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
    print(rec['case'],[(t['same_halfspace_successive_pairs'],round(t['parity_gap_in_reference_rms'],4))for t in tests],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--record',type=Path,required=True);p.add_argument('--cache-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--length',type=int,default=20000000);p.add_argument('--other',type=int,default=262144)
    main(p.parse_args())
