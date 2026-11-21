"""Deduplicate exact float32 operator/parameter tuples across the three GPU screens."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from screen import parameters
from third_growth import anchors,draw
from evolve_three import mutate


def digest(data):
    rows=np.ascontiguousarray(data.transpose(1,0,2))
    return np.frombuffer(b''.join(hashlib.blake2b(row,digest_size=16).digest() for row in rows),dtype='V16')


def main(args):
    out=[]
    groups=[]
    def add(name,parts):
        v=np.concatenate(parts)
        groups.append(dict(stage=name,rows=len(v),distinct_digests=len(np.unique(v))))
        out.append(v)
    for geom in ('single','linked','smooth'):
        add('base-'+geom,[digest(parameters(b,4096,geom)[0]) for b in range(128)])
        parent=anchors(args.results,geom)
        add('perturbed-'+geom,[digest(draw(parent,geom,b,4096)[0]) for b in range(96)])
        print('completed',geom,flush=True)
    with np.load(args.cache/'three-deep'/'deep-single.npz') as f: ids=f['ids'];data=f['data']
    rows=json.loads((args.results/'three-rate-single.json').read_text())['cases']
    parent=np.array([np.flatnonzero(ids==r['index'])[0] for r in rows if r['min_cpu_third']>.001])
    parts=[]
    for b in range(256):
        raw,_,_=mutate(data,parent,b,4096)
        with np.load(args.cache/'three-evolved'/f'evolved-{b:04d}.npz') as f:
            assert np.array_equal(raw[:,f['index'],:],f['data'])
        parts.append(digest(raw))
    add('evolved-single',parts)
    all_keys=np.concatenate(out)
    unique=len(np.unique(all_keys))
    doc=dict(hash='BLAKE2b-128 over contiguous binary32 parameter storage, including surface selector. The common initial state is fixed.',rows=len(all_keys),distinct_digests=unique,duplicated_rows=len(all_keys)-unique,stages=groups)
    (args.results/'distinct-configurations.json').write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
    print(doc,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--results',type=Path,required=True);p.add_argument('--cache',type=Path,required=True)
    main(p.parse_args())
