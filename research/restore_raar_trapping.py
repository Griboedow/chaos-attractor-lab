"""Restore immutable packed-grid inputs without committing large NPY files."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import numpy as np
from audit_raar_trapping_batch import file_hash


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    protocol = json.loads((args.packet/'candidate/protocol.json').read_text())
    shutil.copy2(args.packet/'candidate/protocol.json', args.out/'protocol.json')
    with gzip.open(args.packet/'candidate/mask.npy.gz','rb') as src, (args.out/'mask.npy').open('wb') as dst:
        shutil.copyfileobj(src,dst)
    assert file_hash(args.out/'mask.npy') == protocol['mask_sha256']
    mask=np.load(args.out/'mask.npy',mmap_mode='r')
    pieces=[]
    for begin in range(0,len(mask),1_000_000):
        bits=np.unpackbits(np.asarray(mask[begin:begin+1_000_000]),bitorder='little')
        pieces.append((np.flatnonzero(bits)+8*begin).astype(np.int64))
    keys=np.concatenate(pieces)
    assert len(keys)==protocol['initial_cells'] and np.all(np.diff(keys)>0)
    np.save(args.out/'keys.npy',keys)
    assert file_hash(args.out/'keys.npy') == protocol['keys_sha256']
    print(json.dumps({'status':'restored_exact_bytes','source_cells':len(keys),'mask_sha256':protocol['mask_sha256'],'keys_sha256':protocol['keys_sha256']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,default=Path('results/raar-trapping'));p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
