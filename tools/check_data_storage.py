"""Check tracked-data size policy and optional compact-record provenance."""
from pathlib import Path
import argparse,json,subprocess
import numpy as np
from data_cache import ROOT,registry,records,data_path,digest_file

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--with-cache',action='store_true',help='Hash original cached payloads and compare every retained array byte for byte.')
    args=ap.parse_args();manifest=registry();table=records()
    tracked=set(subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode('utf-8').split('\0'))-{''}
    assert not (tracked & table.keys()),'Optional payloads must not be tracked. Stage intended deletions before this check.'
    ignored=set(subprocess.check_output(['git','check-ignore','--stdin'],cwd=ROOT,input=('\n'.join(table)+'\n').encode()).decode().splitlines())
    assert ignored==set(table),'Every optional payload must remain ignored.'
    compact={r['compact']['path']:r['compact'] for r in table.values() if r['compact']}
    allowed={r['path'] for r in manifest['retained_large_inputs']}|compact.keys()
    total=0
    for name in tracked:
        path=ROOT/name;assert path.is_file(),name
        size=path.stat().st_size;total+=size
        if name.startswith('results/') and path.suffix=='.npz' and size>=manifest['threshold_bytes']:
            assert name in allowed,f'Unregistered large numerical payload: {name}'
    for name,row in compact.items():
        path=ROOT/name;assert name in tracked and path.stat().st_size==row['bytes']
        assert digest_file(path)==row['sha256'],name
        with np.load(path,allow_pickle=False) as saved:assert saved.files==row['keys'],name
    compared=0;hashed=set()
    for name,row in table.items():
        assert sum(a['retained'] for a in row['arrays'])==len(row['compact']['keys'] if row['compact'] else [])
        if not args.with_cache:continue
        original=data_path(ROOT/name)
        if row['sha256'] not in hashed:
            assert original.stat().st_size==row['bytes'] and digest_file(original)==row['sha256'],name
            hashed.add(row['sha256'])
        if row['compact']:
            with np.load(original,allow_pickle=False) as full,np.load(ROOT/row['compact']['path'],allow_pickle=False) as saved:
                for key in saved.files:
                    a,b=full[key],saved[key]
                    assert a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes(),(name,key)
                    compared+=1
    assert sum(r['bytes'] for r in table.values())==manifest['removed_bytes']
    assert sum(r['bytes'] for r in compact.values())==manifest['compact_bytes']
    print(json.dumps({'omitted_payloads':len(table),'omitted_bytes':manifest['removed_bytes'],'compact_files':len(compact),'compact_bytes':manifest['compact_bytes'],'tracked_bytes':total,'retained_original_arrays_compared':compared,'cache_check':'passed' if args.with_cache else 'not requested'}))

if __name__=='__main__':main()
