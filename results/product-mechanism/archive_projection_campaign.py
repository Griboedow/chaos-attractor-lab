"""Archive complete run inventories with large arrays in the optional D: cache."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np


ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    out=ROOT/'results'/args.name;out.mkdir(parents=True,exist_ok=False)
    registry_path=ROOT/'records/generated-data.json';registry=json.loads(registry_path.read_text(encoding='utf-8'))
    omissions=[];inventory=[];ignored=[];runtime_caches=[]
    inputs={'params','points','counts','kinds','active','starts','initial','flags'}
    for folder in args.runs:
        assert (folder/'report.json').is_file(),folder
        for source in sorted(folder.rglob('*')):
            if not source.is_file():continue
            relative=source.relative_to(folder);dest=out/'runs'/folder.name/relative
            dest.parent.mkdir(parents=True,exist_ok=True)
            item=dict(run=folder.name,file=relative.as_posix(),bytes=source.stat().st_size,sha256=sha(source))
            if args.skip_runtime_cache and '__pycache__' in relative.parts:
                runtime_caches.append(item)
                continue
            if source.suffix=='.npz' and source.stat().st_size>=registry['threshold_bytes']:
                checksum=item['sha256'];cache=Path('D:/Codex/chaos-attractor-cache')/checksum[:2]/(checksum+'.npz');cache.parent.mkdir(parents=True,exist_ok=True)
                if not cache.exists():shutil.copyfile(source,cache)
                assert cache.stat().st_size==source.stat().st_size and sha(cache)==checksum
                arrays=[];retained={}
                with np.load(source,allow_pickle=False) as full:
                    for key in full.files:
                        a=full[key];keep=a.nbytes<=65536 or key in inputs
                        arrays.append(dict(key=key,shape=list(a.shape),dtype=str(a.dtype),bytes=a.nbytes,retained=keep))
                        if keep:retained[key]=a
                compact=None
                if retained:
                    cp=ROOT/'results/compact-data'/(checksum+'.npz');np.savez_compressed(cp,**retained)
                    with np.load(cp) as partial,np.load(source) as full:
                        for key in partial.files:assert partial[key].tobytes()==full[key].tobytes()
                    compact=dict(path=cp.relative_to(ROOT).as_posix(),bytes=cp.stat().st_size,sha256=sha(cp),keys=list(retained))
                name=dest.relative_to(ROOT).as_posix()
                row=dict(path=name,bytes=item['bytes'],sha256=checksum,group=args.name,compact=compact,arrays=arrays,recovery='recompute',recompute_command=f'See results/{args.name}/README.md for run {folder.name}. Original large payload {relative.as_posix()} was never committed; reproduce into a new external directory.',original_run_path=str(source))
                registry['files'].append(row);omissions.append(row);ignored.append('/'+name);item['storage']='optional cache';item['compact']=compact
            else:
                shutil.copyfile(source,dest);assert sha(dest)==item['sha256'];item['storage']='retained'
            inventory.append(item)
    (out/'inventory.json').write_text(json.dumps(dict(schema=1,scope='Every research file present in each declared completed run, including intermediate blocks. Explicitly listed runtime bytecode caches are excluded when requested. Optional array bytes remain in the hash-addressed cache and original run directory.',files=inventory,runtime_cache_exclusions=runtime_caches),indent=2)+'\n',encoding='utf-8')
    (out/'omissions.json').write_text(json.dumps(dict(schema=1,files=omissions),indent=2)+'\n',encoding='utf-8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    unique_compact={r['compact']['path']:r['compact'] for r in registry['files'] if r['compact']}
    registry['removed_bytes']=sum(r['bytes'] for r in registry['files']);registry['compact_bytes']=sum(r['bytes'] for r in unique_compact.values())
    registry_path.write_text(json.dumps(registry,indent=2)+'\n',encoding='utf-8')
    with (ROOT/'.gitignore').open('a',encoding='utf-8') as f:f.write('\n# Optional '+args.name+' arrays\n'+'\n'.join(ignored)+'\n')
    print(json.dumps(dict(name=args.name,files=len(inventory),optional=len(omissions),optional_bytes=sum(r['bytes'] for r in omissions))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--skip-runtime-cache',action='store_true',help='Inventory but do not copy __pycache__ files; preserve the external files.')
    main(p.parse_args())
