"""Locate or explicitly recover optional numerical payloads outside Git."""
from pathlib import Path
from functools import lru_cache
import argparse,hashlib,json,os,re,shutil,subprocess,tempfile
from urllib.parse import quote
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parents[1]

def digest_file(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

@lru_cache(maxsize=1)
def registry():
    return json.loads((ROOT/'provenance/generated-data.json').read_text(encoding='utf-8'))

@lru_cache(maxsize=1)
def records():return {r['path']:r for r in registry()['files']}

def cache_root():
    if os.environ.get('CHAOS_DATA_CACHE'):return Path(os.environ['CHAOS_DATA_CACHE']).expanduser().resolve()
    if os.name=='nt' and Path('D:/').exists():return Path('D:/Codex/chaos-attractor-cache')
    return Path(os.environ.get('XDG_CACHE_HOME',str(Path.home()/'.cache')))/'chaos-attractor-lab'

def cached_path(row,cache=None):
    checksum=row['sha256'];assert re.fullmatch('[0-9a-f]{64}',checksum)
    return (cache or cache_root())/checksum[:2]/(checksum+'.npz')

def repository_path(name):
    path=(ROOT/name).resolve()
    if not path.is_relative_to(ROOT):raise ValueError('Path is outside the repository.')
    return path

def data_path(path):
    """Resolve a local input; never download implicitly."""
    path=Path(path)
    if not path.is_absolute():path=ROOT/path
    path=path.resolve()
    if path.is_file():return path
    if not path.is_relative_to(ROOT):raise FileNotFoundError(path)
    name=path.relative_to(ROOT).as_posix();row=records().get(name)
    if row:
        cached=cached_path(row)
        if cached.is_file():return cached
        raise FileNotFoundError(f'Optional numerical data: {name}. Run python tools/data_cache.py fetch {name}')
    raise FileNotFoundError(path)

def recorded_hash(path):
    """Check retained bytes, or use the explicit receipt for an omitted payload."""
    path=repository_path(path)
    if path.is_file():return digest_file(path)
    return records()[path.relative_to(ROOT).as_posix()]['sha256']

def portrait_inputs():
    catalog=json.loads((ROOT/'provenance/portrait-catalog.json').read_text(encoding='utf-8'))
    names={p['source'] for row in catalog['rows'] for p in row['panels']}
    animation=json.loads((ROOT/'docs/assets/portraits/animation.json').read_text(encoding='utf-8'))
    for row in animation['rows']:
        for panel in row['panels']:
            names.add(panel['source'])
            if panel.get('camera_source'):names.add(panel['camera_source'])
    names.add('results/icosahedral-flow/icosa-long/sample.npz')
    return names

def fetch(row,cache,network_only=False):
    dest=cached_path(row,cache)
    if dest.is_file():
        if digest_file(dest)!=row['sha256']:raise ValueError(f'Corrupt cache file: {dest}')
        return dest,'cache'
    dest.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='fetch-',suffix='.part',dir=dest.parent)
    os.close(fd);temporary=Path(name)
    try:
        available=False
        if not network_only:
            try:
                available=subprocess.run(['git','cat-file','-e',row['git_blob']+'^{blob}'],cwd=ROOT,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode==0
            except FileNotFoundError:pass
        if available:
            with temporary.open('wb') as output:
                subprocess.run(['git','cat-file','blob',row['git_blob']],cwd=ROOT,stdout=output,check=True)
            method='local Git object'
        else:
            commit=registry()['source_commit'];assert re.fullmatch('[0-9a-f]{40}',commit)
            url=f"https://raw.githubusercontent.com/Griboedow/chaos-attractor-lab/{commit}/{quote(row['path'],safe='/')}"
            with urlopen(url,timeout=120) as response,temporary.open('wb') as output:shutil.copyfileobj(response,output,1024*1024)
            method='pinned historical download'
        if temporary.stat().st_size!=row['bytes'] or digest_file(temporary)!=row['sha256']:
            raise ValueError(f"Payload verification failed: {row['path']}")
        temporary.replace(dest)
        return dest,method
    finally:
        if temporary.exists():temporary.unlink()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command',choices=['list','fetch'])
    ap.add_argument('paths',nargs='*',help='Exact repository-relative payload paths.')
    ap.add_argument('--group',action='append',default=[],help='Packet directory below results; repeat for several groups.')
    ap.add_argument('--portraits',action='store_true',help='Inputs required by the current portrait renderers.')
    ap.add_argument('--all',action='store_true',help='Select every omitted payload; several gigabytes.')
    ap.add_argument('--cache',type=Path,default=cache_root())
    ap.add_argument('--materialize',action='store_true',help='Also copy into ignored original paths for legacy scripts.')
    ap.add_argument('--network-only',action='store_true',help='Do not recover from local Git objects.')
    args=ap.parse_args();table=records();names=set(args.paths)
    if args.portraits:names.update(portrait_inputs() & table.keys())
    for group in args.group:
        matches={p for p,r in table.items() if r['group']==group}
        if not matches:ap.error(f'Unknown payload group: {group}')
        names.update(matches)
    if args.all or (args.command=='list' and not names):names=set(table)
    if not names:ap.error('Select paths, --group, --portraits or --all.')
    if names-table.keys():ap.error('Unknown payloads: '+', '.join(sorted(names-table.keys())))
    chosen=[table[p] for p in sorted(names)];cache=args.cache.resolve()
    print(json.dumps({'files':len(chosen),'bytes':sum(r['bytes'] for r in chosen),'cache':str(cache)}),flush=True)
    for row in chosen:
        if args.command=='list':
            print(json.dumps({'path':row['path'],'bytes':row['bytes'],'cached':cached_path(row,cache).is_file()}));continue
        path,method=fetch(row,cache,args.network_only)
        if args.materialize:
            destination=repository_path(row['path'])
            if destination.exists():
                if digest_file(destination)!=row['sha256']:raise ValueError(f'Refusing to overwrite changed file: {destination}')
            else:
                destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,destination)
        print(json.dumps({'path':row['path'],'method':method,'materialized':args.materialize}),flush=True)

if __name__=='__main__':main()
