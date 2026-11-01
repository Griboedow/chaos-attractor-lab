"""Verify archived hashes, full-spectrum volume checks and local report links."""
import argparse,ast,hashlib,json,re
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    a=p.parse_args();root=a.root.resolve();packet=root/'results/adabelief-quadratic'
    manifest=json.loads((packet/'manifest.json').read_text(encoding='utf-8'))
    for entry in manifest['files']:
        target=root/entry['path']
        assert target.is_file(),str(target)
        assert target.stat().st_size==entry['bytes'],str(target)
        assert hashlib.sha256(target.read_bytes()).hexdigest()==entry['sha256'],str(target)
    checked=0;worst=0.
    for path in packet.rglob('run-*.npz'):
        data=np.load(path)
        assert all(np.isfinite(data[k]).all() for k in data.files),str(path)
        expected=np.log(np.prod(data['parameters'][:2]))
        error=float(np.max(abs(data['windows'].sum(axis=-1)-expected)))
        assert error<1e-8,(str(path),error)
        assert np.allclose(data['windows'].mean(0),data['spectrum'],rtol=0,atol=1e-14)
        assert np.all(data['minimum']<=data['maximum'])
        checked+=1;worst=max(worst,error)
    original=np.load(packet/'pilot.npz');replay=np.load(packet/'screen-replay/pilot.npz')
    assert all(np.array_equal(original[k],replay[k]) for k in original.files)
    source=(packet/'source/upstream-AdaBelief.py').read_bytes().removesuffix(b'\n')
    blob=hashlib.sha1(b'blob '+str(len(source)).encode()+b'\0'+source).hexdigest()
    assert blob=='cf31bd8ac0fbc7b22e87fa71512452103c1d097c'
    link_count=0
    for report in [root/'docs/adabelief-quadratic.md',packet/'README.md']:
        for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            assert (report.parent/target.split('#')[0]).exists(),(str(report),target)
            link_count+=1
    for path in (root/'research').glob('*adabelief*.py'):
        ast.parse(path.read_text(encoding='utf-8'))
    print(json.dumps(dict(hashes=len(manifest['files']),full_spectrum_runs=checked,
        maximum_window_volume_error=worst,local_links=link_count,
        upstream_original_blob=blob,screen_arrays_identical=True),indent=2))


if __name__=='__main__':main()
