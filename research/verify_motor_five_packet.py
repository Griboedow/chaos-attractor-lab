"""Verify the frozen five-row packet without repeating chaotic integrations."""
import hashlib,json
from pathlib import Path
import numpy as np


def main():
    root=Path(__file__).resolve().parents[1];packet=root/'results/motor-five-geometry'
    manifest=json.loads((packet/'manifest.json').read_text());assets=json.loads((packet/'asset-manifest.json').read_text())
    for base,entries in [(packet,manifest['files']),(root,assets['files'])]:
        for entry in entries:
            path=base/entry['path'];assert path.stat().st_size==entry['bytes']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'],str(path)
    scout=np.load(packet/'scout/complete-or-budgeted.npz');assert len(scout['status'])==4032 and np.all(scout['status']==0)
    for name,count in [('checks',12),('long',4),('refined',1),('starts',16)]:
        analysis=json.loads((packet/name/'analysis.json').read_text());assert analysis['declared']==analysis['completed']==count
        for row in analysis['runs']:
            data=np.load(packet/name/(Path(row['file']).stem+'.compact.npz'))
            assert np.isfinite(data['states']).all() and (np.diff(data['time'])>0).all()
            assert len(data['time'])==len(data['original_sample_indices'])
            assert data['original_sample_indices'][-1]==int(data['original_sample_count'])-1
            if name!='starts':
                spectrum=np.array(row['spectrum']);assert len(spectrum)==14
                assert np.max(abs(spectrum-np.mean([w['spectrum'] for w in row['windows']],axis=0)))<1e-10
                assert row['trace_residual_per_tau']<1e-7
    symmetry=json.loads((packet/'symmetry/audit.json').read_text())
    close={tuple(r['permutation']) for r in symmetry['permutations'] if r['within_three_times_identity']}
    generator=(1,4,3,0,2);powers=set();current=tuple(range(5))
    for _ in range(5):powers.add(current);current=tuple(generator[i] for i in current)
    assert current==tuple(range(5)) and powers==close and len(close)==5
    display=np.load(packet/'rotation/display.npz');meta=json.loads((packet/'rotation/render.json').read_text());h=np.array(meta['contrast_matrix'])
    assert np.max(abs(h@h.T-np.eye(4)))<1e-14 and np.max(abs(h.sum(1)))<1e-14
    x=display['states'];positions=np.c_[x[:,10:],-x[:,10:].sum(1)]
    assert np.max(abs(positions@h.T-display['contrasts']))<1e-12
    print(json.dumps(dict(files=len(manifest['files']),assets=len(assets['files']),scouts=4032,full_spectra=17,fresh_starts=16,checks='passed')))


if __name__=='__main__':main()
