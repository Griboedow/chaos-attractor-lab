"""Verify the slender-rod packet; optionally replay rotation-symmetry controls."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
PACKET=ROOT/'results/hydrogel-slender-screen'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--numerical',action='store_true')
    args=parser.parse_args()
    manifest=read(PACKET/'manifest.json')
    for row in manifest['files']:
        path=PACKET/row['path']
        assert path.stat().st_size==row['bytes'],path
        assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],path
    for row in manifest['dependencies']:
        path=ROOT/row['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],path
    runs={p.stem:read(p) for p in (PACKET/'runs').glob('*.json')}
    assert set(runs)=={f'lambda127-field{x}' for x in (117,130,150)}
    assert sum(r['success'] for r in runs.values())==2
    for name,row in runs.items():
        assert row['segments']==64 and row['dt']==.006
        assert row['parameters']['slenderness']==127 and row['parameters']['eta']==1
        assert row['requested_duration']==96
        assert row['success']==(name!='lambda127-field150')
        assert row['completed_time']==(52.608000000000004 if not row['success'] else 96)
    arrays=read(PACKET/'arrays.json')['files']
    assert {r['run'] for r in arrays}==set(runs)
    for row in arrays:
        assert row['arrays']['q']['shape'][1:]==[64,4]
        assert row['arrays']['k']['shape'][1:]==[64,2]
    clearance={r['run']:r for r in read(PACKET/'clearance.json')}
    for field,count in [(117,0),(130,137),(150,315)]:
        row=clearance[f'lambda127-field{field}']
        assert row['samples']==417 and row['samples_below_threshold']==count
        assert row['minimum_segment_index_gap']==4
        assert row['tube_diameter_threshold']==2/127
    returns={r['run']:r for r in read(PACKET/'returns.json')}
    assert set(returns)==set(runs)
    for row in returns.values():
        for lag in row['return_errors']:
            assert lag['reduced_rms_max']<=lag['ordinary_rms_max']+1e-12
            assert lag['pairs']==len(lag['reduced_rms_values'])
    control=read(PACKET/'rotation-control.json')[0]
    assert control['run']=='rotation-n64'
    assert control['quotient_spread_from_final']['maximum']<.00801
    previous=read(ROOT/'results/hydrogel-spatial-controls/arrays.json')['files']
    previous={r['run']:r for r in previous}['rotation-n64']
    recorded=read(PACKET/'control-input.json')['trajectory']
    assert recorded['sha256']==previous['sha256'] and recorded['bytes']==previous['bytes']
    checks=read(PACKET/'symmetry-check.json')
    assert checks['residual_equivariance_relative_error']<1e-10
    assert checks['recovered_rotation']['rms']<1e-12
    replay='not requested'
    if args.numerical:
        sys.path.insert(0,str(ROOT/'results/hydrogel-spatial-controls/code'))
        sys.path.insert(0,str(PACKET/'code'))
        import symmetry
        replay=symmetry.check()
    print(json.dumps(dict(retained_hashes=len(manifest['files'])+len(manifest['dependencies']),
        runs=3,completed=2,failed=1,numerical=replay,
        scope='Retained packet, earlier control-input identity and optional symmetry algebra. No long integration replay, continuum or chaos certificate.')))


if __name__=='__main__':
    main()
