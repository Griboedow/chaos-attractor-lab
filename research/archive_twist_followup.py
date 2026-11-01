"""Freeze source comparisons, finite mechanism probes and ensemble attraction."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    root=Path(__file__).resolve().parents[1];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    mapping=[('periodic-skeleton','twist-periodic-skeleton'),('splitting','twist-splitting'),('approach','twist-approach'),
        ('approach-outliers','twist-approach-outliers'),('unstable-seeds','twist-unstable-seeds'),
        ('prior-embedding','twist-kicked-top-embedding'),('prior-comparison','twist-kicked-top-comparison'),
        ('prior-replay','twist-prior-replay'),('approach-figure','twist-approach-figure02'),('prior-figure','twist-prior-figure02')]
    inventory=[]
    for label,name in mapping:
        source=args.external/name;target=out/label;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix not in ['.py','.json','.npz','.png']:continue
            record=dict(run=label,name=path.name,bytes=path.stat().st_size,sha256=digest(path))
            if label=='splitting' and path.suffix=='.npz':
                data=np.load(path);contents={k:data[k] for k in data.files};count=len(contents['states'])
                indices=np.unique(np.r_[np.arange(0,count,10),count-1])
                for key in ['states','normal','stable']:contents[key]=contents[key][indices]
                contents['original_sample_indices']=indices;contents['original_sample_count']=np.array(count)
                dest=target/(path.stem+'.compact.npz');np.savez_compressed(dest,**contents)
                record.update(mode='States and two tangent directions sampled every10thpoint, endpoint retained. All200000angles and scalar expansion logs retained.',retained=len(indices),samples=count)
            else:
                dest=target/path.name;shutil.copyfile(path,dest);record['mode']='byte-for-byte'
            record.update(archived=dest.relative_to(out).as_posix(),archived_sha256=digest(dest));inventory.append(record)
    skeleton=json.loads((out/'periodic-skeleton/analysis.json').read_text());sphere_errors=[]
    for row in skeleton['runs']:
        r=row['parameters'][0]
        for orbit in row['periodic_orbits']:
            if orbit['least_period']!=1:continue
            x,y,z=orbit['states'][0];sphere_errors.append(abs((x-r/(1+r))**2+y*y+z*z-1/(1+r)**2))
    assert max(sphere_errors)<1e-9
    approach=json.loads((out/'approach/analysis.json').read_text())['runs'];assert [r['period3_count'] for r in approach]==[12314,0]
    for ident in [1709,3254]:
        d=np.load(out/f'approach/case-{ident}.npz');assert len(d['initial'])==32768
        residual=np.array([np.max(np.linalg.norm(d['final24'][:,lag:]-d['final24'][:,:-lag],axis=2),axis=1) for lag in [1,2,3]]).T
        assert np.max(abs(residual-d['return_residual']))<1e-14
    follow=json.loads((out/'approach-outliers/analysis.json').read_text());assert follow['count']==follow['positive_all_windows']==52
    unstable=json.loads((out/'unstable-seeds/analysis.json').read_text())['runs']
    assert all(g['positive_all_windows']==128 and g['period3']==0 for r in unstable for g in r['groups'])
    replay=json.loads((out/'prior-replay/analysis.json').read_text());assert replay['periods']==[1]*33
    assets=[]
    for folder,name in [('twist-approach-figure02','twist-approach.png'),('twist-prior-figure02','twist-prior-comparison.png')]:
        target=root/'docs/assets'/name;assert not target.exists();shutil.copyfile(args.external/folder/name,target);assets.append(target)
    sources=[]
    for name in ['gonchenko-2014.pdf','gonchenko-2016.pdf']:
        path=args.external/'twist-prior-source'/name;sources.append(dict(name=name,bytes=path.stat().st_size,sha256=digest(path),included=False,scope='Local primary-source copy for equation and figure inspection; not redistributed.'))
    for name in ['gonchenko-2014.json','gonchenko-2016.json','unavailable-2012.json','kicked-top-1988-metadata.json']:
        shutil.copyfile(args.external/'twist-prior-source'/name,out/name)
    write(out/'source-inventory.json',sources)
    write(out/'local-inventory.json',inventory)
    write(out/'audit.json',dict(fixed_point_sphere_identity_maximum_residual=max(sphere_errors),ensemble_starts_per_candidate=32768,
        sail_period3_outcomes=12314,medusa_long_outliers=52,outliers_positive_all_windows=52,linear_unstable_seed_trajectories=512,
        unstable_seed_positive_all_windows=512,prior_caption_replay_starts=33,prior_caption_replay_fixed_outcomes=33,
        scope='Finite numerical controls and explicit source-family embedding; no certified chaos, new topology or historical priority.'))
    write(out/'asset-manifest.json',dict(paths_relative_to='repository root',files=[dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in assets]))
    (out/'README.txt').write_text('This is a separate follow-up packet; results/twist-geometry remains unchanged.\nEvery declared ensemble point, its eleven snapshots, and final24states are retained.\nAll52outlier continuations and all512linear unstable seed histories are retained.\nOnly splitting state/direction arrays are thinned by10, with endpoints; all scalar angles and expansion logs remain complete.\nDistinct compact filenames and original/archived hashes are mapped by local-inventory.json.\nThe initial six-panel prior comparison remains frozen, including its unqualified caption. The reviewed four-panel figure and report explicitly record the unreproduced OCR-caption case.\nPublic author full-text OCR was read for1988; its restricted PDF download was not bypassed. SourcePDFs listed in source-inventory.json remain local.\n',encoding='utf8')
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries));print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
