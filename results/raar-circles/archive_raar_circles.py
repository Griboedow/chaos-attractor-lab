"""Freeze the complete declared spatial-circle RAAR study with copy provenance."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
from verify_rk23_controller import read,sha

FOLDERS={
 'scout':'raar-circles-scout01','follow':'raar-circles-follow01','five-million':'raar-circles-long01',
 'gallery':'raar-circles-gallery01','remaining-gallery':'raar-circles-remaining-gallery01','angles':'raar-circles-angles01',
 'precision':'raar-circles-precision01','upstream':'raar-circles-upstream01','cohort':'raar-circles-cohort01',
 'approach':'raar-circles-approach01','render':'raar-circles-render01','evidence-initial':'raar-circles-evidence01','evidence':'raar-circles-evidence02',
 'relaxation':'raar-circles-relaxation01','slice-figures-initial':'raar-circles-slice-figures01','slice-figures':'raar-circles-slice-figures02',
 'geometry-failure01':'raar-circles-geometry01','geometry-failure02':'raar-circles-geometry02','geometry':'raar-circles-geometry03',
 'fixed-proof':'raar-fixed-proof01','algebra':'raar-algebra01','priority':'raar-priority','upstream-source':'raar-upstream-source'}


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);copies=[];omissions=[]
    for dest,name in FOLDERS.items():
        source=args.external/name;targetdir=out/dest;targetdir.mkdir(parents=True)
        assert source.is_dir(),source
        aggregate={'scout':'screen.npz','follow':'follow.npz','five-million':'follow.npz','relaxation':'scan.npz'}.get(dest)
        joined=read(source/aggregate) if aggregate else None
        for p in sorted(source.rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            record=dict(source=str(p),bytes=p.stat().st_size,sha256=sha(p));reason=None
            if dest=='priority' and p.suffix!='.json':reason='Third-party publication or review rendering remains local; links, hashes and inspection scope are public.'
            elif joined is not None and p.name.startswith('block-'):
                d=read(p)
                if 'ids' in d:indices=[list(joined['ids']).index(i) for i in d['ids']]
                else:
                    offset=int(p.stem.split('-')[1]);indices=list(range(offset,offset+len(d['params'])))
                assert all(np.array_equal(joined[k][indices],v,equal_nan=v.dtype.kind in 'fc') for k,v in d.items()),p
                reason='Every array equals its rows in the archived aggregate.'
            elif dest=='render' and p.suffix=='.npy':reason='Raw eight-million-state portrait retained locally; exact histograms, million-state animations, source, starts and hashes are archived.'
            if reason:omissions.append(dict(**record,reason=reason));continue
            target=targetdir/p.relative_to(source);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
            assert sha(target)==record['sha256'];copies.append(dict(**record,path=target.relative_to(out).as_posix()))
    assets=[]
    for folder,name in [('render','raar-circles-geometry.png'),('render','raar-circles-rotation.gif'),('evidence','raar-circles-evidence.png'),('evidence','raar-circles-approach.png'),('slice-figures','raar-relaxation.png'),('slice-figures','raar-library-clouds.png')]:
        source=out/folder/name;target=args.assets/name;assert not target.exists();shutil.copy2(source,target)
        assets.append(dict(path=name,archive_path=source.relative_to(out).as_posix(),sha256=sha(target),bytes=target.stat().st_size))
    (out/'assets.json').write_text(json.dumps(dict(files=assets),indent=2),encoding='utf-8')
    (out/'copy-provenance.json').write_text(json.dumps(dict(copies=copies,omitted_local_files=omissions),indent=2),encoding='utf-8')
    (out/'README.txt').write_text('Spatial-circle RAAR study, 2026-11-01.\nExisting Luke RAAR algorithm on two noncoplanar circles, a toy inconsistent feasibility benchmark.\nComplete 2051-setting screen; all 53 flagged million-step follow-ups; all 53 gallery clouds; six five-million checks; twelve independent 70/90-digit runs; four original ProxToolbox evaluations; 384 cohort runs; 131072 initial-volume trajectories; 771 relaxation-slice settings.\nInterval fixed-point classification, disjointness and constraint-linking certificates are distinct from finite numerical chaos evidence. Five selected maps have no fixed points; R1679 has exactly two, both unstable. Historical uniqueness and asymptotic strange attraction are NOT established.\nEvery declared negative, weak or periodic outcome is retained. Failed geometry/figure audits remain with failure metadata. Third-party papers are linked and hashed, not redistributed; pinned original ProxToolbox source includes its CC-BY 4.0 license and attribution.\n',encoding='utf-8')
    for n in ['archive_raar_circles.py','verify_raar_circles.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.rglob('*')) if p.is_file()]
    assert max(f['bytes'] for f in files)<100000000
    manifest=dict(count=len(files),bytes=sum(f['bytes'] for f in files),files=files)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps({k:manifest[k] for k in ['count','bytes']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--external',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--assets',type=Path,required=True);main(p.parse_args())
