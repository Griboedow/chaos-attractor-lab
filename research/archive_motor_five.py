"""Freeze the five-row search, thinned histories, full events and visual evidence."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np
from PIL import Image


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    root=Path(__file__).resolve().parents[1];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    mapping=[('scout','motor-five-asymmetry'),('checks','motor-five-checks'),('long','motor-five-long'),
        ('refined','motor-five-refined'),('starts','motor-five-starts'),('start-comparison','motor-five-start-comparison'),
        ('symmetry','motor-five-symmetry'),('gallery','motor-five-gallery'),('rotation','motor-five-rotation')]
    inventory=[]
    for label,name in mapping:
        source=args.external/name;target=out/label;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix in ['.json','.py','.cl'] or (label=='gallery' and path.suffix=='.png'):
                shutil.copyfile(path,target/path.name)
            elif path.suffix=='.npz':
                record=dict(run=label,name=path.name,bytes=path.stat().st_size,sha256=digest(path))
                if path.name.startswith('case-') and label in ['checks','long','refined','starts']:
                    data=np.load(path);stride=10 if label=='checks' else 5
                    count=len(data['time']);indices=np.unique(np.r_[np.arange(0,count,stride),count-1])
                    contents={k:data[k] for k in data.files}
                    assert np.isfinite(contents['states']).all() and len(contents['states'])==count
                    for key in ['time','states']:contents[key]=contents[key][indices]
                    contents['original_sample_indices']=indices;contents['original_sample_count']=np.array(count)
                    dest=target/(path.stem+'.compact.npz');np.savez_compressed(dest,**contents)
                    record.update(archived=dest.relative_to(out).as_posix(),mode='Only time/states are thinned; all other arrays and section events are unchanged.',stride=stride,samples=count,retained=len(indices))
                else:
                    dest=target/path.name;shutil.copyfile(path,dest);record.update(archived=dest.relative_to(out).as_posix(),mode='byte-for-byte')
                record['archived_sha256']=digest(dest);inventory.append(record)
    scout=np.load(out/'scout/complete-or-budgeted.npz')
    assert len(scout['status'])==4032 and np.all(scout['status']==0)
    primary=(scout['lle']>.005)&(scout['lle_windows'].min(1)>.001);strict=scout['lle_windows'].min(1)>.01
    assert primary.sum()==751 and strict.sum()==724
    for label,count in [('checks',12),('long',4),('refined',1),('starts',16)]:
        a=json.loads((out/label/'analysis.json').read_text());assert a['declared']==a['completed']==count
        for row in a['runs']:
            original=args.external/dict(mapping)[label]/row['file'];assert digest(original)==row['sha256']
    symmetry=json.loads((out/'symmetry/audit.json').read_text())
    close=[r['permutation'] for r in symmetry['permutations'] if r['within_three_times_identity']]
    assert close==[[0,1,2,3,4],[1,4,3,0,2],[2,3,1,4,0],[3,0,4,2,1],[4,2,0,1,3]]
    # A post-hoc rotating-wave diagnostic in the empirically selected row order.
    original=args.external/'motor-five-long/case-2816.npz';d=np.load(original);x=d['states'];t=d['time']
    order=[0,1,4,2,3];positions=np.c_[x[:,10:],-x[:,10:].sum(1)][:,order]
    angle=2*np.pi*np.arange(5)/5;h=np.sqrt(2/5)*np.array([np.cos(angle),np.sin(angle),np.cos(2*angle),np.sin(2*angle)])
    q=positions@h.T;z=np.c_[q[:,0]+1j*q[:,1],q[:,2]+1j*q[:,3]];phases=np.unwrap(np.angle(z),axis=0)
    combination=phases[:,0]+2*phases[:,1]
    write(out/'phase-diagnostic.json',dict(source_sha256=digest(original),post_hoc=True,row_order=order,contrast_matrix=h.tolist(),
        mean_angular_speeds=((phases[-1]-phases[0])/(t[-1]-t[0])).tolist(),
        complex_mode_amplitude_ranges=np.c_[abs(z).min(0),abs(z).max(0)].tolist(),
        theta1_plus_2theta2_range=[float(combination.min()),float(combination.max())],
        theta1_plus_2theta2_circular_coherence=float(abs(np.exp(1j*combination).mean())),
        scope='Finite sampled phase-locking diagnostic after choosing the cyclic row order. Not an exact synchronization, invariant-torus or rotating-wave theorem.'))
    assets=[]
    for name in ['motor-five-ribbon.png','motor-five-ribbon.gif']:
        target=root/'docs/assets'/name;assert not target.exists();shutil.copyfile(args.external/'motor-five-rotation'/name,target);assets.append(target)
    gif=Image.open(root/'docs/assets/motor-five-ribbon.gif');assert gif.n_frames==72 and gif.size==(1200,660)
    frames=[]
    for i in [0,18,36,54]:gif.seek(i);frames.append(np.array(gif.convert('RGB')))
    differences=[float(np.mean(abs(a.astype(float)-b.astype(float)))) for a,b in zip(frames,frames[1:])];assert min(differences)>.2
    write(out/'audit.json',dict(scout_completed=4032,scout_primary=int(primary.sum()),scout_strict=int(strict.sum()),
        short_full_spectra=12,long_full_spectra=4,refined=1,fresh_starts=16,empirical_close_permutations=close,
        quarter_turn_image_changes=differences,scope='Numerical finite-time evidence and rendering provenance; no certified chaos, topology, basin probability or historical priority.'))
    write(out/'local-inventory.json',inventory)
    write(out/'asset-manifest.json',dict(paths_relative_to='repository root',files=[dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in assets]))
    (out/'README.txt').write_text('Original run JSON files retain the SHA-256 and filename of the full local histories.\nThe compact NPZ files have distinct names; local-inventory.json maps both hashes.\nOnly time and state sampling is thinned: short checks by10, long/refined/fresh starts by5, always retaining the endpoint.\nSection events, final tangent bases, accumulated exponents and all other arrays are retained unchanged.\nThe rotation/display.npz contains the entire unthinned displayed800-unit tail.\nOriginal full local files are preserved. No old result packet is modified.\n',encoding='utf8')
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries))
    print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
