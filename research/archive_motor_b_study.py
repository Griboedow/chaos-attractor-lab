"""Freeze the motor b slice, periodic branch, symmetry and dispersion evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    folders=['motor-four-b-slice','motor-four-b-checks','motor-four-b-analysis','motor-four-b-analysis02','motor-four-b-analysis03',
        'motor-four-b-cycle-branch','motor-four-b-cycle-branch02','motor-four-b-cycle-fold',
        'motor-four-b-cycle-symmetry','motor-four-b-cycle-symmetry02','motor-four-b-saddle-symmetry',
        'motor-four-b-fold-figure','motor-four-b-fold-figure02','motor-four-b-regimes',
        'motor-four-b-symmetry-0280','motor-four-b-symmetry-0296','motor-four-b-symmetry-0298','motor-four-b-symmetry-0300',
        'motor-four-symmetric-long','motor-four-long-symmetry-0296','motor-four-long-symmetry-0298','motor-four-long-symmetry-0300',
        'motor-four-pair-dispersion','motor-four-pair-solver','motor-four-pair-figure']
    long=json.loads((args.external/'motor-four-symmetric-long/analysis.json').read_text())
    assert long['completed']==long['declared']==3
    inventory=[]
    for folder in folders:
        source=args.external/folder;target=out/folder;target.mkdir()
        assert source.exists()
        for path in sorted(source.iterdir()):
            if not path.is_file() or path.suffix not in ['.json','.npz','.py','.cl','.png','.gif']:continue
            compact=False
            if path.suffix=='.npz':
                data=np.load(path)
                compact=('time' in data.files and 'states' in data.files and data['states'].ndim==2 and len(data['time'])>50000)
                inventory.append(dict(folder=folder,file=path.name,bytes=path.stat().st_size,sha256=sha(path),included_in_full=not compact))
            if compact:
                mask=data['time']>=data['time'][-1]-1200
                values={k:data[k] for k in data.files if k not in ['time','states']}
                values.update(time=data['time'][mask],states=data['states'][mask],first_state=data['states'][0],
                    original_samples=len(data['time']),original_start=data['time'][0],original_end=data['time'][-1])
                np.savez_compressed(target/(path.stem+'-compact.npz'),**values)
            else:shutil.copyfile(path,target/path.name)
    write(out/'selection-and-corrections.json',dict(
        short_checks='IDs280296298300302304320 selected to bracket the local chaotic/cycle regimes after inspecting the b slice; generic runner wording refers to an earlier gallery and is not the precise rationale for these IDs.',
        long_checks='IDs296298300 selected after6000tau full spectra and24-permutation audits suggested a larger near-S4 sampled set; all three declared60000tau runs are retained.',
        final_analysis='motor-four-b-analysis03; replaces misleading distinguished_row output label by minimum_mean_row and clarifies display-only relabeling.',
        cycle_branch='First branch attempt failed before integration because Numba did not support np.r_; complete corrected branch02 is the connected numerical continuation.',
        cycle_symmetry='First scalar-minimization-only attempt left one quarter-period error1.4e-7. Final symmetry02 adds analytic-tangent phase refinement; all original results remain.',
        final_fold_figure='motor-four-b-fold-figure02; layout revision only, after original figure inspection.',
        pair_solver='Original independent-controls.json tested one base trajectory per model, despite its source comment saying pair. The separate pair-solver packet tests both members at four phases,all three deltas,two models,two GPU settings and2/20tau.',
        history_extraction='Large2D state histories: all original solver section events,integrated logs,terminal basis and final1200tau at original sampling retained. Full originals hash-indexed locally. All pair-dispersion checkpoints retained.',
        scope='Floating-point dynamical evidence, not an asymptotic chaos, invariant-set symmetry, global bifurcation or historical novelty theorem.'))
    write(out/'local-inventory.json',inventory)
    write(out/'summary.json',dict(long_runs=[dict(id=r['id'],parameters=r['parameters'],spectrum=r['spectrum'],windows=r['windows'],
        trace_residual_per_tau=r['trace_residual_per_tau']) for r in long['runs']],
        fold=json.loads((args.external/'motor-four-b-cycle-fold/analysis.json').read_text()),
        pair_control={k:v for k,v in json.loads((args.external/'motor-four-pair-solver/analysis.json').read_text()).items() if k!='records'}))
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(files),bytes=sum(r['bytes'] for r in files),files=files));print('files',len(files),'bytes',sum(r['bytes'] for r in files))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
