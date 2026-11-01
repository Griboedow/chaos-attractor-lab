"""Freeze the negative whistle screen, including failed and corrected attempts."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    folders=['whistle-source','whistle-source-cycle','whistle-gpu-controls','whistle-gpu-controls02','whistle-pilot','whistle-returns','whistle-figure','whistle-figure02']
    returns=json.loads((args.external/'whistle-returns/analysis.json').read_text());assert returns['completed']==returns['declared']==6
    inventory=[]
    for folder in folders:
        source=args.external/folder;target=out/folder;target.mkdir()
        for p in sorted(source.iterdir()):
            if not p.is_file() or p.suffix not in ['.json','.npz','.py','.cl','.png']:continue
            compact=p.name=='source-history.npz'
            inventory.append(dict(folder=folder,file=p.name,bytes=p.stat().st_size,sha256=sha(p),included_in_full=not compact))
            if compact:
                d=np.load(p);mask=(d['time']<=1000)|(d['time']>=d['time'][-1]-1000)
                a={k:d[k] for k in d.files if k not in ['time','states']};a.update(time=d['time'][mask],states=d['states'][mask],
                    original_samples=len(d['time']),original_start=d['time'][0],original_end=d['time'][-1])
                np.savez_compressed(target/'source-history-compact.npz',**a)
            else:shutil.copyfile(p,target/p.name)
    write(out/'local-inventory.json',inventory)
    write(out/'scope-and-corrections.json',dict(
        source='Literal printed equations;source pressure sign retained. Declared ambient constants and missing initial chamber pressure are recorded,not fitted.',
        gpu_control='First validation attempt stopped on abs(list) before the first integration;failure and executed source are retained. Corrected controls02 passes.',
        source_discrepancy='Computed source-cycle mean pressure44.951Pa differs from half supply46.5Pa mentioned in source prose. The literal equation is not silently corrected.',
        negative_result='510of514lanes completed;none passes weak threshold. Four step-budget exits remain original incomplete intervals even though their independent endpoints close to numerical period-one returns.',
        extraction='Full one-cycle shooting data,all514scout starts/endpoints,all6independent return histories and all crossing events retained. Original source30000tau history: first1000 and final1000tau at original sampling,full local SHA recorded.',
        final_figure='whistle-figure02: replaces crowded3Dpanels with fixed2Dviews and corrects a frequency label; originalfigure retained as superseded.',
        novelty='No first-chaos,topological novelty,global absence theorem or experimentally calibrated sensitivity claim.'))
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(files),bytes=sum(f['bytes'] for f in files),files=files))
    print('files',len(files),'bytes',sum(f['bytes'] for f in files))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
