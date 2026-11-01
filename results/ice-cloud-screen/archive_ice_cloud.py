"""Freeze all ice-cloud experiments, including the failed difference control."""
import argparse,hashlib,json,shutil
from pathlib import Path


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    folders=['ice-cloud-controls','ice-cloud-algebra','ice-cloud-gpu-controls','ice-cloud-pilot','ice-cloud-returns']
    result=json.loads((args.external/'ice-cloud-returns/analysis.json').read_text())
    assert result['completed']==result['declared']==4
    assert all(r['numerical_period1'] for r in result['runs'])
    for name in folders:
        target=out/name;target.mkdir()
        for p in sorted((args.external/name).iterdir()):
            if p.is_file() and p.suffix in ['.py','.cl','.json','.npz']:shutil.copyfile(p,target/p.name)
    scope=dict(source='Bergner and Spichtinger, arXiv:2507.03475v1 Eq.55--57, final DOI10.1063/5.0297531.',
        source_pdf_sha256='145cd553d7aec8fe6fbaa19fbd7e6db1a66f5acd96a9ea351360c80703360c5d',
        source_scope='Positive-interior full nucleation-mass model, frozen pressure and temperature with updraft source retained.',
        conventions='cp=1005 is declared, not a recovered source constant; final printed sedimentation coefficients retained without inserting an earlier air-drag factor.',
        failure='Original five-point Jacobian assertion failed after all three source integrations. Independent dimensional complex-step control replays the same smooth-side states and passes; original records remain.',
        extension='The coefficient function gained an optional r0 argument after original controls. Default r0=3 and equations are unchanged. Executed source snapshots are frozen in every experiment.',
        selection='All 2048 scout lanes complete; zero primary and two weak. Both starts at both weak tuples receive independent physical-field Radau returns, all numerical period one.',
        archive='All saved source-control histories, scout inputs/endpoints, independent return histories and events are included in full. The original control stores its final ten tau, not every internal solver step.',
        limitation='No global absence theorem, new attractor, first cloud chaos or observational validation claim.')
    (out/'scope-and-corrections.json').write_text(json.dumps(scope,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    files=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(out.rglob('*')) if p.is_file()]
    (out/'manifest.json').write_text(json.dumps(dict(count=len(files),bytes=sum(f['bytes'] for f in files),files=files),indent=2)+'\n')
    print('Archived',len(files),'files,',sum(f['bytes'] for f in files),'bytes')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
